"""Stage 1b calibration: null distribution of the coherence-survey statistics when there is NO periodic component.

For random real stars (v4 series, preprocessing as in coherence_survey.py): fit H_RN (smooth red noise), simulate the star's
delays from its own fitted H_RN (same epochs, errors, bands, band offsets at their prior means), and run exactly the survey
fits (white, H_RN, H_QP per coherence). Gives the null distributions of
  dQP = lnL_QP - lnL_RN          (is a periodic component needed beyond red noise?)
  dCoh = lnL_QP(c = inf) - max_{c <= 1} lnL_QP(c)   (preference for strict coherence)
Output: results/stage1/coherence_null.parquet.
Usage: python scripts/coherence_null.py --n 300 --reps 2 --workers 3
"""
from __future__ import annotations

import argparse
import sys

sys.path.insert(0, __import__("os").path.dirname(__file__))
import chunked  # noqa: E402  (single-threaded BLAS)
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from coherence_survey import CM, LAG  # noqa: E402
from level2 import load, priors_for, to_series  # noqa: E402
from rrlbin.models import fit_qp, fit_rn, fit_white, se_kernel  # noqa: E402
from rrlbin.oc import align_bands, apply_common_mode, design, robust_unwrap  # noqa: E402
from rrlbin.timing import year_labels  # noqa: E402

DAY = 86400.0


def work(job):
    rec, seed, reps = job
    out = []
    try:
        s = apply_common_mode(to_series(rec), CM, year_labels)
        t, y, e, band = s["t"], s["tau"], s["err"], s["band"].astype(int)
        o = np.argsort(t)
        t, y, e, band = t[o], y[o], e[o], band[o]
        pri = priors_for(s["P"], LAG)
        y = align_bands(t, robust_unwrap(t, y, band, s["P"]), e, band, s["P"], pri)
        rn = fit_rn(t, y, e, band, pri, n_amp=10)
        rng = np.random.default_rng(seed)
        X, names = design(t, band, np.average(t, weights=e ** -2))
        off = np.array([pri[b][0] if (pri and b in pri) else 0.0 for b in band])
        C = se_kernel(t, t, rn["A"], rn["ell"]) + np.diag(e ** 2 + rn["s"] ** 2)
        for r in range(reps):
            ys = off + rng.multivariate_normal(np.zeros(t.size), C)
            w = fit_white(t, ys, e, band, pri)["lnl"]
            frn = fit_rn(t, ys, e, band, pri, n_amp=10)
            fqp = fit_qp(t, ys, e, band, pri, n_amp=10)
            row = dict(ogle_id=rec["ogle_id"], rep=r, n=int(t.size), A_rn_true_s=rn["A"] * DAY, ell_true=rn["ell"], s_true_s=rn["s"] * DAY,
                       lnl_white=w, lnl_rn=frn["lnl"], lnl_qp=fqp["lnl"], Pq=fqp["best"]["Pq"] if fqp["best"] else np.nan,
                       c_best=fqp["best"]["c"] if fqp["best"] else np.nan)
            row.update({f"lnl_qp_c{c}": v for c, v in fqp["lnl_by_c"].items()})
            out.append(row)
    except Exception as ex:
        out.append(dict(ogle_id=rec["ogle_id"], err_msg=str(ex)[:200]))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=300)
    ap.add_argument("--reps", type=int, default=2)
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--seed", type=int, default=12)
    ap.add_argument("--out", default="results/stage1/coherence_null")
    a = ap.parse_args()
    d = load("results/real/series_v4")
    d = d[d.ok.astype(bool) & (d.t.apply(len) >= 12)].sample(a.n, random_state=a.seed)
    jobs = [(r, a.seed * 100000 + k, a.reps) for k, r in enumerate(d.to_dict("records"))]
    chunked.run_chunked(work, jobs, a.out, workers=a.workers, chunk=30, flatten=True,
                        manifest=dict(script="coherence_null", ids=[j[0]["ogle_id"] for j in jobs], reps=a.reps, seed=a.seed))
    R = chunked.merge_parts(a.out)
    R.to_parquet(a.out + ".parquet")
    print(f"written {a.out}.parquet: {len(R)} rows; failures {R.err_msg.notna().sum() if 'err_msg' in R else 0}")


if __name__ == "__main__":
    main()
