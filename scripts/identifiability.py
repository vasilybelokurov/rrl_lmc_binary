"""Stage 1 (docs/PLAN.md): can a Keplerian LTTE orbit be distinguished from intrinsic quasi-periodic (QP) phase modulation
with the real season cadences (MACHO B/R + OGLE I, 1992-2026)?

For N real stars (bands I+MB+MR; real season epochs, errors, band-lag offsets), simulate season delays under
  LTTE: Keplerian, P in {2000, 4000} d, e in {0, 0.5}, random omega and periastron; semi-amplitude K = snr x median error;
  QP:   GP with qp_kernel, P_q in {2000, 4000} d, coherence c = l/P_q in {0.5, 1, 2, inf}, rms = K/sqrt(2) (same rms as LTTE);
each plus measurement noise and white jitter 0.5 x median error, snr in {4, 8}. Fit both models (rrlbin.models.fit_ltte,
fit_qp with lnl per coherence value) -> Lambda(c_max) = lnL_LTTE - max_{c <= c_max} lnL_QP.

Output: results/stage1/ident.parquet (one row per simulated data set). Analysis: scripts/identifiability_summary.py.
Usage: python scripts/identifiability.py --n-stars 30 --workers 4
"""
from __future__ import annotations

import argparse
import itertools
import json
import sys

sys.path.insert(0, __import__("os").path.dirname(__file__))
import chunked  # noqa: E402,F401  (single-threaded BLAS)
from multiprocessing import Pool  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from level2 import load  # noqa: E402
from rrlbin.ltte import ltte_delay  # noqa: E402
from rrlbin.models import COHERENCE, fit_ltte, fit_qp, qp_kernel  # noqa: E402

DAY = 86400.0


def configs():
    out = []
    for P, e, snr in itertools.product((2000.0, 4000.0), (0.0, 0.5), (4.0, 8.0)):
        out.append(dict(kind="ltte", P=P, e=e, c=np.nan, snr=snr))
    for P, c, snr in itertools.product((2000.0, 4000.0), (0.5, 1.0, 2.0, np.inf), (4.0, 8.0)):
        out.append(dict(kind="qp", P=P, e=np.nan, c=c, snr=snr))
    return out


def work(job):
    k, cfg, star, seed = job
    try:
        rng = np.random.default_rng(seed)
        t, err, band, off, pri = star["t"], star["err"], star["band"], star["off"], star["pri"]
        e0 = np.median(err)
        K = cfg["snr"] * e0
        if cfg["kind"] == "ltte":
            sig = ltte_delay(t, cfg["P"], K, cfg["e"], rng.uniform(0, 2 * np.pi), rng.uniform(0, cfg["P"]))
        else:
            C = qp_kernel(t, t, K / np.sqrt(2), cfg["P"], cfg["c"] * cfg["P"]) + 1e-12 * np.eye(t.size)
            sig = rng.multivariate_normal(np.zeros(t.size), C)
        y = sig + off + rng.normal(0, 1, t.size) * np.sqrt(err ** 2 + (0.5 * e0) ** 2)
        fl = fit_ltte(t, y, err, band, pri)
        fq = fit_qp(t, y, err, band, pri)
        row = dict(job=k, star=star["id"], **cfg, n=int(t.size), lnl_ltte=fl["lnl"], P_ltte=fl["P"], A_ltte_s=fl["A_s"], e_ltte=fl["e"],
                   lnl_qp=fq["lnl"], lnl_white=fq["lnl_white"], Pq_fit=fq["best"]["Pq"] if fq["best"] else np.nan,
                   c_fit=fq["best"]["c"] if fq["best"] else np.nan)
        row.update({f"lnl_qp_c{c}": v for c, v in fq["lnl_by_c"].items()})
        return row
    except Exception as ex:
        return dict(job=k, star=star["id"], **cfg, err_msg=str(ex)[:200])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-stars", type=int, default=30)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--out", default="results/stage1/ident.parquet")
    a = ap.parse_args()
    d = load("results/real/series_v4")
    d = d[d.ok.astype(bool) & d.bands.eq("IMBMR")].sample(a.n_stars, random_state=a.seed)
    lag = json.loads(Path("results/calib/band_lag_v4.json").read_text())
    stars = []
    for r in d.itertuples():
        t, err, band = np.asarray(r.t), np.asarray(r.err), np.asarray(r.band).astype(int)
        o = np.argsort(t)
        t, err, band = t[o], err[o], band[o]
        pri = {int(b): ((v["slope"] * r.P + v["intercept"]) / DAY, v["sd"] / DAY) for b, v in lag.items()}
        off = np.array([pri[b][0] if b in pri else 0.0 for b in band])
        stars.append(dict(id=r.ogle_id, t=t, err=err, band=band, off=off, pri=pri))
    jobs = [(k, cfg, s, a.seed * 100000 + k) for k, (cfg, s) in enumerate(itertools.product(configs(), stars))]
    print(f"stars {len(stars)}; configs {len(configs())}; jobs {len(jobs)}; coherence grid {COHERENCE}", flush=True)
    with Pool(a.workers) as pool:
        rows = []
        for i, r in enumerate(pool.imap_unordered(work, jobs, chunksize=2)):
            rows.append(r)
            if (i + 1) % 60 == 0:
                print(f"{i + 1}/{len(jobs)} done", flush=True)
    R = pd.DataFrame(rows)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    R.to_parquet(a.out)
    print(f"written {a.out}: {len(R)} rows; failures {R.err_msg.notna().sum() if 'err_msg' in R else 0}")


if __name__ == "__main__":
    main()
