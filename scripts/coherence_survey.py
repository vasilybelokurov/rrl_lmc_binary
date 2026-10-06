"""Stage 1b (docs/PLAN.md): how coherent is the intrinsic timing modulation of real RR Lyrae?

For a random sample of real stars (+ all v4 candidates and frozen candidates), the v4 season delays (common mode applied,
robust unwrap and band alignment as in Level 2) are fitted with
  white noise only (trend + offsets + jitter), H_RN (smooth red noise, SE kernel), H_QP (quasi-periodic, lnL for each
  coherence c = l/P_q in {0.5, 1, 2, 4, inf}).
The distribution of the preferred coherence among stars WITH a significant periodic component tells how coherent intrinsic
modulation is (the bound needed by the hypothesis test, stage 1 result in JOURNAL.md).

Output (chunked, resumable): results/stage1/coherence_survey/ -> merged results/stage1/coherence_survey.parquet.
Usage: python scripts/coherence_survey.py --n-random 3000 --workers 3
"""
from __future__ import annotations

import argparse
import json
import sys

sys.path.insert(0, __import__("os").path.dirname(__file__))
import chunked  # noqa: E402  (single-threaded BLAS)
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from level2 import load, priors_for, to_series  # noqa: E402
from rrlbin.models import COHERENCE, fit_qp, fit_rn, fit_white  # noqa: E402
from rrlbin.oc import align_bands, apply_common_mode, robust_unwrap  # noqa: E402
from rrlbin.timing import year_labels  # noqa: E402

DAY = 86400.0
CM = {tuple(map(int, k.split(","))): v / DAY for k, v in json.loads(Path("results/calib/common_mode_v4.json").read_text()).items()}
LAG = json.loads(Path("results/calib/band_lag_v4.json").read_text())


def work(rec):
    try:
        s = apply_common_mode(to_series(rec), CM, year_labels)
        t, y, e, band = s["t"], s["tau"], s["err"], s["band"].astype(int)
        o = np.argsort(t)
        t, y, e, band = t[o], y[o], e[o], band[o]
        pri = priors_for(s["P"], LAG)
        y = align_bands(t, robust_unwrap(t, y, band, s["P"]), e, band, s["P"], pri)
        e0 = np.median(e)
        lw = fit_white(t, y, e, band, pri)["lnl"]
        rn = fit_rn(t, y, e, band, pri, n_amp=10)
        qp = fit_qp(t, y, e, band, pri, n_amp=10)
        out = dict(ogle_id=rec["ogle_id"], n=int(t.size), baseline=float(np.ptp(t)), err_med_s=e0 * DAY, has_M=bool(np.any(band > 0)),
                   lnl_white=lw, lnl_rn=rn["lnl"], A_rn_s=rn["A"] * DAY, ell_rn=rn["ell"], lnl_qp=qp["lnl"],
                   Pq=qp["best"]["Pq"] if qp["best"] else np.nan, c_best=qp["best"]["c"] if qp["best"] else np.nan,
                   A_qp_s=qp["best"]["A"] * DAY if qp["best"] else np.nan)
        out.update({f"lnl_qp_c{c}": v for c, v in qp["lnl_by_c"].items()})
        return out
    except Exception as ex:
        return dict(ogle_id=rec["ogle_id"], err_msg=str(ex)[:200])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-random", type=int, default=3000)
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--seed", type=int, default=11)
    ap.add_argument("--out", default="results/stage1/coherence_survey")
    a = ap.parse_args()
    d = load("results/real/series_v4")
    d = d[d.ok.astype(bool)]
    d = d[d.t.apply(len) >= 12]
    extra = set(pd.read_csv("results/partC_v4/kepler_targets.csv").ogle_id)
    rnd = set(d.ogle_id.sample(a.n_random, random_state=a.seed))
    sel = d[d.ogle_id.isin(rnd | extra)].copy()
    sel["in_random"] = sel.ogle_id.isin(rnd)
    recs = sel.to_dict("records")
    print(f"stars {len(recs)} (random {len(rnd)}, candidates {len(extra)}, overlap {len(rnd & extra)}); coherence grid {COHERENCE}", flush=True)
    chunked.run_chunked(work, recs, a.out, workers=a.workers, chunk=100,
                        manifest=dict(script="coherence_survey", ids=[r["ogle_id"] for r in recs], seed=a.seed))
    R = chunked.merge_parts(a.out).merge(sel[["ogle_id", "in_random"]], on="ogle_id", how="left")
    R["is_candidate"] = R.ogle_id.isin(extra)
    R.to_parquet(a.out + ".parquet")
    print(f"written {a.out}.parquet: {len(R)} rows; failures {R.err_msg.notna().sum() if 'err_msg' in R else 0}")


if __name__ == "__main__":
    main()
