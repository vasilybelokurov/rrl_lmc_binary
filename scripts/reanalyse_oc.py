"""Re-run the O-C statistics on saved season-delay series after a common-mode (per survey, per year) correction.

Common mode: median over well-behaved stars (jitter < 500 s) of the residual after each star's own quadratic + MACHO
offset (as scripts/common_mode.py; its null test on white-noise delays at the real epochs is consistent with zero).
It is subtracted from every star's delays; then D, P_best, amplitude, jitter, the H0 red-noise fit and the predictive
score are recomputed exactly as in real_oc.py. Light curves are not refitted.

Usage
-----
    python scripts/reanalyse_oc.py --src results/real/oc_all.parquet --out results/real/oc_all_cm.parquet
"""
from __future__ import annotations

import argparse
import sys

sys.path.insert(0, __import__("os").path.dirname(__file__))
import chunked  # noqa: E402,F401  (single-threaded BLAS)
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from rrlbin.ltte import _profile_lnl, fit_red_null, oc_search, period_grid, predictive_score  # noqa: E402
from rrlbin.timing import year_labels  # noqa: E402

DAY = 86400.0
CM = {}


def residuals(t, tau, err, flag):
    X = np.vander((t - t.mean()) / 1000.0, 3)
    if flag.any():
        X = np.column_stack([X, flag])
    s = np.r_[0.0, np.geomspace(0.1, 30, 20) * np.median(err)]
    _, _, b = _profile_lnl(X, tau, err ** 2, s ** 2)
    return tau - X @ b


def common_mode(d: pd.DataFrame) -> dict:
    rows = []
    for _, r in d[d.jit0_s < 500].iterrows():
        t, tau, err, flag = map(np.asarray, (r.t, r.tau, r.err, r.flag))
        res = residuals(t, tau, err, flag)
        rows += list(zip(year_labels(t), flag.astype(int), res))
    R = pd.DataFrame(rows, columns=["year", "flag", "res"])
    g = R.groupby(["flag", "year"]).res.agg(["median", "size"])
    return {k: v for k, v in g["median"].items() if g.loc[k, "size"] >= 100}


def one(r):
    t, tau, err, flag = (np.asarray(x) for x in (r["t"], r["tau"], r["err"], r["flag"]))
    cm = np.array([CM.get((int(f), int(y)), 0.0) for f, y in zip(flag, year_labels(t))])
    tau = tau - cm
    X = flag[:, None] if flag.any() else None
    o = oc_search(t, tau, err, period_grid(np.ptp(t)), X_extra=X)
    rn = fit_red_null(t, tau, err, X_extra=X)
    out = dict(ogle_id=r["ogle_id"], D=o["D"], P_best=o["P_best"], amp_s=o["amp"] * DAY, jit0_s=o["jit0"] * DAY,
               s_red_s=rn["s"] * DAY, rw_rms_s=float(np.sqrt(rn["q"] * np.ptp(t)) * DAY), cm_rms_s=float(np.std(cm) * DAY),
               tau=tau)
    test = (flag == 1) & (t < 450)
    if flag.any() and test.sum() >= 3 and ((flag == 1) & ~test).sum() >= 2 and r["has_O2"]:
        ps = predictive_score(t, tau, err, test, X_extra=X)
        out.update(pred_score=ps["score"], P_train=ps["P_train"], D_train=ps["D_train"])
    return out


def _init(cm):
    CM.update(cm)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", default="results/real/oc_all.parquet")
    ap.add_argument("--out", default="results/real/oc_all_cm.parquet")
    ap.add_argument("--workers", type=int, default=6)
    a = ap.parse_args()
    d = pd.read_parquet(a.src)
    d = d[d.ok]
    cm = common_mode(d)
    print("common mode [s] (flag 1 = MACHO, year):", {k: round(v * DAY, 1) for k, v in sorted(cm.items())}, flush=True)
    from multiprocessing import Pool
    recs = d[["ogle_id", "t", "tau", "err", "flag", "has_O2"]].to_dict("records")
    with Pool(a.workers, initializer=_init, initargs=(cm,)) as pool:
        rows = pool.map(one, recs, chunksize=50)
    new = pd.DataFrame(rows)
    keep = d.drop(columns=["D", "P_best", "amp_s", "jit0_s", "s_red_s", "rw_rms_s", "tau", "pred_score", "P_train", "D_train"],
                  errors="ignore")
    out = keep.merge(new, on="ogle_id")
    out.to_parquet(a.out)
    print(f"written {a.out}: {len(out)} stars; median |common mode| rms per star {out.cm_rms_s.median():.1f} s")


if __name__ == "__main__":
    main()
