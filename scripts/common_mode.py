"""Common-mode timing test: stack all stars' O-C residuals (after each star's own quadratic + MACHO offset, white-jitter
WLS) by observing year. A time-system error common to all stars (survey/season timestamps) appears as a coherent
residual per year; independent orbits and noise average out.

Usage
-----
    python scripts/common_mode.py --src results/real/oc_all.parquet --fig plots/common_mode.png
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from rrlbin.ltte import _profile_lnl  # noqa: E402
from rrlbin.timing import year_labels  # noqa: E402

DAY = 86400.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", default="results/real/oc_all.parquet")
    ap.add_argument("--fig", default="plots/common_mode.png")
    ap.add_argument("--null", action="store_true",
                    help="replace each star's delays by white noise (its own errors + fitted jitter): tests whether the fit alone creates a pattern")
    a = ap.parse_args()
    rng = np.random.default_rng(0)
    d = pd.read_parquet(a.src)
    d = d[d.ok & (d.jit0_s < 500)]           # well-behaved stars only
    rows = []
    for _, r in d.iterrows():
        t, tau, err, flag = map(np.asarray, (r.t, r.tau, r.err, r.flag))
        if a.null:
            tau = rng.normal(0, 1, t.size) * np.sqrt(err ** 2 + (r.jit0_s / DAY) ** 2)
        X = np.vander((t - t.mean()) / 1000.0, 3)
        if flag.any():
            X = np.column_stack([X, flag])
        s = np.r_[0.0, np.geomspace(0.1, 30, 20) * np.median(err)]
        _, s2, b = _profile_lnl(X, tau, err ** 2, s ** 2)
        res = (tau - X @ b) * DAY
        w = 1 / (err ** 2 + s2) / DAY ** 2
        for yr, rr, ww, fl in zip(year_labels(t), res, w, flag):
            rows.append((yr, "MACHO" if fl else "OGLE", rr, ww))
    R = pd.DataFrame(rows, columns=["year", "survey", "res", "w"])
    g = R.groupby(["survey", "year"]).apply(lambda x: pd.Series(dict(
        n=len(x), wmean=np.sum(x.w * x.res) / x.w.sum(), err=x.w.sum() ** -0.5,
        median=np.median(x.res), sem_rob=1.2533 * 1.4826 * np.median(np.abs(x.res - np.median(x.res))) / np.sqrt(len(x)))),
        include_groups=False).reset_index()
    g["hjd_mid"] = 245 + 365.25 * (g.year + 0.5)
    print(g.round(1).to_string(index=False))
    fig, ax = plt.subplots(figsize=(9, 4))
    for sv, c in [("MACHO", "C1"), ("OGLE", "C0")]:
        x = g[g.survey == sv]
        ax.errorbar(x.hjd_mid, x["median"], x.sem_rob, fmt="o", color=c, label=f"{sv} (median ± robust s.e.)")
    ax.axhline(0, c="k", lw=0.5)
    ax.set(xlabel="HJD − 2450000", ylabel="stacked O−C residual [s]", title=f"Common-mode timing test ({len(d)} stars)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(a.fig, dpi=110)


if __name__ == "__main__":
    main()
