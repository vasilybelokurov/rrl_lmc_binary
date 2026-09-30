"""Plot O-C (season delays) of the strongest LTTE-like real stars, with the best circular-orbit + quadratic fit.

Selection (from results/real/oc_all parts or the merged file): D > D_min, alpha veto, amp/err_med > 3, >= 1.5 cycles,
ranked by amp / sqrt(err_med^2 + jit0^2) (signal relative to the total timing noise).

Usage
-----
    python scripts/plot_oc_candidates.py --n 12 --fig figures/oc_candidates.png
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from chunked import merge_parts  # noqa: E402
from rrlbin.ltte import _profile_lnl, a1sini_over_c  # noqa: E402

DAY = 86400.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", default="results/real/oc_all")
    ap.add_argument("--n", type=int, default=12)
    ap.add_argument("--D-min", type=float, default=40)
    ap.add_argument("--fig", default="figures/oc_candidates.png")
    a = ap.parse_args()
    d = merge_parts(a.src) if Path(a.src).is_dir() else pd.read_parquet(a.src)
    d = d[d.ok].copy()
    d["snr"] = d.amp_s / d.err_med_s
    d["snr_tot"] = d.amp_s / np.hypot(d.err_med_s, d.jit0_s)
    # physical ceiling: an LTTE amplitude cannot exceed a1/c for M1 = 0.65 and M2 = 2 Msun at the fitted period
    d["amp_max_s"] = a1sini_over_c(d.P_best.to_numpy(), 0.65, 2.0)
    c = d[(d.D > a.D_min) & (d.alpha_chi2nu < 2) & (d.snr > 3) & (d.baseline / d.P_best > 1.5) & (d.amp_s < d.amp_max_s)]
    print(f"stars {len(d)}; above physical LTTE ceiling: {(d.amp_s > d.amp_max_s).sum()}; candidates: {len(c)}")
    c.drop(columns=["t", "tau", "err", "flag", "alpha", "alpha_err"]).to_csv(Path(a.fig).with_suffix(".csv"), index=False)
    c = c.sort_values("snr_tot", ascending=False).head(a.n)
    cols = ["ogle_id", "has_M", "D", "P_best", "amp_s", "err_med_s", "jit0_s", "snr_tot", "alpha_chi2nu", "pred_score"]
    print(c[cols].round(2).to_string(index=False))
    nc = 3
    nr = int(np.ceil(len(c) / nc))
    fig, axs = plt.subplots(nr, nc, figsize=(15, 3.2 * nr), squeeze=False)
    for ax, (_, r) in zip(axs.ravel(), c.iterrows()):
        t, tau, err, flag = map(np.asarray, (r.t, r.tau, r.err, r.flag))
        X = np.column_stack([np.vander((t - np.average(t, weights=err ** -2)) / 1000.0, 3)] +
                            ([flag] if flag.any() else []))
        w = 2 * np.pi * t / r.P_best
        X1 = np.column_stack([X, np.sin(w), np.cos(w)])
        s_grid = np.r_[0.0, np.geomspace(0.1, 30, 30) * np.median(err)]
        _, _, b = _profile_lnl(X1, tau, err ** 2, s_grid ** 2)
        off = b[3] * flag if flag.any() else 0.0        # remove the fitted MACHO offset
        quad = X[:, :3] @ b[:3]
        y = (tau - off - quad) * DAY
        tt = np.linspace(t.min(), t.max(), 600)
        orbit = (b[-2] * np.sin(2 * np.pi * tt / r.P_best) + b[-1] * np.cos(2 * np.pi * tt / r.P_best)) * DAY
        for fl, col, lab in [(1, "C1", "MACHO B"), (0, "C0", "OGLE I")]:
            m = flag == fl
            if m.any():
                ax.errorbar(t[m], y[m], err[m] * DAY, fmt="o", ms=3, color=col, label=lab)
        ax.plot(tt, orbit, "k-", lw=1)
        ax.set_title(f"{r.ogle_id[-5:]}  P={r.P_best:.0f} d  A={r.amp_s:.0f} s  jit={r.jit0_s:.0f} s  D={r.D:.0f}", fontsize=9)
        ax.set_xlabel("HJD − 2450000")
        ax.set_ylabel("O−C − quadratic [s]")
    axs.ravel()[0].legend(fontsize=8)
    for ax in axs.ravel()[len(c):]:
        ax.axis("off")
    fig.tight_layout()
    fig.savefig(a.fig, dpi=110)
    print("figure:", a.fig)


if __name__ == "__main__":
    main()
