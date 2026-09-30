"""Summarize an injection-recovery run (scripts/inject_recover.py output).

Reports: the null distribution of D and thresholds at per-star false-alarm probabilities; false-positive rates of the
nuisance classes; LTTE recovery vs (P_orb, semi-amplitude) and vs (P_orb, M2); period accuracy of recovered orbits;
and the real stars compared with the nulls (D, jitter). Writes a figure.

Usage
-----
    python scripts/analyze_inject.py results/inject/run1.parquet --fig plots/inject_run1.png
"""
from __future__ import annotations

import argparse

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("file")
    ap.add_argument("--fig", default=None)
    ap.add_argument("--stat", default="D", help="detection statistic column (e.g. D, or D_O for OGLE-only)")
    ap.add_argument("--pcol", default="P_best", help="best-period column matching --stat")
    ap.add_argument("--alpha-veto", type=float, default=2.0)
    a = ap.parse_args()
    d = pd.read_parquet(a.file)
    d = d[d.kind != "error"].copy()
    d["D"], d["P_best"] = d[a.stat], d[a.pcol]
    print("rows:", d.kind.value_counts().to_dict())
    null = d[d.kind == "null"]
    thr = {f: float(np.quantile(null.D, 1 - f)) for f in [0.05, 0.01, 0.001]}
    print("null D quantiles (per-star FAP -> threshold):", {k: round(v, 1) for k, v in thr.items()})
    D1 = thr[0.01]

    print(f"\nfraction with D > D(FAP 1%) = {D1:.1f}:")
    for k in ["null", "blazhko", "jump", "rwalk", "ltte", "real"]:
        s = d[d.kind == k]
        print(f"  {k:8s} N={len(s):5d}  frac={np.mean(s.D > D1):.3f}  with alpha veto={np.mean((s.D > D1) & (s.alpha_chi2nu < a.alpha_veto)):.3f}  median D={s.D.median():.1f}  "
              f"median jitter0={s.jit0_s.median():.0f} s  median alpha_chi2nu={s.alpha_chi2nu.median():.2f}")

    b = d[d.kind == "blazhko"]
    print("\nBlazhko false positives vs P_B:")
    for lo, hi in [(20, 100), (100, 300), (300, 1000), (1000, 3000)]:
        s = b[(b.P_B >= lo) & (b.P_B < hi)]
        print(f"  P_B [{lo},{hi}) d: N={len(s)}, frac D>D1 = {np.mean(s.D > D1):.2f}, "
              f"median alpha_chi2nu = {s.alpha_chi2nu.median():.2f}")

    L = d[d.kind == "ltte"].copy()
    L["det"] = (L.D > D1) & (L.alpha_chi2nu < a.alpha_veto)
    print(f"\n(recovery below uses D > D1 AND alpha_chi2nu < {a.alpha_veto})")
    L["snr"] = L.amp_s / L.err_med_s
    print("\nLTTE recovery vs amplitude / median season error (amp_s / err_med_s):")
    for lo, hi in [(0, 0.5), (0.5, 1), (1, 2), (2, 4), (4, 100)]:
        s = L[(L.snr >= lo) & (L.snr < hi)]
        print(f"  [{lo},{hi}): N={len(s)}, recovered {s.det.mean():.2f}")
    print("\nLTTE recovery vs P_orb and M2 (fraction; N):")
    Pb, Mb = [300, 1000, 3000, 10000], [0.05, 0.15, 0.4, 1.5]
    for i in range(3):
        row = []
        for j in range(3):
            s = L[(L.P_orb >= Pb[i]) & (L.P_orb < Pb[i + 1]) & (L.M2 >= Mb[j]) & (L.M2 < Mb[j + 1])]
            row.append(f"{s.det.mean():.2f} ({len(s)})")
        print(f"  P {Pb[i]}-{Pb[i + 1]} d:  " + "   ".join(f"M2 {Mb[j]}-{Mb[j + 1]}: {row[j]}" for j in range(3)))
    det = L[L.det & (L.P_orb < 6000)]
    if len(det):
        rel = np.abs(det.P_best / det.P_orb - 1)
        print(f"\nrecovered, P_orb < 6000 d: |P_best/P_orb - 1| < 0.2 for {np.mean(rel < 0.2):.2f} of {len(det)}")

    r = d[d.kind == "real"]
    print(f"\nreal: N={len(r)}, frac D > D1 {np.mean(r.D > D1):.3f}, > D(0.1%) {np.mean(r.D > thr[0.001]):.3f}; "
          f"jitter0 p50/p90 = {r.jit0_s.median():.0f}/{r.jit0_s.quantile(.9):.0f} s vs null p90 {null.jit0_s.quantile(.9):.0f} s")

    if a.fig:
        fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))
        bins = np.linspace(0, 80, 41)
        for k, c in [("null", "k"), ("blazhko", "C1"), ("jump", "C2"), ("rwalk", "C4"), ("ltte", "C0"), ("real", "C3")]:
            s = d[d.kind == k]
            ax[0].hist(np.clip(s.D, 0, 79.9), bins, histtype="step", density=True, color=c, label=k, lw=1.5)
        ax[0].axvline(D1, ls=":", c="k")
        ax[0].set(xlabel="D = 2 max ΔlnL (circular LTTE vs quadratic)", ylabel="density", yscale="log")
        ax[0].legend(fontsize=8)
        sc = ax[1].scatter(L.P_orb, L.amp_s, c=L.det, cmap="coolwarm", s=8, vmin=0, vmax=1)
        ax[1].set(xscale="log", yscale="log", xlabel="P_orb [d]", ylabel="injected a1 sin i / c [s]",
                  title="LTTE injections: red = recovered at FAP 1%")
        ax[2].hist(null.jit0_s, np.linspace(0, 2000, 41), histtype="step", density=True, color="k", label="null")
        ax[2].hist(r.jit0_s, np.linspace(0, 2000, 41), histtype="step", density=True, color="C3", label="real")
        ax[2].set(xlabel="extra timing jitter under H0 [s]", ylabel="density")
        ax[2].legend()
        fig.tight_layout()
        fig.savefig(a.fig, dpi=120)
        print("figure:", a.fig)


if __name__ == "__main__":
    main()
