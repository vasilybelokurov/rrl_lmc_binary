"""Compare our LMC candidates (Part C tiers) with the Galactic-bulge LTTE candidates of Hajdu et al. 2021 (ApJ 915, 50;
table parsed by scripts/parse_hajdu2021.py).

Quantities common to both: P_orb, LTTE semi-amplitude A = a1 sin i / c [s] (1 au = 499.005 s), e, M2,min (M1 = 0.65, i = 90),
K1; precision metrics: relative uncertainty of A (posterior sd for Hajdu; bootstrap half 16-84% width for us), of P_orb,
and the absolute uncertainty of e; orbital cycles covered = baseline / P_orb. Hajdu's per-star baseline is not tabulated:
ASSUMED 6700 d (OGLE-III 2001.4 to OGLE-IV 2019.7). Ours: per-star baseline from stats_v3.

Output: results/partC/hajdu_comparison.csv (summary), plots/partC/hajdu_comparison.png.

Usage: python scripts/compare_hajdu.py
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

AU_S = 499.005
H_BASELINE = 6700.0


def main():
    h = pd.read_csv("data/external/hajdu2021_binprop.csv")
    h["A_s"] = h.a1sini_au * AU_S
    h["A_relerr"] = h.a1sini_err / h.a1sini_au
    h["P_relerr"] = h.Porb_err / h.Porb
    h["cycles"] = H_BASELINE / h.Porb
    h["group"] = "Hajdu " + h.quality
    o = pd.read_csv("results/partC/partC_tiers.csv")
    o = o[np.isfinite(o.P_kep) & ~o.e_at_bound]
    o["A_s"] = o.A_kep_s
    o["A_relerr"] = (o.A_s_p84 - o.A_s_p16) / 2 / o.A_kep_s
    o["P_relerr"] = (o.P_p84 - o.P_p16) / 2 / o.P_kep
    o["e_err"] = (o.e_p84 - o.e_p16) / 2
    o["cycles"] = o.baseline / o.P_kep
    o["group"] = "LMC Tier " + o.tier.astype(str)
    o = o.rename(columns={"P_kep": "Porb", "e_kep": "e", "M2min_kep": "Msmin", "K1_kep": "K1"})
    keep = ["group", "Porb", "A_s", "e", "Msmin", "K1", "A_relerr", "P_relerr", "e_err", "cycles"]
    allc = pd.concat([h[keep], o[keep]], ignore_index=True)
    order = ["Hajdu Q1", "Hajdu Q2", "Hajdu Q3", "LMC Tier 1", "LMC Tier 2", "LMC Tier 3"]
    summ = allc.groupby("group")[keep[1:]].median().reindex(order)
    summ.insert(0, "N", allc.groupby("group").size().reindex(order))
    summ.round(3).to_csv("results/partC/hajdu_comparison.csv")
    print(summ.round(3).to_string())
    for col in ["A_relerr", "P_relerr", "e_err", "cycles"]:
        a = allc[allc.group == "Hajdu Q1"][col]
        b = allc[allc.group == "LMC Tier 1"][col]
        print(f"{col:9s}: Hajdu Q1 p25/50/75 {np.nanpercentile(a, [25, 50, 75]).round(3)}  vs  LMC Tier 1 {np.nanpercentile(b, [25, 50, 75]).round(3)}")

    col = {"Hajdu Q1": "C2", "Hajdu Q2": "C1", "Hajdu Q3": "C3", "LMC Tier 1": "k", "LMC Tier 2": "0.45", "LMC Tier 3": "0.75"}
    mk = {"Hajdu Q1": "o", "Hajdu Q2": "o", "Hajdu Q3": "o", "LMC Tier 1": "*", "LMC Tier 2": "s", "LMC Tier 3": "x"}
    fig, ax = plt.subplots(2, 3, figsize=(16, 9.5))
    ax = ax.ravel()
    pp = np.geomspace(500, 2e4, 100)
    G, MSUN, C = 6.6743e-11, 1.98841e30, 299792458.0
    for m2, ls in [(0.067, ":"), (0.2, "--"), (0.6, "-")]:
        a1 = (G * (0.65 + m2) * MSUN * (pp * 86400) ** 2 / (4 * np.pi ** 2)) ** (1 / 3) * m2 / (0.65 + m2) / C
        ax[0].plot(pp, a1, "k", ls=ls, lw=0.8, label=f"edge-on M2={m2}")
    for g in order:
        d = allc[allc.group == g]
        ax[0].scatter(d.Porb, d.A_s, s=40 if "Tier 1" in g else 18, c=col[g], marker=mk[g], label=f"{g} ({len(d)})")
    ax[0].set(xscale="log", yscale="log", xlabel="P_orb [d]", ylabel="a1 sin i / c [s]", title="(a) amplitude vs period")
    ax[0].legend(fontsize=7)
    panels = [("A_relerr", "relative uncertainty of A", True), ("P_relerr", "relative uncertainty of P_orb", True),
              ("e_err", "uncertainty of e", False), ("cycles", "orbital cycles covered (Hajdu: assumed 6700-d baseline)", False),
              ("Msmin", "M2,min [Msun]", True)]
    for a_, (c_, lab, logx) in zip(ax[1:], panels):
        vals = allc[c_].replace([np.inf, -np.inf], np.nan).dropna()
        bins = np.geomspace(max(vals[vals > 0].min(), 1e-3), vals.max(), 25) if logx else np.linspace(0, vals.max(), 25)
        for g in ["Hajdu Q1", "Hajdu Q2", "Hajdu Q3", "LMC Tier 1", "LMC Tier 2"]:
            d = allc[allc.group == g][c_].dropna()
            a_.hist(d, bins, histtype="step", lw=2 if "Tier 1" in g or "Q1" in g else 1, color=col[g], label=g, density=False)
        a_.set(xlabel=lab, ylabel="N", xscale="log" if logx else "linear")
    ax[1].legend(fontsize=7)
    fig.tight_layout()
    Path("plots/partC").mkdir(parents=True, exist_ok=True)
    fig.savefig("plots/partC/hajdu_comparison.png", dpi=110)


if __name__ == "__main__":
    main()
