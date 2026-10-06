"""Part B summary figure: real MACHO stars vs the v3 simulations (noise model, contamination, completeness).

Panels:
  (a) survival function P(D > x) of the detection statistic: real MACHO stars vs simulated classes;
  (b) fraction of each simulated class passing the candidate cuts and the Tier-1 proxy;
  (c) observed vs noise-predicted numbers of real candidates at three selection stages, the amplitude (alpha) veto corrected
      with the real pass fraction (the empirical class has no amplitude modulation, so it passes the veto ~100%);
  (d) completeness of the candidate cuts for injected LTTE orbits (MACHO cadences) by orbital period and companion mass.

Inputs: results/real/stats_v3.parquet, results/inject/stats_v3_macho.parquet, results/partB/completeness.csv
Output: plots/partB_summary.png; prints the numbers shown in (c).
Usage: python scripts/plot_partB.py
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_summary_stats import cut_flags  # noqa: E402
sys.path.insert(0, __import__('os').path.dirname(__file__))
from version import CM, LAG, PARTB, PARTC, PLOTS_C, PREV, PREV_CANDS, SERIES, SFX, SIM_M, SIM_O, STATS, V  # noqa: E402,F401

OTHERS = ["c1_D", "c3_snr", "c4_cycles", "c5_ceiling"]          # all cuts except the amplitude veto c2_alpha
CLASSES = [("null", "k", "white noise"), ("rwalk", "C4", "random walk"), ("jump", "C8", "period jump"),
           ("blazhko", "C1", "Blazhko"), ("rwalk_big", "C5", "large random walk"), ("jump_big", "C2", "large jump"),
           ("empirical", "C3", "empirical noise"), ("ltte", "C0", "LTTE orbit")]


def main():
    r = pd.read_parquet(STATS)
    r = r[r.ok & r.has_M].copy()
    r = r.join(cut_flags(r, "amp_s", r.baseline))
    s = pd.read_parquet(SIM_M,
                        columns=["kind", "D", "amp_circ_s", "err_med_s", "P_best", "baseline", "alpha_chi2nu", "coh_chi2nu",
                                 "pred_score", "P_orb", "M2"])
    s = s[s.kind != "error"].copy()
    s = s.join(cut_flags(s, "amp_circ_s", s.baseline))
    for d in (r, s):
        d["oth"] = d[OTHERS].all(axis=1)
        d["t1_noalpha"] = d.oth & (d.coh_chi2nu < 1.5) & (d.pred_score > 2)
        d["t1"] = d.t1_noalpha & d.c2_alpha & (d.alpha_chi2nu < 1.5)
    e = s[s.kind == "empirical"]
    N = len(r)
    # (c) numbers: noise prediction = empirical pass rate x N; alpha veto corrected by the real stars' pass fraction
    f_alpha = r.loc[r.oth, "c2_alpha"].mean()
    f_alpha_t1 = (r.loc[r.t1_noalpha, "c2_alpha"] & (r.loc[r.t1_noalpha, "alpha_chi2nu"] < 1.5)).mean()
    stages = [("all cuts\nexcept α veto", r.oth.sum(), e.oth.mean() * N),
              ("all cuts\n(α corrected)", r["all"].sum(), e.oth.mean() * N * f_alpha),
              ("Tier-1 proxy\n(α corrected)", r.t1.sum(), e.t1_noalpha.mean() * N * f_alpha_t1)]
    print(f"N real MACHO = {N}; real alpha-veto pass fraction {f_alpha:.2f} (Tier-1 proxy {f_alpha_t1:.2f})")
    for lab, o, x in stages:
        print(f"  {lab.replace(chr(10), ' ')}: observed {o}, expected from noise {x:.0f} (Poisson sd {np.sqrt(x):.0f})")

    fig, ax = plt.subplots(2, 2, figsize=(12, 9.5))
    ax = ax.ravel()
    # (a) survival functions of D
    x = np.linspace(0, 120, 241)
    sf = lambda v: np.array([(v > xx).mean() for xx in x])
    ax[0].plot(x, sf(r.D.to_numpy()), color="0.3", lw=3, alpha=0.6, label=f"real MACHO stars (N={N})")
    for k, col, lab in CLASSES:
        if k in ("null", "empirical", "ltte", "jump_big", "blazhko"):
            ax[0].plot(x, sf(s[s.kind == k].D.to_numpy()), color=col, lw=1.4, label=f"sim: {lab}")
    ax[0].axvline(40, c="r", ls="--", lw=1)
    ax[0].set(yscale="log", ylim=(1e-4, 1.2), xlabel="D = 2ΔlnL (orbit vs quadratic)", ylabel="P(D > x)",
              title="(a) detection statistic: real vs simulated")
    ax[0].legend(fontsize=8)
    # (b) pass rates
    xi = np.arange(len(CLASSES))
    pc = [max(s[s.kind == k]["all"].mean(), 1e-4) for k, _, _ in CLASSES]
    pt = [max(s[s.kind == k]["t1"].mean(), 1e-4) for k, _, _ in CLASSES]
    ax[1].bar(xi - 0.2, pc, 0.4, color=[c for _, c, _ in CLASSES], label="candidate cuts")
    ax[1].bar(xi + 0.2, pt, 0.4, color=[c for _, c, _ in CLASSES], alpha=0.45, hatch="//", label="Tier-1 proxy")
    ax[1].axhline(1e-4, c="0.5", lw=0.6)
    ax[1].set_xticks(xi, [lab for _, _, lab in CLASSES], rotation=35, ha="right", fontsize=8)
    ax[1].set(yscale="log", ylim=(8e-5, 0.5), ylabel="fraction passing (floor 1e-4 = none)", title="(b) pass rates of simulated classes")
    ax[1].legend(fontsize=8)
    # (c) observed vs expected
    xs = np.arange(len(stages))
    ax[2].bar(xs - 0.2, [o for _, o, _ in stages], 0.4, color="0.4", label="observed (real)")
    ax[2].bar(xs + 0.2, [v for _, _, v in stages], 0.4, yerr=[np.sqrt(v) for _, _, v in stages], color="C3", alpha=0.7,
              capsize=4, label="expected from timing noise (empirical sims)")
    for i, (_, o, v) in enumerate(stages):
        ax[2].text(i - 0.2, o, f"{o}", ha="center", va="bottom", fontsize=9)
        ax[2].text(i + 0.2, v + np.sqrt(v), f"{v:.0f}", ha="center", va="bottom", fontsize=9)
    ax[2].set_xticks(xs, [lab for lab, _, _ in stages], fontsize=9)
    ax[2].set(ylabel="number of MACHO stars", title="(c) candidates: observed vs noise-only expectation")
    ax[2].legend(fontsize=8)
    # (d) completeness
    c = pd.read_csv(f"{PARTB}/completeness.csv")
    c = c[c["sample"] == "MACHO"]
    Ps, Ms = list(dict.fromkeys(c.P)), list(dict.fromkeys(c.M2))
    Z = np.array([[c[(c.P == p) & (c.M2 == m)].comp_cuts.iloc[0] for p in Ps] for m in Ms])
    im = ax[3].imshow(Z, origin="lower", cmap="viridis", vmin=0, vmax=0.6, aspect="auto")
    for i in range(len(Ms)):
        for j in range(len(Ps)):
            ax[3].text(j, i, f"{Z[i, j]:.2f}", ha="center", va="center", color="w" if Z[i, j] < 0.35 else "k", fontsize=10)
    ax[3].set_xticks(range(len(Ps)), [f"{p} d" for p in Ps])
    ax[3].set_yticks(range(len(Ms)), [f"{m} Msun" for m in Ms])
    ax[3].set(xlabel="injected orbital period", ylabel="injected companion mass M2",
              title="(d) completeness of the cuts (MACHO cadences, isotropic)")
    plt.colorbar(im, ax=ax[3], label="fraction recovered")
    fig.tight_layout()
    fig.savefig(f"plots/partB_summary{SFX}.png", dpi=110)
    print(f"figure: plots/partB_summary{SFX}.png")


if __name__ == "__main__":
    main()
