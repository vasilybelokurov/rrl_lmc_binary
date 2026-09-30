"""Summary statistics of the O-C fits for all real RRab, the candidate selection, and simulated classes for reference.

Candidate cuts (applied in this order; the same function is applied to sims):
  1. D > 40                              (circular LTTE vs quadratic, per-star FAP << 0.1% for white noise)
  2. chi2_nu(alpha) < 2                  (no significant season-to-season amplitude change: Blazhko veto)
  3. amp / sigma_season > 3
  4. baseline / P_best > 1.5             (at least 1.5 orbital cycles observed)
  5. amp < a1/c(M1 = 0.65, M2 = 2 Msun)  (physical LTTE ceiling at the fitted period)

Writes plots/summary_stats.png and results/real/candidates.csv (with f(M), M2,min).

Usage
-----
    python scripts/plot_summary_stats.py --src results/real/oc_all_cm.parquet
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import brentq

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from rrlbin.ltte import a1sini_over_c, mass_function  # noqa: E402

D_MIN, ALPHA_MAX, SNR_MIN, NCYC_MIN = 40.0, 2.0, 3.0, 1.5


def cut_flags(x: pd.DataFrame, amp: str, baseline) -> pd.DataFrame:
    ceil = a1sini_over_c(x.P_best.to_numpy(float), 0.65, 2.0)
    f = pd.DataFrame(index=x.index)
    f["c1_D"] = x.D > D_MIN
    f["c2_alpha"] = x.alpha_chi2nu < ALPHA_MAX
    f["c3_snr"] = x[amp] / x.err_med_s > SNR_MIN
    f["c4_cycles"] = baseline / x.P_best > NCYC_MIN
    f["c5_ceiling"] = x[amp] < ceil
    f["all"] = f.all(axis=1)
    return f


def m2_min(f):
    return brentq(lambda m: m ** 3 / (0.65 + m) ** 2 - f, 1e-4, 100)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", default="results/real/oc_all_cm.parquet")
    ap.add_argument("--sims", nargs="+", default=["results/inject/macho_v2.parquet"])
    ap.add_argument("--fig", default="plots/summary_stats.png")
    a = ap.parse_args()
    r = pd.read_parquet(a.src)
    r = r[r.ok].drop(columns=["t", "tau", "err", "flag", "alpha", "alpha_err"], errors="ignore").copy()
    fr = cut_flags(r, "amp_s", r.baseline)
    r = r.join(fr)
    s = pd.concat([pd.read_parquet(f) for f in a.sims], ignore_index=True)
    s = s[~s.kind.isin(["error", "real"])].drop(columns=["t", "tau", "err", "flag"], errors="ignore").copy()
    s = s.join(cut_flags(s, "amp_best_s", s.baseline if "baseline" in s else 8350.0))

    c = r[r["all"]].copy()
    c["fM"] = mass_function(c.P_best, c.amp_s)
    c["M2min"] = [m2_min(v) for v in c.fM]
    c["K1_kms"] = 2 * np.pi * c.amp_s * 299792.458 / (c.P_best * 86400)
    c["snr_tot"] = c.amp_s / np.hypot(c.err_med_s, c.jit0_s)
    c = c.sort_values("snr_tot", ascending=False)
    c.to_csv("results/real/candidates.csv", index=False)

    # funnel
    print("selection funnel (real; all / with MACHO):")
    m = np.ones(len(r), bool)
    for k in ["c1_D", "c2_alpha", "c3_snr", "c4_cycles", "c5_ceiling"]:
        m &= r[k].to_numpy()
        print(f"  after {k:10s}: {m.sum():6d} / {(m & r.has_M.to_numpy()).sum():5d}")
    print("sim pass fractions (all cuts):", s.groupby("kind")["all"].mean().round(4).to_dict())

    grey, red = "0.55", "C3"
    fig, ax = plt.subplots(2, 4, figsize=(20, 9.5))
    ax = ax.ravel()
    # (a) D distributions
    bins = np.linspace(0, 100, 51)
    ax[0].hist(np.clip(r.D, 0, 99.9), bins, density=True, histtype="stepfilled", color=grey, alpha=0.5, label=f"real (N={len(r)})")
    for k, col in [("null", "k"), ("jump_big", "C2"), ("blazhko", "C1"), ("ltte", "C0")]:
        ax[0].hist(np.clip(s[s.kind == k].D, 0, 99.9), bins, density=True, histtype="step", color=col, lw=1.4, label=f"sim {k}")
    ax[0].axvline(D_MIN, c=red, ls="--")
    ax[0].set(yscale="log", xlabel="D = 2ΔlnL (orbit vs quadratic)", ylabel="density", title="(a) detection statistic")
    ax[0].legend(fontsize=8)
    # (b) D vs P_best
    ax[1].scatter(r.P_best, r.D, s=2, c=grey, alpha=0.4, rasterized=True)
    ax[1].scatter(c.P_best, c.D, s=14, c=red, label=f"candidates ({len(c)})")
    ax[1].axhline(D_MIN, c=red, ls="--")
    ax[1].set(xscale="log", yscale="log", xlabel="best orbital period P [d]", ylabel="D", title="(b) D vs period")
    ax[1].legend(fontsize=8)
    # (c) amplitude vs period with physical curves
    pp = np.geomspace(250, 2e4, 200)
    ax[2].scatter(r.P_best, r.amp_s, s=2, c=grey, alpha=0.3, rasterized=True)
    for M2, ls in [(2.0, "-"), (0.5, "--"), (0.15, ":")]:
        ax[2].plot(pp, a1sini_over_c(pp, 0.65, M2), "k", ls=ls, lw=1, label=f"edge-on, M2 = {M2} Msun")
    ax[2].scatter(c.P_best, c.amp_s, s=14, c=red)
    ax[2].set(xscale="log", yscale="log", ylim=(20, 3e4), xlabel="P [d]", ylabel="fitted semi-amplitude a1 sin i / c [s]",
              title="(c) amplitude vs period (solid = ceiling)")
    ax[2].legend(fontsize=8)
    # (d) signal-to-noise vs jitter
    ax[3].scatter(r.jit0_s / r.err_med_s, r.amp_s / r.err_med_s, s=2, c=grey, alpha=0.3, rasterized=True)
    ax[3].scatter(c.jit0_s / c.err_med_s, c.amp_s / c.err_med_s, s=14, c=red)
    ax[3].axhline(SNR_MIN, c=red, ls="--")
    ax[3].set(xscale="symlog", yscale="log", xlabel="extra jitter under H0 / σ_season", ylabel="amplitude / σ_season",
              title="(d) signal vs timing noise")
    # (e) amplitude constancy
    b2 = np.geomspace(0.05, 100, 50)
    ax[4].hist(r.alpha_chi2nu.clip(0.05, 99), b2, histtype="stepfilled", color=grey, alpha=0.5, label="real")
    ax[4].hist(r[r.c1_D].alpha_chi2nu.clip(0.05, 99), b2, histtype="step", color="k", label="real, D > 40")
    ax[4].hist(c.alpha_chi2nu, b2, histtype="stepfilled", color=red, label="candidates")
    ax[4].axvline(ALPHA_MAX, c=red, ls="--")
    ax[4].set(xscale="log", yscale="log", xlabel="χ²_ν of season amplitudes α_j", ylabel="stars", title="(e) Blazhko veto")
    ax[4].legend(fontsize=8)
    # (f) cycles covered
    b3 = np.geomspace(0.3, 40, 50)
    ax[5].hist(r.baseline / r.P_best, b3, histtype="stepfilled", color=grey, alpha=0.5, label="real")
    ax[5].hist(c.baseline / c.P_best, b3, histtype="stepfilled", color=red, label="candidates")
    ax[5].axvline(NCYC_MIN, c=red, ls="--")
    ax[5].set(xscale="log", yscale="log", xlabel="baseline / P_best (orbital cycles covered)", ylabel="stars", title="(f) cycles")
    ax[5].legend(fontsize=8)
    # (g) funnel incl. sims
    steps = ["c1_D", "c2_alpha", "c3_snr", "c4_cycles", "c5_ceiling"]
    labels = ["D>40", "+α veto", "+S/N>3", "+≥1.5 cyc", "+ceiling"]
    for k, col, lab in [(None, red, "real"), ("ltte", "C0", "sim LTTE"), ("jump_big", "C2", "sim large jump"),
                        ("rwalk_big", "C5", "sim large rwalk"), ("blazhko", "C1", "sim Blazhko"), ("rwalk", "C4", "sim rwalk"),
                        ("jump", "C8", "sim jump"), ("null", "k", "sim null")]:
        x = r if k is None else s[s.kind == k]
        m = np.ones(len(x), bool)
        fr_ = []
        for st in steps:
            m &= x[st].to_numpy()
            fr_.append(max(m.mean(), 1e-5))
        ax[6].plot(labels, fr_, "o-", color=col, label=lab)
    ax[6].set(yscale="log", ylabel="fraction passing", title="(g) selection funnel")
    ax[6].legend(fontsize=8)
    # (h) candidates: M2,min vs P
    sc = ax[7].scatter(c.P_best / 365.25, c.M2min, c=np.log10(c.D), s=30, cmap="viridis")
    ax[7].set(xlabel="P_orb [yr]", ylabel="M2,min [Msun] (M1 = 0.65, i = 90°)", yscale="log",
              title="(h) candidates: minimum companion mass")
    for m2, lab in [(0.6, "0.6"), (0.2, "0.2"), (0.067, "0.067")]:
        ax[7].axhline(m2, c="0.6", lw=0.8, ls=":")
        ax[7].text(ax[7].get_xlim()[0] if False else c.P_best.min() / 365.25, m2 * 1.05, f"Hajdu+21 {lab}", fontsize=7, color="0.4")
    plt.colorbar(sc, ax=ax[7], label="log10 D")
    fig.tight_layout()
    fig.savefig(a.fig, dpi=110)
    print("figure:", a.fig, "; candidates:", len(c), "-> results/real/candidates.csv")


if __name__ == "__main__":
    main()
