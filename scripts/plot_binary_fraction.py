"""Figures for the binary-fraction write-up (docs/binary_fraction/). PNG in plots/binary_fraction/.

fig_oc_examples      O-C (season delays, 1992-2026) of a typical star, a star with strong intrinsic wander, and 13854 (candidate)
fig_D_survival       distribution of the orbit-search statistic D: real MACHO+OGLE stars vs simulated stars with measurement noise
                     only, with realistic intrinsic (red) timing noise, and with injected orbits
fig_funnel           candidates surviving each cut: real stars vs the number expected from intrinsic timing noise alone
fig_amp_period       best-fit O-C amplitude vs period for all stars, with LTTE amplitudes of companions (edge-on)
fig_blind_test       the frozen (2026-10-01) orbit predictions vs the 2017-2026 OGLE seasons for 13854 (passed) and 03269 (failed)
fig_efficiency       detection efficiency vs (P_orb, M2) from orbits injected into the real stars (needs inject_real.parquet)
fig_upper_limit      95% upper limit on the binary fraction vs P_orb for three companion-mass ranges, with literature values
Usage: python scripts/plot_binary_fraction.py [--only name ...]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from level2 import load, priors_for, to_series  # noqa: E402
from plot_summary_stats import cut_flags  # noqa: E402
from rrlbin.ltte import a1sini_over_c, ltte_delay  # noqa: E402
from rrlbin.oc import align_bands, apply_common_mode, design, profile_lnl, default_s_grid, robust_unwrap  # noqa: E402
from rrlbin.timing import year_labels  # noqa: E402

DAY = 86400.0
OUT = Path("plots/binary_fraction")
# reference palette (light mode), categorical slots in fixed order
C_REAL, C_NOISE, C_LTTE, C_NULL = "#2a78d6", "#eb6834", "#1baf7a", "#8a8984"
INK, INK2 = "#0b0b0b", "#52514e"
BAND_STYLE = {0: ("OGLE I", C_REAL, "o"), 1: ("MACHO B", C_NOISE, "s"), 2: ("MACHO R", C_LTTE, "^")}

plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": "#e6e5e0",
                     "grid.linewidth": 0.6, "legend.frameon": False, "lines.linewidth": 2, "savefig.dpi": 200})
CM = {tuple(map(int, k.split(","))): v / DAY for k, v in json.loads(Path("results/calib/common_mode_v4.json").read_text()).items()}
LAG = json.loads(Path("results/calib/band_lag_v4.json").read_text())


def prepared(rec):
    """v4 season delays as used by the search: common mode, unwrapping, band alignment; MACHO offsets removed for display."""
    s = apply_common_mode(to_series(rec), CM, year_labels)
    t, y, e, b = s["t"], s["tau"], s["err"], s["band"].astype(int)
    o = np.argsort(t)
    t, y, e, b = t[o], y[o], e[o], b[o]
    pri = priors_for(s["P"], LAG)
    y = align_bands(t, robust_unwrap(t, y, b, s["P"]), e, b, s["P"], pri)
    X, names = design(t, b, np.average(t, weights=e ** -2))
    _, _, beta, _ = profile_lnl(X, y, e, names, default_s_grid(e), pri)
    for k in names:
        if k.startswith("off"):
            y = y - beta[names.index(k)] * (b == int(k[3:]))
    return t, y * DAY, e * DAY, b


def year(t):
    """HJD - 2450000 -> calendar year (approximate, for axes)."""
    return 2000.0 + (np.asarray(t) - 1544.5) / 365.25


def save(fig, name):
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{name}.png", bbox_inches="tight")
    plt.close(fig)
    print("wrote", OUT / f"{name}.png")


def fig_oc_examples(ser, st, kep):
    m = st[st.has_M & st.err_msg.isna()] if "err_msg" in st else st[st.has_M]
    quiet = m[(m.D.between(4, 7)) & (m.err_med_s < 150) & (m.jit0_s < 50)].sort_values("err_med_s").index[0]
    noisy = m[(m.D > 150) & (m.amp_s > 3000)].sort_values("D").index[len(m[(m.D > 150) & (m.amp_s > 3000)]) // 2]
    ids = [(quiet, "a typical star: no significant timing variation"), (noisy, "intrinsic period wander (not a binary: amplitude far above any orbit)"),
           ("OGLE-LMC-RRLYR-13854", "candidate 13854 with its Keplerian fit")]
    fig, axs = plt.subplots(3, 1, figsize=(7.2, 7.6), sharex=True)
    for ax, (oid, title) in zip(axs, ids):
        t, y, e, b = prepared(dict(ser.loc[oid].to_dict(), ogle_id=oid))
        ko = np.zeros_like(t)
        if oid in kep.index:
            r = kep.loc[oid]
            ko = ltte_delay(t, r.P_kep, r.A_kep_s, r.e_kep, r.omega_kep, r.tp_kep)
        q = np.polyfit(t, y - ko, 2, w=1 / e)        # constant period change (quadratic O-C) removed for display
        y = y - np.polyval(q, t)
        for bb, (lab, col, mk) in BAND_STYLE.items():
            k = b == bb
            if k.any():
                ax.errorbar(year(t[k]), y[k], e[k], fmt=mk, ms=4, color=col, ecolor=col, elinewidth=1, label=lab)
        if oid in kep.index:
            tt = np.linspace(t.min() - 100, t.max() + 100, 800)
            ax.plot(year(tt), ltte_delay(tt, r.P_kep, r.A_kep_s, r.e_kep, r.omega_kep, r.tp_kep), color=INK, lw=1.2,
                    label=f"Keplerian orbit: P = {r.P_kep / 365.25:.1f} yr, e = {r.e_kep:.2f}")
        ax.set_title(f"{oid}: {title}", fontsize=8.5, loc="left", color=INK)
        ax.set_ylabel("O − C − quadratic [s]")
        lim = np.max(np.abs(y)) * 1.35
        ax.set_ylim(-lim, lim)
    axs[0].legend(ncol=3, fontsize=7.5, loc="lower left")
    axs[2].legend(fontsize=7.5, loc="lower left", ncol=2)
    axs[-1].set_xlabel("year")
    save(fig, "fig_oc_examples")


def fig_D_survival(real, sims):
    fig, ax = plt.subplots(figsize=(5.4, 3.8))
    Dg = np.logspace(0, 3, 200)
    for lab, x, col, ls in [("real stars (6,612)", real.D, C_REAL, "-"),
                            ("simulated: measurement noise only", sims[sims.kind == "null"].D, C_NULL, "--"),
                            ("simulated: + realistic intrinsic timing noise", sims[sims.kind == "empirical"].D, C_NOISE, "-"),
                            ("simulated: + binary orbit (broad prior)", sims[sims.kind == "ltte"].D, C_LTTE, ":")]:
        x = np.asarray(x.dropna())
        ax.plot(Dg, [(x > d).mean() for d in Dg], color=col, ls=ls, label=lab)
    ax.axvline(40, color=INK2, lw=0.8)
    ax.text(42, 0.5, "candidate\nthreshold D = 40", fontsize=7.5, color=INK2)
    ax.set(xscale="log", yscale="log", ylim=(1e-4, 1.1), xlabel="orbit-search statistic D (2 Δ ln L)", ylabel="fraction of stars with D larger")
    ax.legend(fontsize=7.5, loc="lower left")
    save(fig, "fig_D_survival")


def fig_funnel(real, sims):
    stages = [("D > 40", ["c1_D"]), ("+ S/N > 3", ["c1_D", "c3_snr"]), ("+ ≥ 1.5 cycles", ["c1_D", "c3_snr", "c4_cycles"]),
              ("+ below orbit\nceiling", ["c1_D", "c3_snr", "c4_cycles", "c5_ceiling"]),
              ("+ steady\namplitude", ["c1_D", "c3_snr", "c4_cycles", "c5_ceiling", "c2_alpha"])]
    emp = sims[sims.kind == "empirical"]
    N = len(real)
    obs, exp = [], []
    for lab, cols in stages:
        obs.append(int(real[cols].all(axis=1).sum()))
        if "c2_alpha" in cols:   # simulated noise has no amplitude variability: use the real alpha pass rate of such stars
            base = [c for c in cols if c != "c2_alpha"]
            a_pass = real[real[base].all(axis=1)].c2_alpha.mean()
            exp.append(emp[base].all(axis=1).mean() * N * a_pass)
        else:
            exp.append(emp[cols].all(axis=1).mean() * N)
    x = np.arange(len(stages))
    fig, ax = plt.subplots(figsize=(6.2, 3.6))
    ax.bar(x - 0.2, obs, 0.38, color=C_REAL, label="real stars")
    ax.bar(x + 0.2, exp, 0.38, color=C_NOISE, label="expected from intrinsic timing noise alone")
    for i in range(len(x)):
        ax.text(x[i] - 0.2, obs[i] * 1.08, f"{obs[i]}", ha="center", fontsize=7.5, color=INK)
        ax.text(x[i] + 0.2, exp[i] * 1.08, f"{exp[i]:.0f}", ha="center", fontsize=7.5, color=INK2)
    ax.set(xticks=x, yscale="log", ylabel="number of stars (of 6,612)", ylim=(30, 4000))
    ax.set_xticklabels([s[0] for s in stages], fontsize=7.5)
    ax.legend(fontsize=7.5, loc="upper right")
    save(fig, "fig_funnel")
    return dict(zip([s[0].replace("\n", " ") for s in stages], zip(obs, [round(v) for v in exp])))


def fig_amp_period(real):
    cand = real["all"]
    fig, ax = plt.subplots(figsize=(5.6, 4.0))
    ax.scatter(real.P_best[~cand] / 365.25, real.amp_s[~cand], s=2, color=C_NULL, alpha=0.25, lw=0, label="all stars: best-fit O − C sinusoid")
    ax.scatter(real.P_best[cand] / 365.25, real.amp_s[cand], s=14, color=C_REAL, lw=0, label="the 93 candidates")
    P = np.geomspace(300, 2e4, 200)
    for m2, ls in ((0.1, ":"), (0.4, "--"), (1.5, "-")):
        ax.plot(P / 365.25, a1sini_over_c(P, 0.65, m2), color=INK, lw=1, ls=ls)
        ax.text(P[-1] / 365.25 * 1.03, a1sini_over_c(P[-1], 0.65, m2), f"{m2} M$_\\odot$", fontsize=7, va="center", color=INK)
    ax.set(xscale="log", yscale="log", xlim=(0.8, 70), ylim=(10, 1e5), xlabel="best-fit period [yr]", ylabel="best-fit O − C semi-amplitude [s]")
    ax.text(1.0, 2.5e4, "lines: light-travel amplitude of an edge-on\ncompanion of the given mass (M$_1$ = 0.65 M$_\\odot$)", fontsize=7, color=INK2)
    ax.legend(fontsize=7.5, loc="lower right", markerscale=2)
    save(fig, "fig_amp_period")


def fig_blind_test():
    sys.argv = sys.argv[:1]
    import test_frozen_predictions as tf
    ser = load("results/real/series_v3").set_index("ogle_id")
    tf.G.update(cm={tuple(map(int, k.split(","))): v / DAY for k, v in json.loads(Path("results/calib/common_mode_v3.json").read_text()).items()},
                lag=json.loads(Path("results/calib/band_lag_v3.json").read_text()))
    S = pd.read_csv("results/predictions/test_2026-10-06_seasons.csv")
    R = pd.read_csv("results/predictions/test_2026-10-06.csv").set_index("ogle_id")
    fig, axs = plt.subplots(1, 2, figsize=(7.4, 3.4), sharey=False)
    for ax, oid, verdict in zip(axs, ["OGLE-LMC-RRLYR-13854", "OGLE-LMC-RRLYR-03269"], ["prediction confirmed", "prediction failed"]):
        t, y, e, band, offs = tf.frozen_frame(dict(ser.loc[oid].to_dict(), ogle_id=oid))
        z = np.load(f"results/predictions/cov_2026-10-01/cov_{oid}.npz")
        for bb, (lab, col, mk) in BAND_STYLE.items():
            k = band == bb
            if k.any():
                ax.errorbar(year(t[k]), (y[k] - offs.get(bb, 0.0)) * DAY, e[k] * DAY, fmt=mk, ms=3, color=col, ecolor=col, elinewidth=0.8)
        tg = z["t_pred"]
        s1 = np.sqrt(np.diag(z["h1_cov"]) + float(z["white_h1_s"]) ** 2)
        ax.fill_between(year(tg), z["h1_mean"] - 2 * s1, z["h1_mean"] + 2 * s1, color=INK, alpha=0.15, lw=0, label="orbit prediction (±2σ)")
        ax.plot(year(tg), z["h1_mean"], color=INK, lw=1)
        g = S[S.ogle_id == oid]
        ax.errorbar(year(g.t), g.y_s, g.err_s, fmt="D", ms=5, color=C_NOISE, ecolor=C_NOISE, label="new OGLE data 2017–2026")
        ax.axvline(year(7516), color=INK2, lw=0.8, ls=":")
        ax.set_title(f"{oid[-5:]}: {verdict}  (ln BF = {R.loc[oid, 'lnBF_all']:+.0f})", fontsize=8.5, loc="left")
        ax.set_xlabel("year")
        lo = min(np.min(y * DAY), g.y_s.min()) - 1500
        hi = max(np.max(y * DAY), g.y_s.max()) + 1500
        ax.set_ylim(lo, hi)
    axs[0].set_ylabel("O − C [s]")
    axs[0].legend(fontsize=7, loc="lower left")
    fig.text(0.5, -0.08, "dotted line: end of the data used for the predictions (frozen 1 Oct 2026, before the new data were received)",
             ha="center", fontsize=7, color=INK2)
    save(fig, "fig_blind_test")




def fig_efficiency():
    inj = pd.read_parquet("results/partB_v4/inject_real.parquet")
    inj = inj[inj.err_msg.isna()] if "err_msg" in inj else inj
    pe = np.geomspace(300, 1e4, 9)
    me = np.geomspace(0.05, 1.5, 8)
    H = np.full((len(me) - 1, len(pe) - 1), np.nan)
    for i in range(len(me) - 1):
        for j in range(len(pe) - 1):
            g = inj[(inj.M2 >= me[i]) & (inj.M2 < me[i + 1]) & (inj.P_orb >= pe[j]) & (inj.P_orb < pe[j + 1])]
            if len(g) >= 20:
                H[i, j] = g["all"].mean()
    fig, ax = plt.subplots(figsize=(5.6, 3.9))
    im = ax.pcolormesh(pe / 365.25, me, H, cmap="Blues", vmin=0, vmax=0.7, shading="flat")
    for i in range(len(me) - 1):
        for j in range(len(pe) - 1):
            if np.isfinite(H[i, j]):
                ax.text(np.sqrt(pe[j] * pe[j + 1]) / 365.25, np.sqrt(me[i] * me[i + 1]), f"{H[i, j]:.2f}", ha="center", va="center",
                        fontsize=6.5, color="white" if H[i, j] > 0.4 else INK)
    ax.set(xscale="log", yscale="log", xlabel="orbital period [yr]", ylabel="companion mass M$_2$ [M$_\\odot$]")
    ax.grid(False)
    cb = fig.colorbar(im, ax=ax)
    cb.set_label("probability that the binary is flagged")
    save(fig, "fig_efficiency")


def fig_upper_limit():
    """Final limits from results/partB_v4/upper_limits.csv (efficiency from the real-star injections with the smearing
    correction), drawn as horizontal segments over each period bin."""
    T = pd.read_csv("results/partB_v4/upper_limits.csv")
    fig, ax = plt.subplots(figsize=(5.6, 3.9))
    for m2, col in [("0.4-1.5", C_REAL), ("0.15-0.4", C_NOISE), ("0.05-0.15", C_LTTE)]:
        g = T[(T.M2 == m2) & (T.P != "1000-10000")]
        first = True
        for r in g.itertuples():
            p0, p1 = (float(x) / 365.25 for x in r.P.split("-"))
            f = r.f_max_final
            if not np.isfinite(f):
                continue
            ax.plot([p0, p1], [f, f], color=col, lw=2.5, solid_capstyle="butt",
                    label=f"M$_2$ = {m2.replace('-', '–')} M$_\\odot$" if first else None)
            ax.annotate("", xy=(np.sqrt(p0 * p1), f * 0.8), xytext=(np.sqrt(p0 * p1), f), arrowprops=dict(arrowstyle="->", color=col, lw=1))
            first = False
    ax.axhline(1.0, color=INK2, lw=0.8)
    ax.text(28, 1.08, "100%", fontsize=7, color=INK2, ha="right")
    ax.plot([9.5], [0.04], marker="^", color=INK2, ms=6, ls="none")
    ax.annotate("Hajdu+2015 (bulge, all masses):\nestimate ≳ 4% (a lower bound)", xy=(9.5, 0.042), xytext=(0.9, 0.012), fontsize=7,
                color=INK2, arrowprops=dict(arrowstyle="->", color=INK2, lw=0.8))
    ax.plot([2.5, 5.5], [0.3, 0.3], color=INK2, lw=2)
    ax.text(1.0, 0.255, "Iorio+2026: ≲ 30% (metal-rich MW disc, Gaia DR3)", fontsize=7, color=INK2)
    ax.set(xscale="log", yscale="log", xlim=(0.7, 35), ylim=(0.01, 3), xlabel="orbital period [yr]",
           ylabel="95% upper limit on the binary fraction")
    ax.legend(fontsize=7.5, loc="upper right")
    save(fig, "fig_upper_limit")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", nargs="*", default=None)
    a = ap.parse_args()
    want = lambda n: a.only is None or n in a.only  # noqa: E731
    st = pd.read_parquet("results/real/stats_v4.parquet").set_index("ogle_id")
    real = st[st.has_M].copy()
    real = real.join(cut_flags(real, "amp_s", real.baseline))
    sims = pd.read_parquet("results/inject/stats_v4_macho.parquet")
    sims = sims.join(cut_flags(sims, "amp_s", sims.baseline))
    if want("fig_oc_examples"):
        ser = load("results/real/series_v4").set_index("ogle_id")
        kep = pd.read_csv("results/partC_v4/kepler.csv").set_index("ogle_id")
        fig_oc_examples(ser, real, kep)
    if want("fig_D_survival"):
        fig_D_survival(real, sims)
    if want("fig_funnel"):
        print(fig_funnel(real, sims))
    if want("fig_amp_period"):
        fig_amp_period(real)
    if want("fig_blind_test"):
        fig_blind_test()
    if want("fig_efficiency") and Path("results/partB_v4/inject_real.parquet").exists():
        fig_efficiency()
    if want("fig_upper_limit") and Path("results/partB_v4/upper_limits_meta.json").exists():
        fig_upper_limit()


if __name__ == "__main__":
    main()
