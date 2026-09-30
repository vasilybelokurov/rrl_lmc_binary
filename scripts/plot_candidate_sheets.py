"""One diagnostic sheet per candidate (results/real/candidates.csv), from the light curves.

Panels:
  A  OGLE I light curve folded with the linear ephemeris + per-season delays (delay-corrected), template overlaid
  B  MACHO B, the same (if available)
  C  the companion signal seen directly in the light curve: OGLE I epochs of the seasons with the most positive and most
     negative O-C, folded with the linear ephemeris only (no delay correction), rising branch zoom; the template shifted
     by each group's mean delay is overlaid
  D  O-C (common-mode corrected; the quadratic and MACHO offset removed) with the circular-orbit fit; residuals below
  E  orbit-period search D(P) (white-jitter statistic) with the D = 40 threshold
  F  per-season amplitude factor alpha_j and mean-magnitude offset (Blazhko / blending check)

Output: plots/candidates/<rank>_<ogle_id>.png (PNG only)

Usage
-----
    python scripts/plot_candidate_sheets.py --n 67
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from rrlbin.io import read_macho, read_ogle4_ident  # noqa: E402
from rrlbin.ltte import _profile_lnl, oc_search, period_grid  # noqa: E402
from rrlbin.timing import fit_timing, fourier_eval, year_labels  # noqa: E402
from timing_sample import ogle_I  # noqa: E402

DAY = 86400.0
SEGC = {"O2": "C2", "O3": "C0", "O4": "C9", "M": "C1"}


def fold_panel(ax, t, m, seg, f, lab, title):
    """Delay-corrected fold with template."""
    tau_ep = np.zeros_like(t)
    lab_e = year_labels(t)
    for s, tau in zip(f.season, f.tau):
        tau_ep[lab_e == s] = tau
    keep = np.isin(lab_e, f.season)
    ph = np.mod((t - tau_ep - f.T0) / f.P, 1.0)
    zp = dict(zip(f.seg_names, f.zp))
    for sname in np.unique(seg[keep]):
        mm = keep & (seg == sname)
        ax.plot(ph[mm], m[mm] - zp[sname], ".", ms=1.5, color=SEGC.get(sname, "k"), alpha=0.5, label=sname)
    x = np.linspace(0, 1, 400)
    ax.plot(x, fourier_eval(f.coef, x), "k-", lw=1.2)
    ax.invert_yaxis()
    ax.set(xlabel="pulsation phase (delay-corrected)", ylabel=f"{lab} − zero point [mag]", title=title)
    ax.legend(fontsize=7, markerscale=5, loc="lower left")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cands", default="results/real/candidates.csv")
    ap.add_argument("--oc", default="results/real/oc_all_cm.parquet")
    ap.add_argument("--n", type=int, default=67)
    ap.add_argument("--outdir", default="plots/candidates")
    a = ap.parse_args()
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    c = pd.read_csv(a.cands).head(a.n)
    oc = pd.read_parquet(a.oc).set_index("ogle_id")
    inv = pd.read_parquet("data/lc_inventory.parquet").set_index("ogle_id")
    ident = read_ogle4_ident("data/raw/ogle4_lmc_rrlyr/ident.dat").set_index("ogle_id")
    par = pd.read_fwf("data/raw/ogle4_lmc_rrlyr/RRab.dat", colspecs=[(0, 20), (22, 28), (37, 46), (58, 68), (70, 75)],
                      names=["ogle_id", "I", "P", "T0", "amp_I"], header=None).set_index("ogle_id")
    for rank, (_, cr) in enumerate(c.iterrows(), 1):
        oid = cr.ogle_id
        P, T0 = float(par.loc[oid, "P"]), float(par.loc[oid, "T0"])
        tO, mO, eO, segO = ogle_I(oid)
        fO = fit_timing(tO, mO, eO, segO, P, T0, K=8, labels=year_labels(tO))
        fM, tM = None, None
        mid = inv.loc[oid, "macho_id"]
        if isinstance(mid, str) and mid.count(".") == 2:
            f_, t_, s_ = map(int, mid.split("."))
            tile = Path(f"data/raw/macho/{f_}.{t_}.parquet")
            if tile.exists():
                tM, mM, eM = read_macho(tile, s_, ident.loc[oid, "ra"], ident.loc[oid, "dec"], band="b")
                if tM.size >= 100:
                    fM = fit_timing(tM, mM, eM, np.full(tM.size, "M"), P, T0, K=8, labels=year_labels(tM))
        r = oc.loc[oid]
        t, tau, err, flag = (np.asarray(x) for x in (r.t, r.tau, r.err, r.flag))

        fig = plt.figure(figsize=(17, 10))
        gs = fig.add_gridspec(3, 3, height_ratios=[1, 1, 0.45], hspace=0.45, wspace=0.25)
        axA, axB, axC = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[0, 2])
        axD = fig.add_subplot(gs[1, :2])
        axDr = fig.add_subplot(gs[2, :2], sharex=axD)
        axE, axF = fig.add_subplot(gs[1, 2]), fig.add_subplot(gs[2, 2])

        fold_panel(axA, tO, mO, segO, fO, "I", f"A. OGLE I, delay-corrected ({fO.mask.sum()} ep.)")
        if fM is not None:
            fold_panel(axB, tM, mM, np.full(tM.size, "M"), fM, "MACHO B", f"B. MACHO B, delay-corrected ({fM.mask.sum()} ep.)")
        else:
            axB.axis("off")
            axB.text(0.5, 0.5, "no MACHO light curve", ha="center")

        # D: O-C with fit
        X = np.vander((t - np.average(t, weights=err ** -2)) / 1000.0, 3)
        if flag.any():
            X = np.column_stack([X, flag])
        nq = X.shape[1]
        Pb = float(cr.P_best)
        X1 = np.column_stack([X, np.sin(2 * np.pi * t / Pb), np.cos(2 * np.pi * t / Pb)])
        sg = np.r_[0.0, np.geomspace(0.1, 30, 30) * np.median(err)]
        _, s21, b = _profile_lnl(X1, tau, err ** 2, sg ** 2)
        nuis = X1[:, :nq] @ b[:nq]
        y = (tau - nuis) * DAY
        tt = np.linspace(t.min() - 200, t.max() + 200, 800)
        orb = (b[nq] * np.sin(2 * np.pi * tt / Pb) + b[nq + 1] * np.cos(2 * np.pi * tt / Pb)) * DAY
        model_pts = (X1[:, nq:] @ b[nq:]) * DAY
        e_tot = np.sqrt(err ** 2 + s21) * DAY
        for fl, col, lab in [(1, "C1", "MACHO B"), (0, "C0", "OGLE I")]:
            mm = flag == fl
            if mm.any():
                axD.errorbar(t[mm], y[mm], err[mm] * DAY, fmt="o", ms=4, color=col, label=lab)
                axDr.errorbar(t[mm], y[mm] - model_pts[mm], e_tot[mm], fmt="o", ms=3, color=col)
        axD.plot(tt, orb, "k-", lw=1.2, label=f"circular orbit, P = {Pb:.0f} d")
        axD.set(ylabel="O−C − (quadratic + offset) [s]", title="D. light-travel-time signal (season delays)")
        axD.legend(fontsize=8)
        axDr.axhline(0, c="k", lw=0.8)
        axDr.set(xlabel="HJD − 2450000", ylabel="resid. [s]")
        axDr.text(0.01, 0.8, f"error bars include fitted jitter {np.sqrt(s21) * DAY:.0f} s", transform=axDr.transAxes, fontsize=8)

        # C: the orbit signal in the light curve. Fold OGLE I with the linear ephemeris + the fitted period-change part
        # (quadratic; the MACHO offset does not apply to OGLE) only, and compare the seasons with the largest and smallest
        # ORBIT delay. The remaining shift between the two groups is the fitted light-travel-time signal (about 2A).
        oO = np.flatnonzero(flag == 0)
        orb_s = (X1[oO, nq:] @ b[nq:])                 # orbit delay per OGLE season [d]
        nuis_s = (X1[oO, :3] @ b[:3])                  # quadratic part per OGLE season [d]
        seasO = year_labels(t[oO])
        k = max(2, len(oO) // 4)
        hi, lo = np.argsort(orb_s)[-k:], np.argsort(orb_s)[:k]
        labO = year_labels(tO)
        nuis_ep = np.full(tO.size, np.nan)
        for sname, v in zip(seasO, nuis_s):
            nuis_ep[labO == sname] = v
        x = np.linspace(0, 1, 2000)
        x_rise = x[np.argmin(np.gradient(fourier_eval(fO.coef, x)))]
        zpO = dict(zip(fO.seg_names, fO.zp))
        mzp = mO - np.array([zpO.get(sg_, 0.0) for sg_ in segO])
        for idx, col, name in [(hi, "C3", "seasons of max orbit delay"), (lo, "C0", "seasons of min orbit delay")]:
            mm = np.isin(labO, seasO[idx]) & np.isfinite(nuis_ep)
            ph = (tO[mm] - nuis_ep[mm] - T0) / P
            dphi = np.mod(ph - x_rise + 0.5, 1) - 0.5
            axC.plot(dphi * P * 1440, mzp[mm], ".", ms=3, color=col, alpha=0.6, label=f"{name} ({orb_s[idx].mean() * DAY:+.0f} s)")
            xx = np.linspace(-0.25, 0.25, 600)
            axC.plot(xx * P * 1440, fourier_eval(fO.coef, x_rise + xx - orb_s[idx].mean() / P), "-", color=col, lw=1.3)
        dtau = (orb_s[hi].mean() - orb_s[lo].mean()) * DAY
        w = max(0.06, 3 * abs(dtau) / DAY / P)
        axC.set(xlim=(-w * P * 1440, w * P * 1440), xlabel="time from rising branch (period change removed) [min]",
                ylabel="I − zero point [mag]", title=f"C. rising branch, orbit-extreme seasons (Δ = {dtau:.0f} s = {dtau / 60:.1f} min)")
        axC.invert_yaxis()
        axC.legend(fontsize=7)

        # E: period search
        pg = period_grid(np.ptp(t))
        o = oc_search(t, tau, err, pg, X_extra=flag[:, None] if flag.any() else None)
        axE.plot(pg, o["Dp"], "k-", lw=1)
        axE.axhline(40, c="C3", ls="--")
        axE.axvline(Pb, c="C3", lw=0.8)
        axE.set(xscale="log", xlabel="trial orbital period [d]", ylabel="D(P)", title="E. orbit-period search")

        # F: amplitude factors
        for f_, col, lab in [(fO, "C0", "OGLE I"), (fM, "C1", "MACHO B")]:
            if f_ is not None:
                axF.errorbar(f_.t_season, f_.alpha, f_.alpha_err, fmt="o", ms=3, color=col, label=lab)
        axF.axhline(1, c="k", lw=0.8)
        axF.set(xlabel="HJD − 2450000", ylabel="α_j", title=f"F. season amplitude (χ²_ν = {cr.alpha_chi2nu:.2f})")
        axF.legend(fontsize=7)

        fig.suptitle(f"#{rank} {oid}   P_puls = {P:.5f} d, I = {par.loc[oid, 'I']}, A_I = {par.loc[oid, 'amp_I']}   |   "
                     f"P_orb = {Pb:.0f} d ({Pb / 365.25:.1f} yr), a1 sin i/c = {cr.amp_s:.0f} s, f(M) = {cr.fM:.3f} Msun, "
                     f"M2,min = {cr.M2min:.2f} Msun, K1 = {cr.K1_kms:.1f} km/s   |   D = {cr.D:.0f}, "
                     f"pred. score = {cr.pred_score if pd.notna(cr.pred_score) else 'n/a'}", fontsize=11)
        fig.savefig(out / f"{rank:02d}_{oid}.png", dpi=90)
        plt.close(fig)
        print(f"#{rank} {oid} done", flush=True)


if __name__ == "__main__":
    main()
