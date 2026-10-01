"""Part C6: v3 candidate sheets (PNG, plots/partC/sheets/) and the combined Part-C table (results/partC/partC_table.csv).

Per candidate: (A) O-C (common mode removed, bands aligned) minus the fitted quadratic and band offsets, with the Keplerian
and circular LTTE curves; (B) residuals; (C) orbit-period search D(P); (D) per-season amplitude factors alpha_j and
harmonic coherence (tau_1 - tau_h) by band; header with P, A, e, f(M), M2,min, K1, V, predictive score, crowding flags.
Fits with e >= 0.9 are flagged 'e at bound (unreliable)'.

Usage
-----
    python scripts/partC_sheets.py [--n 75]
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

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from level2 import load, priors_for, to_series  # noqa: E402
from rrlbin.kepler_fit import fit_keplerian  # noqa: E402
from rrlbin.ltte import ltte_delay  # noqa: E402
from rrlbin.oc import align_bands, apply_common_mode, orbit_search, period_grid, robust_unwrap  # noqa: E402
from rrlbin.timing import year_labels  # noqa: E402

DAY = 86400.0
BC = {0: ("C0", "OGLE I"), 1: ("C1", "MACHO B"), 2: ("C2", "MACHO R")}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=200)
    a = ap.parse_args()
    out = Path("plots/partC/sheets")
    out.mkdir(parents=True, exist_ok=True)
    for f in out.glob("*.png"):
        f.unlink()
    cand = pd.read_csv("results/partC/candidates_v3_prov.csv")
    kep = pd.read_csv("results/partC/kepler.csv").set_index("ogle_id")
    crowd = pd.read_csv("results/partC/crowding.csv").set_index("ogle_id")
    fol = pd.read_csv("results/partC/followup.csv").set_index("ogle_id")
    st = pd.read_parquet("results/real/stats_v3.parquet", columns=["ogle_id", "jit1_s"]).set_index("ogle_id")
    ser = load("results/real/series_v3")
    ser = ser[ser.ogle_id.isin(set(cand.ogle_id))].set_index("ogle_id")
    cm = {tuple(map(int, k.split(","))): v / DAY for k, v in json.loads(Path("results/calib/common_mode_v3.json").read_text()).items()}
    lag = json.loads(Path("results/calib/band_lag_v3.json").read_text())
    rows = []
    for rank, c in enumerate(cand.head(a.n).itertuples(), 1):
        oid = c.ogle_id
        r = dict(ser.loc[oid].to_dict(), ogle_id=oid)
        s = apply_common_mode(to_series(r), cm, year_labels)
        t, y, e, band = s["t"], s["tau"], s["err"], s["band"].astype(int)
        al, ale, coh, cohe = (np.asarray(s[k]) for k in ("alpha", "alpha_err", "coh", "coh_err"))
        o = np.argsort(t)
        t, y, e, band, al, ale, coh, cohe = (x[o] for x in (t, y, e, band, al, ale, coh, cohe))
        pri = priors_for(s["P"], lag)
        y = robust_unwrap(t, y, band, s["P"])
        y = align_bands(t, y, e, band, s["P"], pri)
        k = kep.loc[oid] if oid in kep.index else None
        jit = float(st.loc[oid, "jit1_s"]) / DAY
        P0 = float(k.P_kep) if (k is not None and np.isfinite(k.P_kep)) else float(c.P_best)
        try:
            fk = fit_keplerian(t, y, e, band, P0, s_jit=jit, priors=pri, n_boot=0)
            fc = fit_keplerian(t, y, e, band, c.P_best, s_jit=jit, priors=pri, n_boot=0, circular=True)
        except Exception as ex:
            rows.append(dict(rank=rank, ogle_id=oid, flags=f"Keplerian fit failed: {str(ex)[:60]}"))
            continue
        lin = np.asarray(fk["beta"])[:-1]
        from rrlbin.oc import design
        X, _ = design(t, band, np.average(t, weights=np.sqrt(e ** 2 + jit ** 2) ** -2))
        resid_nuis = (y - X @ lin) * DAY
        tt = np.linspace(t.min() - 300, t.max() + 300, 900)
        kc = ltte_delay(tt, fk["P"], fk["A_s"] / DAY, fk["e"], fk["omega"], fk["t_p"]) * DAY
        cc = ltte_delay(tt, fc["P"], fc["A_s"] / DAY, 0.0, fc["omega"], fc["t_p"]) * DAY
        model_pts = ltte_delay(t, fk["P"], fk["A_s"] / DAY, fk["e"], fk["omega"], fk["t_p"]) * DAY
        fig = plt.figure(figsize=(15, 9))
        gs = fig.add_gridspec(3, 2, height_ratios=[1.2, 0.45, 1], hspace=0.45, wspace=0.22)
        aA, aB = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[1, 0])
        aC, aD, aE = fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[2, 0]), fig.add_subplot(gs[2, 1])
        stot = np.sqrt(e ** 2 + jit ** 2) * DAY
        for b in np.unique(band):
            m = band == b
            col, lab = BC[b]
            aA.errorbar(t[m], resid_nuis[m], e[m] * DAY, fmt="o", ms=4, color=col, label=lab)
            aB.errorbar(t[m], resid_nuis[m] - model_pts[m], stot[m], fmt="o", ms=3, color=col)
            aD.errorbar(t[m], al[m], ale[m], fmt="o", ms=3, color=col, label=lab)
            aE.errorbar(t[m], coh[m] * DAY, cohe[m] * DAY, fmt="o", ms=3, color=col, label=lab)
        aA.plot(tt, kc, "k-", lw=1.3, label=f"Keplerian P={fk['P']:.0f} d, e={fk['e']:.2f}")
        aA.plot(tt, cc, "k:", lw=1, label=f"circular P={fc['P']:.0f} d")
        aA.set(ylabel="O−C − (quadratic + offsets) [s]", title="A. LTTE signal (common mode removed, bands aligned)")
        aA.legend(fontsize=7, ncol=2)
        aB.axhline(0, c="k", lw=0.8)
        aB.set(xlabel="HJD − 2450000", ylabel="resid. [s]")
        pg = period_grid(np.ptp(t))
        sr = orbit_search(t, y, e, band, pg, 1, pri)
        aC.plot(pg, sr["Dp"], "k-")
        aC.axhline(40, c="C3", ls="--")
        aC.set(xscale="log", xlabel="trial P_orb [d]", ylabel="D(P)", title="C. orbit-period search")
        aD.axhline(1, c="k", lw=0.8)
        aD.set(xlabel="HJD − 2450000", ylabel="α_j", title=f"D. season amplitude (χ²_ν = {c.alpha_chi2nu:.2f})")
        aD.legend(fontsize=7)
        aE.axhline(0, c="k", lw=0.8)
        aE.set(xlabel="HJD − 2450000", ylabel="τ₁ − τ_h [s]", title=f"E. harmonic coherence (χ²_ν = {c.coh_chi2nu:.2f})")
        cr = crowd.loc[oid] if oid in crowd.index else None
        fo = fol.loc[oid] if oid in fol.index else None
        flags = []
        if k is not None and k.e_kep >= 0.9:
            flags.append("e at bound (unreliable)")
        if cr is not None and bool(cr.flag_crowded):
            flags.append("crowding flag")
        if c.alpha_chi2nu > 1.5:
            flags.append("α χ²>1.5")
        ps = "n/a" if pd.isna(c.pred_score) else f"{c.pred_score:.1f}"
        hdr = (f"#{rank} {oid}  V={fo.V if fo is not None else np.nan:.2f}  bands {c.bands}  |  P={fk['P']:.0f} d ({fk['P'] / 365.25:.1f} yr)  "
               f"A={fk['A_s']:.0f} s  e={fk['e']:.2f}  f(M)={k.fM_kep if k is not None else np.nan:.3f}  M2,min={k.M2min_kep if k is not None else np.nan:.2f}  "
               f"K1={k.K1_kep if k is not None else np.nan:.1f} km/s  |  D={c.D:.0f}  pred={ps}  ΔBIC_ecc={k.dBIC_ecc if k is not None else np.nan:.1f}"
               + ("  |  FLAGS: " + ", ".join(flags) if flags else ""))
        fig.suptitle(hdr, fontsize=10)
        fig.savefig(out / f"{rank:02d}_{oid}.png", dpi=90)
        plt.close(fig)
        rows.append(dict(rank=rank, ogle_id=oid, flags="; ".join(flags)))
    T = cand.merge(kep.reset_index(), on="ogle_id", how="left").merge(
        crowd.reset_index()[["ogle_id", "flag_crowded", "n_gaia_within", "pct_sigma_all", "ruwe", "dGI"]], on="ogle_id", how="left").merge(
        fol.reset_index()[["ogle_id", "V", "dv_2027_2030"]], on="ogle_id", how="left").merge(pd.DataFrame(rows), on="ogle_id", how="left")
    T["e_at_bound"] = T.e_kep >= 0.9
    T.to_csv("results/partC/partC_table.csv", index=False)
    print(f"sheets: {len(rows)} -> {out}; table -> results/partC/partC_table.csv; e at bound: {int(T.e_at_bound.sum())}; "
          f"crowding flag: {int(T.flag_crowded.fillna(False).sum())}")


if __name__ == "__main__":
    main()
