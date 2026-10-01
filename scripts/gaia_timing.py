"""Gaia DR3 epoch-photometry timing test of the Part-C candidates.

For each candidate with Gaia G transits: own G-band Fourier template (K = 4) fitted jointly with two delays, one per half of
the DR3 window (split at HJD' 7470), using the OGLE-IV ephemeris (P, T0). The constant G-vs-I band lag cancels in the
DIFFERENCE dtau = tau(2nd half) - tau(1st half), which is compared with the Keplerian prediction from the OGLE+MACHO fit
(results/partC/kepler.csv) at the same Fisher-weighted epochs. Times: Gaia BJD(TCB) - 2455197.5 -> HJD' (BJD-HJD and
TCB-UTC differences are constant to ~1 min over 3 yr and cancel in dtau to < 1 s).

z = (dtau_obs - dtau_pred) / sigma. Report: Gaia timing precision; per-candidate z; fraction with |z| < 2.

Usage: python scripts/gaia_timing.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import sqlutilpy

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from rrlbin.ltte import ltte_delay  # noqa: E402
from rrlbin.timing import delays_fixed_template, fit_timing, year_labels  # noqa: E402
from rrlbin.pipeline import load_star  # noqa: E402

DAY = 86400.0
SPLIT = 7470.0
INFLATE = 1.64   # Gaia error inflation from scripts/gaia_timing_control.py (robust sd of z for 228 quiet non-candidates,
                 # Gaia-own template, OGLE reference errors included; all excess attributed to Gaia = conservative)
P_TAIL3 = 0.136  # fraction of control stars with |z| > 3 (heavy tail) -> expected chance rate of |z| > 3


def main():
    T = pd.read_csv("results/partC/partC_tiers.csv")
    T = T[np.isfinite(T.P_kep) & ~T.e_at_bound]
    g = pd.read_parquet("../rrl_lmc_rotation/data/ogle4_lmc_rrl_gaiadr3.parquet", columns=["ogle_id", "source_id"])
    par = pd.read_fwf("data/raw/ogle4_lmc_rrlyr/RRab.dat", colspecs=[(0, 20), (37, 46), (58, 68)], names=["ogle_id", "P", "T0"], header=None)
    T = T.merge(g, on="ogle_id").merge(par, on="ogle_id")
    T = T[T.source_id > 0]
    q = """select e.source_id, e.g_transit_time, e.g_transit_mag, e.g_transit_flux_over_error, e.variability_flag_g_reject
           from mytab t join gaia_dr3.epoch_photometry e on e.source_id = t.sid"""
    res = sqlutilpy.local_join(q, "mytab", (T.source_id.astype(np.int64).to_numpy(),), ("sid",), asDict=True)
    ep = {int(s): (np.asarray(a, float), np.asarray(b, float), np.asarray(c, float), np.asarray(d, bool))
          for s, a, b, c, d in zip(res["source_id"], res["g_transit_time"], res["g_transit_mag"],
                                    res["g_transit_flux_over_error"], res["variability_flag_g_reject"])}
    rows = []
    for r in T.itertuples():
        if int(r.source_id) not in ep:
            continue
        tt, mm, snr, rej = ep[int(r.source_id)]
        ok = np.isfinite(tt) & np.isfinite(mm) & np.isfinite(snr) & (snr > 0) & ~rej
        t = tt[ok] + 2455197.5 - 2450000.0
        m, e = mm[ok], 1.0857 / snr[ok]
        lab = (t > SPLIT).astype(int)
        if min((lab == 0).sum(), (lab == 1).sum()) < 10:
            continue
        try:
            f = fit_timing(t, m, e, np.full(t.size, "G"), float(r.P), float(r.T0), K=4, labels=lab, min_season=10, clip=5.0)
        except Exception:
            continue
        if f.season.size < 2:
            continue
        d_obs = (f.tau[1] - f.tau[0] + r.P / 2) % r.P - r.P / 2
        s_obs = np.hypot(*f.tau_err[:2]) * INFLATE
        pred = ltte_delay(f.t_season[:2], r.P_kep, r.A_kep_s / DAY, r.e_kep, r.omega_kep, r.tp_kep)
        d_pred = pred[1] - pred[0]
        rows.append(dict(ogle_id=r.ogle_id, tier=r.tier, nG=int(t.size), err1_s=f.tau_err[0] * DAY, err2_s=f.tau_err[1] * DAY,
                         dtau_obs_s=d_obs * DAY, dtau_pred_s=d_pred * DAY, sigma_s=s_obs * DAY, z=(d_obs - d_pred) / s_obs,
                         snr_pred=abs(d_pred) / s_obs, t1=f.t_season[0], t2=f.t_season[1]))
    R = pd.DataFrame(rows)
    R.to_csv("results/partC/gaia_timing.csv", index=False)
    print(f"candidates tested: {len(R)}; Gaia half-window delay error p10/50/90 = {np.percentile(np.r_[R.err1_s, R.err2_s], [10, 50, 90]).round(0)} s")
    print(f"predicted |dtau| / sigma p10/50/90 = {np.percentile(R.snr_pred, [10, 50, 90]).round(2)}  (a test needs this >~ 2)")
    for tr, g_ in R.groupby("tier"):
        inf = g_[g_.snr_pred > 2]
        print(f"tier {tr}: N={len(g_)}, informative (|pred|/sigma > 2): {len(inf)}; of these |z| < 2: {int((inf.z.abs() < 2).sum())}; "
              f"sign agrees: {int((np.sign(inf.dtau_obs_s) == np.sign(inf.dtau_pred_s)).sum())}")
    inf = R[R.snr_pred > 2]
    print(f"informative total {len(inf)}; |z| > 3: {int((inf.z.abs() > 3).sum())} (expected by chance from the control tail ~{P_TAIL3 * len(inf):.1f})")
    print(R[R.snr_pred > 2].sort_values("snr_pred", ascending=False)[["ogle_id", "tier", "nG", "dtau_obs_s", "dtau_pred_s", "sigma_s", "z"]]
          .round(1).to_string(index=False))


if __name__ == "__main__":
    main()
