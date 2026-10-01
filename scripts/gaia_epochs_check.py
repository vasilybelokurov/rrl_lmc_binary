"""Feasibility: Gaia DR3 epoch G photometry for the Part-C candidates (WSDB gaia_dr3.epoch_photometry; source_id from the
companion project's OGLE-IV x Gaia DR3 match). Reports per star the number of usable G transits, the time span and G.

Usage: python scripts/gaia_epochs_check.py
"""
import numpy as np
import pandas as pd
import sqlutilpy

c = pd.read_csv("results/partC/partC_tiers.csv")[["ogle_id", "tier"]]
g = pd.read_parquet("../rrl_lmc_rotation/data/ogle4_lmc_rrl_gaiadr3.parquet", columns=["ogle_id", "source_id", "phot_g_mean_mag"])
c = c.merge(g, on="ogle_id", how="left")
ids = c.source_id[c.source_id > 0].astype(np.int64).to_numpy()
q = """select e.source_id, e.g_transit_time, e.g_transit_mag, e.g_transit_flux_over_error, e.variability_flag_g_reject
       from mytab t join gaia_dr3.epoch_photometry e on e.source_id = t.sid"""
res = sqlutilpy.local_join(q, "mytab", (ids,), ("sid",), asDict=True)
rows = []
for sid, tt, mm, snr, rej in zip(res["source_id"], res["g_transit_time"], res["g_transit_mag"], res["g_transit_flux_over_error"],
                                  res["variability_flag_g_reject"]):
    tt, mm = np.asarray(tt, float), np.asarray(mm, float)
    ok = np.isfinite(tt) & np.isfinite(mm) & ~np.asarray(rej, bool)
    # Gaia transit times are BJD(TCB) - 2455197.5; convert to HJD' approximately (difference HJD vs BJD < 4 s; TCB-UTC ~ 69 s,
    # irrelevant for this feasibility check)
    t_hjdp = tt[ok] + 2455197.5 - 2450000.0
    rows.append(dict(source_id=sid, n_G=int(ok.sum()), t_first=t_hjdp.min() if ok.any() else np.nan,
                     t_last=t_hjdp.max() if ok.any() else np.nan, med_snr=float(np.nanmedian(np.asarray(snr, float)[ok])) if ok.any() else np.nan))
R = c.merge(pd.DataFrame(rows), on="source_id", how="left")
R.to_csv("results/partC/gaia_epochs_check.csv", index=False)
print(f"candidates {len(c)}; with Gaia epoch photometry {R.n_G.notna().sum()}")
print("n_G p10/50/90", np.nanpercentile(R.n_G, [10, 50, 90]), "| span HJD'", np.nanmin(R.t_first).round(0), "-", np.nanmax(R.t_last).round(0),
      "| median G", R.phot_g_mean_mag.median().round(2), "| median flux S/N per transit", np.nanmedian(R.med_snr).round(1))
print("epochs after the public OGLE-IV end (HJD' > 7507): see per-star counts in the csv")
