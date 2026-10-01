"""Forecast: how informative would Gaia DR4 (expected Dec 2026; ~66 months, 2014.6-2020.1; ASSUMED ~2x the DR3 epochs) be
for testing the candidates' Keplerian predictions? Two bins (2014.6-2017.4, 2017.4-2020.1); per-bin delay error scaled from
the measured DR3 half-window errors: sigma_bin = sigma_DR3half * sqrt(n_DR3half / n_DR4bin) * INFLATE, with n_DR4bin ~ n_DR3
(all DR3 epochs) per bin, i.e. sigma_bin ~ sigma_DR3half / sqrt(2) * 1.64. The constant G-I lag cancels in the difference.
Also: the second bin lies entirely beyond the public OGLE-IV end (2016.3), so it is a genuine out-of-sample test.

Usage: python scripts/gaia_dr4_forecast.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from rrlbin.ltte import ltte_delay  # noqa: E402

DAY = 86400.0
INFLATE = 1.64
g = pd.read_csv("results/partC/gaia_timing.csv")
T = pd.read_csv("results/partC/partC_tiers.csv").merge(g[["ogle_id", "err1_s", "err2_s"]], on="ogle_id", how="inner")
T = T[np.isfinite(T.P_kep) & ~T.e_at_bound]
t1, t2 = 7240.0 + 0 * 0, 8380.0                      # bin centres (HJD') of 2014.6-2017.4 and 2017.4-2020.1
sig_half = np.sqrt((T.err1_s ** 2 + T.err2_s ** 2) / 2)
sig_bin = sig_half / np.sqrt(2) * INFLATE
sig_diff = np.sqrt(2) * sig_bin
pred = np.array([ltte_delay(np.array([t1, t2]), r.P_kep, r.A_kep_s / DAY, r.e_kep, r.omega_kep, r.tp_kep) for r in T.itertuples()]) * DAY
T["dtau_pred_s"] = pred[:, 1] - pred[:, 0]
T["sigma_diff_s"] = sig_diff
T["snr"] = np.abs(T.dtau_pred_s) / T.sigma_diff_s
T.to_csv("results/partC/gaia_dr4_forecast.csv", index=False)
for tr, x in T.groupby("tier"):
    print(f"tier {tr}: N={len(x)}; DR4 predicted |dtau|/sigma > 2: {int((x.snr > 2).sum())}, > 3: {int((x.snr > 3).sum())}; "
          f"median snr {x.snr.median():.2f}")
print(f"median per-bin error {np.median(sig_bin):.0f} s (DR3 half-window median {np.median(sig_half):.0f} s, uninflated)")
