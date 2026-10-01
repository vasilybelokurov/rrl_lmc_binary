"""How decisive would the (non-public) OGLE-IV seasons 2017.4-2020.4 be for each candidate?

For each candidate: the Keplerian model (orbit + quadratic + offsets; results/partC/kepler.csv refit here to get the full model)
and the no-orbit model (quadratic + offsets, H0) are both fitted to the existing data and extrapolated to the centres of the
three post-2016 seasons. Discriminating power: S/N = sqrt(sum_k ((tau_orbit_k - tau_H0_k) / sigma_k)^2), sigma_k = per-season
delay error expected from the measured cadence of the candidates' fields (default 130 s; scripts/ogle_post2016_cadence.py).
The uncertainty of the extrapolations themselves (orbit and quadratic) is NOT included -> optimistic.

Output: results/partC/ogle_post2016_forecast.csv

Usage: python scripts/ogle_post2016_forecast.py [--sigma 130]
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from level2 import load, priors_for, to_series  # noqa: E402
from rrlbin.kepler_fit import fit_keplerian  # noqa: E402
from rrlbin.ltte import ltte_delay  # noqa: E402
from rrlbin.oc import align_bands, apply_common_mode, default_s_grid, design, profile_lnl, robust_unwrap  # noqa: E402
from rrlbin.timing import year_labels  # noqa: E402

DAY = 86400.0
T_POST = np.array([8100.0, 8465.0, 8830.0])   # centres of the 2017.4-18.4, 2018.4-19.4, 2019.4-20.4 seasons (HJD')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sigma", type=float, default=130.0)
    a = ap.parse_args()
    T = pd.read_csv("results/partC/partC_tiers.csv")
    T = T[np.isfinite(T.P_kep) & ~T.e_at_bound & (T.tier <= 2)]
    ser = load("results/real/series_v3").set_index("ogle_id")
    st = pd.read_parquet("results/real/stats_v3.parquet", columns=["ogle_id", "jit1_s"]).set_index("ogle_id")
    cm = {tuple(map(int, k.split(","))): v / DAY for k, v in json.loads(Path("results/calib/common_mode_v3.json").read_text()).items()}
    lag = json.loads(Path("results/calib/band_lag_v3.json").read_text())
    rows = []
    for r in T.itertuples():
        s = apply_common_mode(to_series(dict(ser.loc[r.ogle_id].to_dict(), ogle_id=r.ogle_id)), cm, year_labels)
        t, y, e, band = s["t"], s["tau"], s["err"], s["band"].astype(int)
        o = np.argsort(t)
        t, y, e, band = t[o], y[o], e[o], band[o]
        pri = priors_for(s["P"], lag)
        y = align_bands(t, robust_unwrap(t, y, band, s["P"]), e, band, s["P"], pri)
        jit = float(st.loc[r.ogle_id, "jit1_s"]) / DAY
        k = fit_keplerian(t, y, e, band, float(r.P_kep), s_jit=jit, priors=pri, n_boot=0)
        t_ref = np.average(t, weights=np.sqrt(e ** 2 + jit ** 2) ** -2)
        # Keplerian prediction (OGLE band: offsets do not apply)
        bk = np.asarray(k["beta"])
        x = (T_POST - t_ref) / 1000.0
        quad_k = bk[0] * x ** 2 + bk[1] * x + bk[2]
        orb = ltte_delay(T_POST, k["P"], k["A_s"] / DAY, k["e"], k["omega"], k["t_p"])
        pred_k = quad_k + orb
        # no-orbit (H0) prediction
        X, names = design(t, band, t_ref)
        _, _, b0, _ = profile_lnl(X, y, e, names, default_s_grid(e), pri)
        pred_0 = b0[0] * x ** 2 + b0[1] * x + b0[2]
        d = (pred_k - pred_0) * DAY
        rows.append(dict(ogle_id=r.ogle_id, tier=r.tier, P_orb_yr=k["P"] / 365.25, A_s=k["A_s"], diff_2017_s=d[0], diff_2018_s=d[1],
                         diff_2019_s=d[2], snr=float(np.sqrt(np.sum((d / a.sigma) ** 2)))))
    R = pd.DataFrame(rows).sort_values(["tier", "snr"], ascending=[True, False])
    R.to_csv("results/partC/ogle_post2016_forecast.csv", index=False)
    print(R.round(1).to_string(index=False))
    for tr, g in R.groupby("tier"):
        print(f"tier {tr}: N={len(g)}; S/N > 3: {int((g.snr > 3).sum())}, > 5: {int((g.snr > 5).sum())}; median {g.snr.median():.1f}")


if __name__ == "__main__":
    main()
