"""Expected decisiveness of new OGLE seasons for the frozen predictions (results/predictions/cov_<date>/).

New season delays at season centres (2017.4-2020.2: 3 seasons; 2022.8-2026.4: 4 seasons; grid-interpolated), each with
measurement + white jitter variance (sigma_season^2 + s_white^2). Under each hypothesis the predictive distribution is
Gaussian (mean, cov from the red-noise kriging). Expected log Bayes factor if H1 is true = KL(N1 || N0); if H0 is true
= -KL(N0 || N1). |ln BF| > 3 ~ strong.

Usage: python scripts/predictions_decisiveness.py [--date 2026-10-01] [--sigma 130]
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def kl(m1, C1, m0, C0):
    C0i = np.linalg.inv(C0)
    d = m0 - m1
    k = len(m1)
    return 0.5 * (np.trace(C0i @ C1) + d @ C0i @ d - k + np.linalg.slogdet(C0)[1] - np.linalg.slogdet(C1)[1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--date", default="2026-10-01")
    ap.add_argument("--sigma", type=float, default=130.0)
    a = ap.parse_args()
    seasons = {"2017-2020": np.array([8100.0, 8465.0, 8830.0]), "2017-2026": np.array([8100.0, 8465.0, 8830.0, 9950.0, 10315.0, 10680.0, 11045.0])}
    T = pd.read_csv("results/partC/partC_tiers.csv").set_index("ogle_id")
    rows = []
    for f in sorted(Path(f"results/predictions/cov_{a.date}").glob("cov_*.npz")):
        z = np.load(f)
        oid = f.stem[4:]
        tp = z["t_pred"]
        r = dict(ogle_id=oid, tier=int(T.loc[oid, "tier"]))
        for lab, ts in seasons.items():
            idx = np.array([np.argmin(np.abs(tp - x)) for x in ts])
            m1, C1 = z["h1_mean"][idx], z["h1_cov"][np.ix_(idx, idx)] + np.eye(idx.size) * (a.sigma ** 2 + float(z["white_h1_s"]) ** 2)
            m0, C0 = z["h0_mean"][idx], z["h0_cov"][np.ix_(idx, idx)] + np.eye(idx.size) * (a.sigma ** 2 + float(z["white_h0_s"]) ** 2)
            r[f"lnBF_if_orbit_{lab}"] = kl(m1, C1, m0, C0)
            r[f"lnBF_if_noorbit_{lab}"] = -kl(m0, C0, m1, C1)
        rows.append(r)
    R = pd.DataFrame(rows).sort_values(["tier", "lnBF_if_orbit_2017-2026"], ascending=[True, False])
    R.to_csv(f"results/predictions/decisiveness_{a.date}.csv", index=False)
    print(R.round(1).to_string(index=False))
    for lab in seasons:
        for tr, g in R.groupby("tier"):
            print(f"{lab} tier {tr}: median lnBF if orbit {g[f'lnBF_if_orbit_{lab}'].median():.1f}, if no orbit {g[f'lnBF_if_noorbit_{lab}'].median():.1f}; "
                  f"|lnBF|>3 if orbit: {int((g[f'lnBF_if_orbit_{lab}'] > 3).sum())}/{len(g)}, if no orbit: {int((g[f'lnBF_if_noorbit_{lab}'] < -3).sum())}/{len(g)}")


if __name__ == "__main__":
    main()
