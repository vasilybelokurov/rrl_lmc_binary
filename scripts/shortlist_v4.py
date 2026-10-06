"""Short list of the most binary-like LMC RRab for follow-up (RV / Gaia), from the v4 (1992-2026) analysis.

Starting set: the 140 v4 candidates + the 27 frozen-prediction stars (results/partC_v4/kepler.csv: Keplerian fits on all data).
Checks (each reported; a star is short-listed if it passes 1-4):
  1. physical:   Keplerian amplitude below the LTTE ceiling (M1 = 0.65, M2 = 2 Msun, edge-on) at the Keplerian period, and the fit is
                 sane (e below the 0.95 bound, chi2_nu < 2)
  2. cycles:     >= 2 orbital cycles in the 1992-2026 baseline (P_kep <= baseline / 2)
  3. amplitude:  pulsation amplitude steady (alpha chi2_nu < 2 in Level 2 AND no smooth amplitude variability in the survey,
                 d2lnl < 13.8): an orbit does not change the light curve; most intrinsic modulators do
  4. prediction: the orbit predicted later data: frozen test (orbits frozen 2026-10-01, blind) passed, OR the split test (fit to
                 1992-2016, predict 2017-2026) passed: p(chi2|orbit) > 0.01 and ln BF > 3. The split test is weaker: for stars first
                 selected with the post-2016 data it is not blind.
  bonus:         eccentric (dBIC_ecc > 6 and e_p16 > 0.1): a lopsided timing curve is hard to produce intrinsically.
Rank: blind-test passes first, then by ln BF of the split test.
Output: results/partC_v4/shortlist.csv (all checked stars with flags), printed short list.
Usage: python scripts/shortlist_v4.py
"""
import numpy as np
import pandas as pd

G, MSUN, C, DAY = 6.674e-11, 1.989e30, 2.998e8, 86400.0


def ceiling_s(P_d, M1=0.65, M2=2.0):
    a = (G * (M1 + M2) * MSUN * (P_d * DAY) ** 2 / (4 * np.pi ** 2)) ** (1 / 3)
    return a * M2 / (M1 + M2) / C


def main():
    K = pd.read_csv("results/partC_v4/kepler.csv")
    K = K[K.err_msg.isna()] if "err_msg" in K else K
    st = pd.read_parquet("results/real/stats_v4.parquet", columns=["ogle_id", "alpha_chi2nu", "baseline", "has_M", "D"]).set_index("ogle_id")
    sv = pd.read_parquet("results/stage1/survey_all.parquet").set_index("ogle_id")[["av_d2lnl"]]
    sp = pd.read_csv("results/split_v4/real.csv").set_index("ogle_id")[["ok", "lnBF", "p_h1"]].rename(columns={"lnBF": "lnBF_split", "p_h1": "p_split"})
    fz = pd.read_csv("results/predictions/test_2026-10-06.csv").set_index("ogle_id")[["lnBF_all", "p_h1_all"]]
    T = K.set_index("ogle_id").join(st).join(sv).join(sp).join(fz)
    T["phys"] = (T.A_kep_s < ceiling_s(T.P_kep)) & (T.e_kep < 0.94) & (T.chi2nu_kep < 2)
    T["cycles"] = T.baseline / T.P_kep >= 2
    T["amp_steady"] = (T.alpha_chi2nu < 2) & (T.av_d2lnl.fillna(0) < 13.8)
    T["blind_pass"] = (T.p_h1_all > 0.01) & (T.lnBF_all > 3)
    T["split_pass"] = (T.p_split > 0.01) & (T.lnBF_split > 3)
    T["predicts"] = T.blind_pass | T.split_pass
    T["eccentric"] = (T.dBIC_ecc > 6) & (T.e_p16 > 0.1)
    T["shortlist"] = T.phys & T.cycles & T.amp_steady & T.predicts
    T = T.sort_values(["blind_pass", "lnBF_split"], ascending=[False, False])
    T.to_csv("results/partC_v4/shortlist.csv")
    print(f"checked {len(T)}: physical {int(T.phys.sum())}, >= 2 cycles {int(T.cycles.sum())}, amplitude steady {int(T.amp_steady.sum())}, "
          f"predicts {int(T.predicts.sum())} (blind {int(T.blind_pass.sum())}) -> SHORT LIST {int(T.shortlist.sum())}")
    S = T[T.shortlist]
    cols = ["P_kep", "e_kep", "A_kep_s", "M2min_kep", "K1_kep", "K1_p16", "K1_p84", "blind_pass", "lnBF_split", "eccentric", "has_M"]
    print(S[cols].round(2).to_string())
    print("\nfailed only the prediction test (physical, >= 2 cycles, steady amplitude):",
          list(T[T.phys & T.cycles & T.amp_steady & ~T.predicts].index.str[-5:]))


if __name__ == "__main__":
    main()
