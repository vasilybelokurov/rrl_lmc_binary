"""Part C: candidate tiers from the combined table (results/partC/partC_table.csv). Provisional (Part B pending).

Tier 1 (strong): Keplerian fit not at the e bound and >= 1.5 orbital cycles at the Keplerian period, chi2_nu(Kepler, with jitter) < 1.5, alpha chi2_nu < 1.5,
                 coherence chi2_nu < 1.5, no crowding flag, predictive score > 2.
Tier 2 (good):   not at the e bound, chi2_nu(Kepler) < 2, at most one of {crowding flag, alpha 1.5-2, coherence 1.5-2},
                 predictive score > 0 or not available.
Tier 3:          the rest (e at bound, poor fit, negative predictive score, or several flags).

Usage
-----
    python scripts/partC_tiers.py
"""
import numpy as np
import pandas as pd


def tier(r):
    if r.e_at_bound or not np.isfinite(r.P_kep) or r.baseline / r.P_kep < 1.5:   # >= 1.5 cycles at the KEPLERIAN period too
        return 3
    soft = int(bool(r.flag_crowded)) + int(1.5 <= r.alpha_chi2nu < 2) + int(1.5 <= r.coh_chi2nu < 2)
    if (r.chi2nu_kep < 1.5 and r.alpha_chi2nu < 1.5 and r.coh_chi2nu < 1.5 and not bool(r.flag_crowded)
            and np.isfinite(r.pred_score) and r.pred_score > 2):
        return 1
    if r.chi2nu_kep < 2 and soft <= 1 and (not np.isfinite(r.pred_score) or r.pred_score > 0):
        return 2
    return 3


def main():
    T = pd.read_csv("results/partC/partC_table.csv")
    if "baseline" not in T:
        T = T.merge(pd.read_parquet("results/real/stats_v3.parquet", columns=["ogle_id", "baseline"]), on="ogle_id", how="left")
    T["flag_crowded"] = T.flag_crowded.fillna(False).astype(bool)
    T["tier"] = [tier(r) for r in T.itertuples()]
    T.sort_values(["tier", "snr_tot"], ascending=[True, False]).to_csv("results/partC/partC_tiers.csv", index=False)
    print(T.tier.value_counts().sort_index().to_dict())
    cols = ["ogle_id", "bands", "V", "P_kep", "A_kep_s", "e_kep", "M2min_kep", "K1_kep", "D", "pred_score", "alpha_chi2nu",
            "coh_chi2nu", "chi2nu_kep", "dv_2027_2030"]
    t1 = T[T.tier == 1].sort_values("snr_tot", ascending=False)
    print("\nTier 1:")
    print(t1[cols].round(2).to_string(index=False))


if __name__ == "__main__":
    main()
