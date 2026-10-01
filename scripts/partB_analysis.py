"""Part B analysis: thresholds, contamination and completeness from the v3 simulations, per sample (MACHO / OGLE-only), and
the expected number of false candidates in the real data.

Cuts: C1-C5 (plot_summary_stats.cut_flags; fitted amplitude = amp_circ_s for sims, amp_s for real) and a Tier-1 proxy
(C1-C5 + alpha chi2_nu < 1.5 + coherence chi2_nu < 1.5 + predictive score > 2; the Keplerian-fit and crowding criteria of
Part C are not available for sims). The 'empirical' class carries timing noise drawn from the real stars' own REML noise fits,
so N_real * pass_rate(empirical) estimates the number of false candidates expected from the real population's noise.

Output: results/partB/summary.json, results/partB/completeness.csv; prints a summary.

Usage: python scripts/partB_analysis.py
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_summary_stats import cut_flags  # noqa: E402

OUT = Path("results/partB")


def flags(d, amp):
    f = cut_flags(d, amp, d.baseline)
    ps = d["pred_score"] if "pred_score" in d else pd.Series(np.nan, index=d.index)   # absent for OGLE-only (no MACHO)
    f["tier1p"] = f["all"] & (d.alpha_chi2nu < 1.5) & (d.coh_chi2nu < 1.5) & (ps > 2)
    return f


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    real = pd.read_parquet("results/real/stats_v3.parquet")
    real = real.join(flags(real, "amp_s"))
    res = {}
    comp_rows = []
    for name, f_, has_m in [("MACHO", "results/inject/stats_v3_macho.parquet", True), ("OGLE-only", "results/inject/stats_v3_ogle.parquet", False)]:
        s = pd.read_parquet(f_)
        s = s[s.kind != "error"].drop(columns=["t", "tau", "err", "band", "alpha", "alpha_err", "coh", "coh_err"], errors="ignore")
        s = s.join(flags(s, "amp_circ_s"))
        r = real[real.has_M == has_m]
        null = s[s.kind == "null"]
        g = dict(n_real=int(len(r)), n_real_cand=int(r["all"].sum()), n_real_tier1p=int(r["tier1p"].sum()),
                 null_D_q99=float(np.quantile(null.D, 0.99)), null_D_q999=float(np.quantile(null.D, 0.999)), null_D_max=float(null.D.max()),
                 n_null=int(len(null)))
        g["pass_cuts"] = {k: float(x["all"].mean()) for k, x in s.groupby("kind")}
        g["pass_tier1p"] = {k: float(x["tier1p"].mean()) for k, x in s.groupby("kind")}
        g["n_sims"] = {k: int(len(x)) for k, x in s.groupby("kind")}
        emp = s[s.kind == "empirical"]
        for lab, col in [("cuts", "all"), ("tier1p", "tier1p")]:
            k = int(emp[col].sum())
            rate = k / len(emp)
            up = 0.5 * chi2.ppf(0.95, 2 * (k + 1)) / len(emp)
            g[f"expected_false_{lab}"] = dict(rate=rate, k=k, n=len(emp), expected=rate * len(r), expected_95up=up * len(r))
        L = s[s.kind == "ltte"]
        for p0, p1 in [(400, 1000), (1000, 3000), (3000, 10000)]:
            for m0, m1 in [(0.05, 0.15), (0.15, 0.4), (0.4, 1.5)]:
                x = L[(L.P_orb >= p0) & (L.P_orb < p1) & (L.M2 >= m0) & (L.M2 < m1)]
                comp_rows.append(dict(sample=name, P=f"{p0}-{p1}", M2=f"{m0}-{m1}", n=len(x), comp_cuts=x["all"].mean(), comp_tier1p=x["tier1p"].mean()))
        # completeness-corrected upper limit (all candidates treated as binaries), M2 0.4-1.5, P 1-10 kd
        x = L[(L.P_orb >= 1000) & (L.P_orb < 10000) & (L.M2 >= 0.4)]
        k95 = 0.5 * chi2.ppf(0.95, 2 * (g["n_real_cand"] + 1))
        g["upper_limit_M04_P1_10"] = float(k95 / (len(r) * x["all"].mean()))
        x2 = L[(L.P_orb >= 1000) & (L.P_orb < 10000) & (L.M2 >= 0.15) & (L.M2 < 0.4)]
        g["upper_limit_M015_04_P1_10"] = float(k95 / (len(r) * x2["all"].mean()))
        res[name] = g
    Path(OUT / "summary.json").write_text(json.dumps(res, indent=1))
    pd.DataFrame(comp_rows).to_csv(OUT / "completeness.csv", index=False)
    for name, g in res.items():
        print(f"\n=== {name}: real stars {g['n_real']}, candidates (cuts) {g['n_real_cand']}, Tier-1 proxy {g['n_real_tier1p']}")
        print(f"  null D q99/q99.9/max = {g['null_D_q99']:.1f}/{g['null_D_q999']:.1f}/{g['null_D_max']:.1f} (N={g['n_null']})")
        print("  pass rate, cuts:   ", {k: round(v, 4) for k, v in g["pass_cuts"].items()})
        print("  pass rate, Tier-1p:", {k: round(v, 4) for k, v in g["pass_tier1p"].items()})
        for lab in ("cuts", "tier1p"):
            e = g[f"expected_false_{lab}"]
            print(f"  expected false ({lab}) from empirical noise: {e['expected']:.1f} (95% up {e['expected_95up']:.1f}; {e['k']}/{e['n']} sims)")
        print(f"  upper limit f(M2 0.4-1.5, 1-10 kd) < {g['upper_limit_M04_P1_10']:.3%}; f(M2 0.15-0.4) < {g['upper_limit_M015_04_P1_10']:.3%}")
    print("\ncompleteness (cuts / Tier-1 proxy):")
    print(pd.DataFrame(comp_rows).round(3).to_string(index=False))


if __name__ == "__main__":
    main()
