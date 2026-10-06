"""Upper limit on the binary fraction of LMC RRab from the 1992-2026 (v4) timing data, MACHO+OGLE stars.

Logic (model-free, conservative): a binary with companion mass M2 and period P is recovered by the candidate selection with
probability eps (completeness), measured by injecting Keplerian LTTE orbits into simulated light curves of real stars (same
epochs, errors, pipeline, cuts). If a fraction f of the N stars had such companions, about f N eps would be detected. We detect
k candidates in total, and EVERY one of them is counted as a binary (no noise subtraction), so f N eps <= k95 (Poisson 95% upper
limit on k) -> f < k95 / (N eps). The noise-subtracted (model-dependent) limit uses the expected number of noise candidates b
from the empirical red-noise simulations: f < (k95(k, b)) / (N eps), Bayesian with a flat prior on the signal.

Injected population: P ~ logU(300, 1e4) d, M2 ~ logU(0.05, 1.5) Msun, M1 = 0.65, isotropic inclination, half eccentric.
Output: results/partB_v4/upper_limits.csv; prints the table.
Usage: python scripts/upper_limit_v4.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import brentq
from scipy.stats import poisson

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_summary_stats import cut_flags  # noqa: E402


def k95_poisson(k, cl=0.95):
    """Classical upper limit on a Poisson mean given k observed: P(N <= k | mu) = 1 - cl."""
    return brentq(lambda mu: poisson.cdf(k, mu) - (1 - cl), 1e-9, 10 * k + 50)


def k95_bayes(k, b, cl=0.95):
    """Bayesian upper limit on the signal s for k observed with known background b (flat prior s >= 0)."""
    s = np.linspace(0, 10 * k + 50, 200001)
    lik = poisson.pmf(k, s + b)
    cdf = np.cumsum(lik) / lik.sum()
    return float(s[np.searchsorted(cdf, cl)])


def main():
    real = pd.read_parquet("results/real/stats_v4.parquet")
    sims = pd.read_parquet("results/inject/stats_v4_macho.parquet")
    real = real.join(cut_flags(real, "amp_s", real.baseline))
    sims = sims.join(cut_flags(sims, "amp_s", sims.baseline))
    R = real[real.has_M & real.err_msg.isna()] if "err_msg" in real else real[real.has_M]
    N, k = len(R), int(R["all"].sum())
    emp = sims[sims.kind == "empirical"]
    noal = ["c1_D", "c3_snr", "c4_cycles", "c5_ceiling"]
    a_pass = R[R[noal].all(axis=1)].c2_alpha.mean()               # alpha-veto pass rate of real stars (sims carry no amplitude noise)
    b = emp[noal].all(axis=1).mean() * a_pass * N                  # expected noise candidates
    L = sims[sims.kind == "ltte"]
    k95 = k95_poisson(k)
    rows = []
    bins = [("M2 0.4-1.5", 0.4, 1.5), ("M2 0.15-0.4", 0.15, 0.4), ("M2 0.05-0.15", 0.05, 0.15)]
    pbins = [(1000, 3000), (3000, 10000), (1000, 10000), (300, 1000)]
    for name, m0, m1 in bins:
        for p0, p1 in pbins:
            g = L[(L.M2 >= m0) & (L.M2 < m1) & (L.P_orb >= p0) & (L.P_orb < p1)]
            eps = g["all"].mean()
            if len(g) < 20 or eps <= 0:
                rows.append(dict(bin=name, P=f"{p0}-{p1}", n_inj=len(g), eps=eps))
                continue
            rows.append(dict(bin=name, P=f"{p0}-{p1}", n_inj=len(g), eps=eps, N=N, k=k, k95=k95, f_max_modelfree=k95 / (N * eps),
                             b=b, f_max_bkgsub=k95_bayes(k, b) / (N * eps), f_max_bkgsub_sys50=k95_bayes(k, 0.5 * b) / (N * eps)))
    T = pd.DataFrame(rows)
    Path("results/partB_v4").mkdir(parents=True, exist_ok=True)
    T.to_csv("results/partB_v4/upper_limits.csv", index=False)
    print(f"MACHO+OGLE stars N = {N}; candidates k = {k}; Poisson 95% k95 = {k95:.1f}; expected noise candidates b = {b:.0f} "
          f"(empirical red noise, alpha pass {a_pass:.2f})")
    print(T.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
