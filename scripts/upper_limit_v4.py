"""Upper limit on the binary fraction of LMC RRab from the 1992-2026 (v4) timing data, MACHO+OGLE stars (N = 6,612).

Counting argument (model-free in the noise): if a fraction f of the N stars has a companion drawn from the injected orbit
population of a given bin, the expected number of flagged binaries is f N eps, where eps is the probability that such a binary
star passes all candidate cuts. Every one of the k flagged stars is counted as a binary (no subtraction of noise candidates), so
f N eps <= k95 (classical Poisson 95% upper limit for k observed) -> f < k95 / (N eps).

eps, primary: orbits injected into the REAL stars' season delays (scripts/inject_real.py): real intrinsic timing noise, real
amplitude (alpha) behaviour, real cadence and errors; no separate amplitude-veto correction. eps uncertainty: binomial (68%)
from the number of injections in the bin, propagated to the limit.
eps, comparison: light-curve-level injections into simulated stars without intrinsic noise (v4 sims, class 'ltte'), times the
real alpha pass rate 0.609 (the earlier estimate), and series-level injection into the same simulated stars (validation of the
series-level shortcut).
Noise-subtracted (model-dependent) limit: Bayesian, flat prior on the signal, background b = expected noise candidates from the
empirical red-noise simulations (x the real alpha pass rate of such stars), and with b halved.

Injected population: P ~ logU(300, 1e4) d, M2 ~ logU(0.05, 1.5) Msun, M1 = 0.65, isotropic inclination, e = 0 (half) or U(0, 0.7).
Output: results/partB_v4/upper_limits.csv; prints the table.
Usage: python scripts/upper_limit_v4.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import brentq
from scipy.stats import beta as beta_dist
from scipy.stats import poisson

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_summary_stats import cut_flags  # noqa: E402

ALPHA_PASS_ALL = None   # computed below (fraction of the N stars passing the amplitude veto)
M_BINS = [("0.4-1.5", 0.4, 1.5), ("0.15-0.4", 0.15, 0.4), ("0.05-0.15", 0.05, 0.15)]
P_BINS = [(300, 1000), (1000, 3000), (3000, 10000), (1000, 10000)]


def k95_poisson(k, cl=0.95):
    """Classical upper limit on a Poisson mean given k observed: P(N <= k | mu) = 1 - cl."""
    return brentq(lambda mu: poisson.cdf(k, mu) - (1 - cl), 1e-9, 10 * k + 50)


def k95_bayes(k, b, cl=0.95):
    """Bayesian upper limit on the signal s for k observed with known background b (flat prior s >= 0)."""
    s = np.linspace(0, 10 * k + 50, 200001)
    lik = poisson.pmf(k, s + b)
    cdf = np.cumsum(lik) / lik.sum()
    return float(s[np.searchsorted(cdf, cl)])


def eff(df, m0, m1, p0, p1):
    g = df[(df.M2 >= m0) & (df.M2 < m1) & (df.P_orb >= p0) & (df.P_orb < p1)]
    n, x = len(g), int(g["all"].sum())
    lo, hi = (beta_dist.ppf(0.16, x, n - x + 1) if x > 0 else 0.0), beta_dist.ppf(0.84, x + 1, n - x)
    return n, x / max(n, 1), lo, hi


def main():
    real = pd.read_parquet("results/real/stats_v4.parquet")
    real = real.join(cut_flags(real, "amp_s", real.baseline))
    R = real[real.has_M & (real.err_msg.isna() if "err_msg" in real else True)]
    N, k = len(R), int(R["all"].sum())
    k95 = k95_poisson(k)
    a_all = float((R.alpha_chi2nu < 2).mean())
    sims = pd.read_parquet("results/inject/stats_v4_macho.parquet")
    sims = sims.join(cut_flags(sims, "amp_s", sims.baseline))
    emp = sims[sims.kind == "empirical"]
    noal = ["c1_D", "c3_snr", "c4_cycles", "c5_ceiling"]
    a_cond = float(R[R[noal].all(axis=1)].c2_alpha.mean())
    b = float(emp[noal].all(axis=1).mean() * a_cond * N)
    inj_real = pd.read_parquet("results/partB_v4/inject_real.parquet")
    inj_real = inj_real[inj_real.err_msg.isna()] if "err_msg" in inj_real else inj_real
    lc = sims[sims.kind == "ltte"]
    p_ser = Path("results/partB_v4/inject_null_sims.parquet")
    inj_ser = pd.read_parquet(p_ser) if p_ser.exists() else None
    rows = []
    for name, m0, m1 in M_BINS:
        for p0, p1 in P_BINS:
            n, e, elo, ehi = eff(inj_real, m0, m1, p0, p1)
            nl, el, _, _ = eff(lc, m0, m1, p0, p1)
            row = dict(M2=name, P=f"{p0}-{p1}", n_inj=n, eps=e, eps_lo=elo, eps_hi=ehi, eps_lc_sims=el, eps_lc_sims_x_alpha=el * a_all)
            if inj_ser is not None:
                row["eps_series_sims"] = eff(inj_ser, m0, m1, p0, p1)[1]
            if e > 0:
                row.update(f_max=k95 / (N * e), f_max_lo=k95 / (N * ehi), f_max_hi=k95 / (N * elo) if elo > 0 else np.inf,
                           f_max_bkgsub=k95_bayes(k, b) / (N * e), f_max_bkgsub_b_half=k95_bayes(k, 0.5 * b) / (N * e),
                           f_max_old_method=k95 / (N * el * a_all) if el > 0 else np.inf)
            rows.append(row)
    T = pd.DataFrame(rows)
    T.attrs = {}
    Path("results/partB_v4").mkdir(parents=True, exist_ok=True)
    T.to_csv("results/partB_v4/upper_limits.csv", index=False)
    meta = dict(N=N, k=k, k95=k95, b=b, alpha_pass_all=a_all, alpha_pass_given_other_cuts=a_cond, n_inj_real=len(inj_real))
    pd.Series(meta).to_json("results/partB_v4/upper_limits_meta.json", indent=1)
    print(f"N = {N}; k = {k}; k95 = {k95:.1f}; noise background b = {b:.0f}; alpha pass: all {a_all:.3f}, given other cuts {a_cond:.3f}; "
          f"real-star injections {len(inj_real)}")
    print(T.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
