"""Background-subtracted upper limit on the number of binaries among the MACHO candidates.

N_obs real stars passing all cuts except the amplitude veto (where the empirical noise class is a valid background model:
it has no amplitude modulation); background b = empirical pass rate x N_real, with a fractional systematic uncertainty
(default 25%: noise-model realism, circularity of fitting the noise on the real stars). Bayesian: flat prior on s >= 0,
Gaussian prior on b (truncated at 0), Poisson likelihood; 95% upper limit on s; converted to a fraction with the LTTE
completeness of the same cuts (M2 0.4-1.5 Msun, P_orb 1-10 kd).

Usage: python scripts/partB_bkg_limit.py [--sys 0.25]
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import poisson

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_summary_stats import cut_flags  # noqa: E402
sys.path.insert(0, __import__('os').path.dirname(__file__))
from version import CM, LAG, PARTB, PARTC, PLOTS_C, PREV, PREV_CANDS, SERIES, SFX, SIM_M, SIM_O, STATS, V  # noqa: E402,F401

ap = argparse.ArgumentParser()
ap.add_argument("--sys", type=float, default=0.25)
a = ap.parse_args()
others = ["c1_D", "c3_snr", "c4_cycles", "c5_ceiling"]
r = pd.read_parquet(STATS)
r = r[r.has_M]
r = r.join(cut_flags(r, "amp_s", r.baseline))
s = pd.read_parquet(SIM_M)
s = s[s.kind.isin(["empirical", "ltte"])].copy()
s = s.join(cut_flags(s, "amp_circ_s", s.baseline))
e = s[s.kind == "empirical"]
L = s[(s.kind == "ltte") & (s.P_orb >= 1000) & (s.P_orb < 10000) & (s.M2 >= 0.4)]
n_obs = int(r[others].all(axis=1).sum())
b0 = e[others].all(axis=1).mean() * len(r)
comp = L[others].all(axis=1).mean()
sg = np.linspace(0, 200, 2001)
bg = np.linspace(max(b0 * (1 - 5 * a.sys), 0), b0 * (1 + 5 * a.sys), 801)
wb = np.exp(-0.5 * ((bg - b0) / (a.sys * b0)) ** 2)
post = np.array([np.sum(wb * poisson.pmf(n_obs, s_ + bg)) for s_ in sg])
cdf = np.cumsum(post) / post.sum()
s95 = sg[np.searchsorted(cdf, 0.95)]
print(f"N_obs (all cuts except alpha) = {n_obs}; background b = {b0:.0f} +- {a.sys * b0:.0f}; completeness (M2 0.4-1.5, 1-10 kd) = {comp:.3f}")
print(f"95% upper limit on binaries s < {s95:.1f}  ->  f(M2 0.4-1.5, P 1-10 kd) < {s95 / (len(r) * comp):.3%}  (model-dependent)")
