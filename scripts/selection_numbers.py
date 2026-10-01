"""Numbers behind the candidate selection (for the write-up): null thresholds of D, per-cut pass fractions for real stars
and each simulated class, LTTE completeness in (P_orb, M2) bins, and the completeness-corrected upper limit.

Same cuts as scripts/plot_summary_stats.py (cut_flags). Writes results/real/selection_numbers.json and prints a summary.

Usage
-----
    python scripts/selection_numbers.py --src results/real/oc_all_cm.parquet --sims results/inject/macho_v2.parquet
    python scripts/selection_numbers.py --src results/real/stats_v3.parquet --sims results/inject/stats_v3_macho.parquet --sim-amp amp_circ_s --out results/real/selection_numbers_v3.json
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_summary_stats import cut_flags  # noqa: E402

STEPS = ["c1_D", "c2_alpha", "c3_snr", "c4_cycles", "c5_ceiling"]


def funnel(x):
    m = np.ones(len(x), bool)
    out = []
    for st in STEPS:
        m &= x[st].to_numpy()
        out.append(float(m.mean()))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", default="results/real/oc_all_cm.parquet")
    ap.add_argument("--sims", nargs="+", default=["results/inject/macho_v2.parquet"])
    ap.add_argument("--out", default="results/real/selection_numbers.json")
    ap.add_argument("--sim-amp", default="amp_best_s", help="fitted-amplitude column of the sims (v3: amp_circ_s)")
    a = ap.parse_args()
    r = pd.read_parquet(a.src)
    r = r[r.ok].drop(columns=["t", "tau", "err", "flag", "alpha", "alpha_err"], errors="ignore")
    r = r.join(cut_flags(r, "amp_s", r.baseline))
    s = pd.concat([pd.read_parquet(f) for f in a.sims], ignore_index=True)
    s = s[~s.kind.isin(["error", "real"])].drop(columns=["t", "tau", "err", "flag"], errors="ignore")
    s = s.join(cut_flags(s, a.sim_amp, s.baseline if "baseline" in s else 8350.0))
    res = {}
    null = s[s.kind == "null"]
    res["null_D_quantiles"] = {q: float(np.quantile(null.D, 1 - q)) for q in ["0.05", "0.01", "0.001"] for q in [q]} \
        if False else {str(q): float(np.quantile(null.D, 1 - q)) for q in [0.05, 0.01, 0.001]}
    res["null_D_max"] = float(null.D.max())
    res["n_null"] = int(len(null))
    rm = r[r.has_M]
    res["funnel_real_all"] = dict(zip(STEPS, [int(round(f * len(r))) for f in funnel(r)]))
    res["funnel_real_macho"] = dict(zip(STEPS, [int(round(f * len(rm))) for f in funnel(rm)]))
    res["funnel_sim_fraction"] = {k: dict(zip(STEPS, funnel(g))) for k, g in s.groupby("kind")}
    res["alpha_median"] = {k: float(g.alpha_chi2nu.median()) for k, g in s.groupby("kind")}
    res["alpha_median"]["real"] = float(r.alpha_chi2nu.median())
    L = s[s.kind == "ltte"]
    comp = {}
    for p0, p1 in [(400, 1000), (1000, 3000), (3000, 10000), (1000, 10000)]:
        for m0, m1 in [(0.05, 0.15), (0.15, 0.4), (0.4, 1.5)]:
            g = L[(L.P_orb >= p0) & (L.P_orb < p1) & (L.M2 >= m0) & (L.M2 < m1)]
            comp[f"P{p0}-{p1}_M{m0}-{m1}"] = [float(g["all"].mean()), int(len(g))]
    res["completeness"] = comp
    N, k = len(rm), int(rm["all"].sum())
    up = 0.5 * chi2.ppf(0.95, 2 * (k + 1))
    res["upper_limit"] = dict(N=N, k=k, k95=float(up),
                              f_M04_15_P1_10=float(up / (N * comp["P1000-10000_M0.4-1.5"][0])),
                              f_M015_04_P1_10=float(up / (N * comp["P1000-10000_M0.15-0.4"][0])))
    Path(a.out).write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
