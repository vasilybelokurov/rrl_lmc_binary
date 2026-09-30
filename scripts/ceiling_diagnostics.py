"""Why do many stars have a fitted O-C 'orbit' amplitude above the physical LTTE ceiling?

Ceiling: a1/c for M1 = 0.65, M2 = 2 Msun, edge-on, at the fitted period. For each real star compute diagnostics and
assign the dominant cause:
  slip      : a jump between consecutive seasons within one survey with |dtau| > 0.35 P (cycle slip / mis-wrap)
  boundary  : the largest jump is at the MACHO -> OGLE switch (absorbed by the free offset only if constant)
  short_P   : P_best < 1000 d (ceiling < ~500 s; the sinusoid fits noise/jitter)
  large_pc  : smooth O-C with total range > P/4 over the baseline (large period change, beyond a quadratic)
  other
and compare with the simulated classes (fraction above the ceiling per class; sims processed identically).

Usage
-----
    python scripts/ceiling_diagnostics.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from rrlbin.ltte import a1sini_over_c  # noqa: E402

DAY = 86400.0


def main():
    d = pd.read_parquet("results/real/oc_all.parquet")
    d = d[d.ok].copy()
    par = pd.read_fwf("data/raw/ogle4_lmc_rrlyr/RRab.dat", colspecs=[(0, 20), (37, 46)], names=["ogle_id", "P"], header=None)
    d = d.merge(par, on="ogle_id")
    d["ceil_s"] = a1sini_over_c(d.P_best.to_numpy(), 0.65, 2.0)
    d["above"] = d.amp_s > d.ceil_s
    cause, jump_max, rng = [], [], []
    for _, r in d.iterrows():
        t, tau, flag = map(np.asarray, (r.t, r.tau, r.flag))
        o = np.argsort(t)
        t, tau, flag = t[o], tau[o], flag[o]
        dt = np.diff(tau)
        same = flag[1:] == flag[:-1]
        jm = np.max(np.abs(dt[same])) / r.P if same.any() else 0.0
        bnd = np.max(np.abs(dt[~same])) / r.P if (~same).any() else 0.0
        rg = np.ptp(tau[flag == 0]) / r.P
        jump_max.append(jm)
        rng.append(rg)
        if not r.above:
            cause.append("-")
        elif jm > 0.35:
            cause.append("slip")
        elif bnd > 0.35 and bnd > jm:
            cause.append("boundary")
        elif r.P_best < 1000:
            cause.append("short_P")
        elif rg > 0.25:
            cause.append("large_pc")
        else:
            cause.append("other")
    d["cause"], d["jump_max_P"], d["range_P"] = cause, jump_max, rng
    a = d[d.above]
    print(f"stars {len(d)}; above ceiling {len(a)} ({len(a) / len(d):.1%})")
    print("\ncause breakdown (above ceiling):")
    for k, s in a.groupby("cause"):
        print(f"  {k:9s} N={len(s):5d} ({len(s) / len(a):.1%}); median D {s.D.median():.1f}, median P_best {s.P_best.median():.0f} d, "
              f"median amp {s.amp_s.median():.0f} s vs ceiling {s.ceil_s.median():.0f} s, median jitter {s.jit0_s.median():.0f} s, "
              f"has_M {s.has_M.mean():.2f}")
    print("\nfraction above ceiling vs P_best (real):")
    for lo, hi in [(300, 600), (600, 1000), (1000, 2000), (2000, 4000), (4000, 2e4)]:
        s = d[(d.P_best >= lo) & (d.P_best < hi)]
        print(f"  [{lo},{hi}) d: N={len(s)}, above {s.above.mean():.2f}, median ceiling {s.ceil_s.median():.0f} s")
    sims = pd.concat([pd.read_parquet(f) for f in ["results/inject/macho_run1.parquet", "results/inject/macho_pred_run1.parquet"]])
    sims = sims[sims.kind != "error"].copy()
    sims["above"] = sims.amp_best_s > a1sini_over_c(sims.P_best.to_numpy(), 0.65, 2.0)
    print("\nsimulations processed identically, fraction above ceiling:",
          sims.groupby("kind").above.mean().round(3).to_dict())
    L = sims[sims.kind == "ltte"]
    print("  LTTE sims above ceiling with injected amplitude < ceiling (estimator overshoot):",
          f"{np.mean(L.above & (L.amp_s < a1sini_over_c(L.P_orb.to_numpy(), 0.65, 2.0))):.3f}")
    d.drop(columns=["t", "tau", "err", "flag", "alpha", "alpha_err"]).to_parquet("results/real/ceiling_diag.parquet")


if __name__ == "__main__":
    main()
