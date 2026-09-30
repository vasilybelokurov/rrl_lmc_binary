"""How much does the irregular season sampling suppress annual aliases?

(1) Spectral window |W(f)| = |sum_j w_j exp(2 pi i f t_j)| / sum w_j of real season epochs (Fisher-weighted), and the
    scatter of season epochs about a strict one-per-year grid.
(2) Injection: circular LTTE with P_orb in {500, 700, 1000, 3000} d, amplitude A, on real season epochs and errors of random
    MACHO+OGLE stars; O-C search on a grid starting at 300 d; fraction with P_best within 15% of the truth vs at the
    annual alias |1/P - k/yr|^-1.
(3) Within-season smearing: the season delay is a weighted mean over a ~240-d season, which suppresses short periods.

Usage
-----
    python scripts/alias_test.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from rrlbin.ltte import oc_search, period_grid  # noqa: E402

DAY = 86400.0
YR = 365.25


def main():
    rng = np.random.default_rng(0)
    d = pd.read_parquet("results/real/oc_all_cm.parquet", columns=["ogle_id", "has_M", "t", "err", "flag"])
    d = d[d.has_M].sample(200, random_state=1)
    # (1) epoch scatter and spectral window
    dev = []
    for t in d.t:
        t = np.asarray(t)
        dev += list(np.mod(t - 245, YR) - YR / 2)
    dev = np.array(dev)
    print(f"season epochs: position within the year (days from mid-year) p10/50/90 = {np.percentile(dev, [10, 50, 90]).round(0)}; "
          f"rms {dev.std():.0f} d")
    f = np.linspace(1e-5, 1 / 300, 4000)
    W = []
    for t, e in zip(d.t, d.err):
        t, w = np.asarray(t), np.asarray(e) ** -2
        W.append(np.abs(np.exp(2j * np.pi * np.outer(f, t)) @ w) / w.sum())
    W = np.median(W, axis=0)
    k = np.argmin(np.abs(f - 1 / YR))
    print(f"median spectral window at f = 1/yr: |W| = {W[k]:.2f} (1 = perfect alias; 0 = no alias)")
    # (2) injection with a grid from 300 d
    print("\ninjection, grid from 300 d (fraction of P_best: true / alias / other):")
    for P in [500, 700, 1000, 3000]:
        alias = [1 / abs(1 / P - 1 / YR), 1 / (1 / P + 1 / YR)]
        for A_over_err in [2, 4]:
            res = []
            for t, e, fl in zip(d.t, d.err, d.flag):
                t, e, fl = map(np.asarray, (t, e, fl))
                amp = A_over_err * np.median(e)
                y = amp * np.sin(2 * np.pi * t / P + rng.uniform(0, 6.3)) + rng.normal(0, 1, t.size) * e
                r = oc_search(t, y, e, period_grid(np.ptp(t), p_min=300.0), X_extra=fl[:, None] if fl.any() else None)
                pb = r["P_best"]
                res.append("true" if abs(pb / P - 1) < 0.15 else ("alias" if min(abs(pb / a - 1) for a in alias) < 0.15 else "other"))
            res = pd.Series(res).value_counts(normalize=True)
            print(f"  P = {P:5d} d (alias {alias[0]:.0f} d), A = {A_over_err} sigma: true {res.get('true', 0):.2f}, "
                  f"alias {res.get('alias', 0):.2f}, other {res.get('other', 0):.2f}")
    # (3) smearing of a season mean over a 240-d uniform season
    print("\nwithin-season smearing (amplitude kept by a 240-d season mean, |sinc|):")
    for P in [300, 500, 700, 1000, 2000]:
        x = np.pi * 240 / P
        print(f"  P = {P} d: {abs(np.sin(x) / x):.2f}")


if __name__ == "__main__":
    main()
