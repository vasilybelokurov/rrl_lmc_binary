"""Mock tests of the population mixture fit (rrlbin.mixture) using the MACHO+OGLE injection runs.

Simulations are split by STAR into a density half (class densities p_k) and a mock half (draw mock 'real' samples
with known class fractions). Tests:
  A. recovery of f_LTTE in {0, 0.01, 0.03, 0.1} at a fixed nuisance mix, mock size N = 3000 (bias, bootstrap coverage);
  B. robustness: the mock nuisance mix differs from anything assumed (all fractions are free anyway) and the mock random-walk
     class is restricted to larger amplitudes (misspecified nuisance shape) -> bias of f_LTTE.

Usage
-----
    python scripts/mixture_mock.py --runs results/inject/macho_run1.parquet results/inject/macho_pred_run1.parquet
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from rrlbin.mixture import ClassDensities, bootstrap_fractions, cv_bandwidth, features, fit_fractions  # noqa: E402

CLASSES = ["null", "blazhko", "jump", "rwalk", "jump_big", "rwalk_big", "ltte"]


def draw_mock(pool: dict, fr: dict, N: int, rng) -> np.ndarray:
    n = rng.multinomial(N, [fr[k] for k in CLASSES])
    return np.vstack([pool[k][rng.integers(0, len(pool[k]), m)] for k, m in zip(CLASSES, n) if m > 0])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", nargs="+", default=["results/inject/macho_big.parquet", "results/inject/macho_tail.parquet",
                                                  "results/inject/macho_run1.parquet", "results/inject/macho_pred_run1.parquet"])
    ap.add_argument("--N", type=int, default=3000)
    ap.add_argument("--n-mock", type=int, default=20)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    d = pd.concat([pd.read_parquet(f) for f in a.runs], ignore_index=True)
    d = d[d.kind.isin(CLASSES)].dropna(subset=["D", "P_best", "alpha_chi2nu", "jit0_s", "err_med_s"])
    stars = d.ogle_id.unique()
    half = set(rng.choice(stars, stars.size // 2, replace=False))
    dens_df, mock_df = d[d.ogle_id.isin(half)], d[~d.ogle_id.isin(half)]
    Xall = features(dens_df)
    mu, sd = Xall.mean(0), Xall.std(0)
    bw = {k: cv_bandwidth(features(dens_df[dens_df.kind == k]), dens_df[dens_df.kind == k].ogle_id.to_numpy(), mu, sd)
          for k in CLASSES}
    print("cross-validated bandwidths (held-out stars):", bw)
    dens = ClassDensities({k: features(dens_df[dens_df.kind == k]) for k in CLASSES}, bw=bw)
    pool = {k: features(mock_df[mock_df.kind == k]) for k in CLASSES}
    print("density sims per class:", dens_df.kind.value_counts().to_dict())

    base = dict(null=0.45, blazhko=0.15, jump=0.08, rwalk=0.17, jump_big=0.05, rwalk_big=0.10)
    print(f"\nA. recovery, N = {a.N}, nuisance mix {base}")
    rows = []
    for fl in [0.0, 0.01, 0.03, 0.10]:
        fr = {k: v * (1 - fl) for k, v in base.items()}
        fr["ltte"] = fl
        est, cov = [], []
        for j in range(a.n_mock):
            X = draw_mock(pool, fr, a.N, rng)
            lp = dens.logpdf(X)
            f = fit_fractions(lp)[-1]
            bs = bootstrap_fractions(lp, n_boot=40, seed=j)[:, -1]
            lo, hi = np.percentile(bs, [16, 84])
            est.append(f)
            cov.append(lo <= fl <= hi)
        est = np.array(est)
        rows.append(dict(f_true=fl, f_mean=est.mean(), f_sd=est.std(), coverage68=np.mean(cov)))
        print(f"  f_LTTE true {fl:.3f}: fitted mean {est.mean():.4f} +- {est.std():.4f} (sd over mocks); 68% coverage {np.mean(cov):.2f}")

    print("\nB. misspecified mock: random-walk class restricted to rw_rms > 400 s; different nuisance mix")
    rw_sel = features(mock_df[(mock_df.kind == "rwalk") & (mock_df.rw_rms_s > 400)])
    pool_b = dict(pool, rwalk=rw_sel)
    for fl in [0.0, 0.03]:
        fr = dict(null=0.3 * (1 - fl), blazhko=0.1 * (1 - fl), jump=0.1 * (1 - fl), rwalk=0.3 * (1 - fl),
                  jump_big=0.1 * (1 - fl), rwalk_big=0.1 * (1 - fl), ltte=fl)
        est = [fit_fractions(dens.logpdf(draw_mock(pool_b, fr, a.N, rng)))[-1] for _ in range(a.n_mock)]
        print(f"  f_LTTE true {fl:.3f}: fitted mean {np.mean(est):.4f} +- {np.std(est):.4f}")
    Path("results/mixture").mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv("results/mixture/mock_recovery.csv", index=False)


if __name__ == "__main__":
    main()
