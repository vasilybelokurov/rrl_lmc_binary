"""Population mixture fit of per-star O-C summary statistics.

Each star i has a feature vector x_i (functions of the O-C search statistic D, best period, amplitude-constancy
statistic, and the fitted extra jitter). The population is modelled as a mixture of classes k (null, Blazhko, jump,
random walk, LTTE) whose feature densities p_k(x) are estimated from simulations on the real cadences:

    p(x_i | f) = sum_k f_k p_k(x_i),     sum_k f_k = 1.

The class fractions f are fitted by EM (maximum likelihood) and their uncertainty by bootstrap over stars.
The LTTE fraction is conditional on the orbital prior used in the simulations.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

FEATURES = ["lnD", "logP", "lnalpha", "jit"]


def features(df: pd.DataFrame) -> np.ndarray:
    """Feature matrix (n, 4) from a table with columns D, P_best, alpha_chi2nu, jit0_s, err_med_s."""
    return np.column_stack([
        np.log1p(np.clip(df.D.to_numpy(float), 0, None)),
        np.log10(df.P_best.to_numpy(float)),
        np.log(np.clip(df.alpha_chi2nu.to_numpy(float), 1e-3, None)),
        np.arcsinh(df.jit0_s.to_numpy(float) / df.err_med_s.to_numpy(float)),
    ])


class ClassDensities:
    """Gaussian-KDE density per class in feature space (features standardized with the pooled scale)."""

    def __init__(self, X_by_class: dict[str, np.ndarray], bw: float | str | dict = "scott", jitter: float = 0.02,
                 seed: int = 0):
        rng = np.random.default_rng(seed)
        allX = np.vstack(list(X_by_class.values()))
        self.mu, self.sd = allX.mean(0), allX.std(0)
        self.classes = list(X_by_class)
        self.kde = {}
        for k, X in X_by_class.items():
            Z = (X - self.mu) / self.sd
            Z = Z + rng.normal(0, jitter, Z.shape)   # break ties (e.g. jitter = 0 exactly) so the KDE is non-singular
            self.kde[k] = gaussian_kde(Z.T, bw_method=bw[k] if isinstance(bw, dict) else bw)

    def logpdf(self, X: np.ndarray) -> np.ndarray:
        """(n, K) log densities (in standardized coordinates; the Jacobian is common to all classes)."""
        Z = ((X - self.mu) / self.sd).T
        return np.column_stack([self.kde[k].logpdf(Z) for k in self.classes])


def fit_fractions(logp: np.ndarray, n_iter: int = 2000, tol: float = 1e-10, f0=None) -> np.ndarray:
    """EM for mixture weights with fixed component densities; logp is (n, K)."""
    n, K = logp.shape
    f = np.full(K, 1 / K) if f0 is None else np.asarray(f0, float)
    m = logp.max(1, keepdims=True)
    P = np.exp(logp - m)
    for _ in range(n_iter):
        R = P * f
        R /= R.sum(1, keepdims=True)
        fn = R.mean(0)
        if np.max(np.abs(fn - f)) < tol:
            f = fn
            break
        f = fn
    return f


def bootstrap_fractions(logp: np.ndarray, n_boot: int = 200, seed: int = 0) -> np.ndarray:
    """Bootstrap over stars: (n_boot, K) fitted fractions."""
    rng = np.random.default_rng(seed)
    n = logp.shape[0]
    f_hat = fit_fractions(logp)
    return np.array([fit_fractions(logp[rng.integers(0, n, n)], f0=f_hat) for _ in range(n_boot)])


def cv_bandwidth(X: np.ndarray, groups: np.ndarray, mu, sd, grid=(0.2, 0.25, 0.3, 0.35, 0.45, 0.6, 0.8),
                 n_fold: int = 4, jitter: float = 0.02, seed: int = 0) -> float:
    """KDE bandwidth (scipy bw_method factor) maximizing the log-likelihood of held-out GROUPS (stars), so that the
    density generalizes to stars not used to build it. X: (n, d) features of one class; groups: star id per row."""
    rng = np.random.default_rng(seed)
    Z = (X - mu) / sd + rng.normal(0, jitter, X.shape)
    ug = np.unique(groups)
    fold = dict(zip(ug, rng.integers(0, n_fold, ug.size)))
    fo = np.array([fold[g] for g in groups])
    best = (-np.inf, None)
    for bw in grid:
        ll = sum(gaussian_kde(Z[fo != f].T, bw_method=bw).logpdf(Z[fo == f].T).sum() for f in range(n_fold))
        if ll > best[0]:
            best = (ll, bw)
    return best[1]
