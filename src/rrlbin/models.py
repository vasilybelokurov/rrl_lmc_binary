"""Explicit timing models for hypothesis testing (docs/PLAN.md, stages 1-3).

All models share the trend (quadratic O-C: P, dP/dt) and the MACHO band offsets (Gaussian priors as pseudo-rows), solved by
GLS, and a Gaussian likelihood with covariance  C = K(theta) + diag(err^2 + s^2)  (K: correlated noise of the model, s: white
jitter). Log-likelihoods are maximum-likelihood values with the linear parameters profiled (ltte._gls_lnl), so that models
can be compared by likelihood ratio on the same footing; their null distributions are calibrated by simulation.

Models implemented here
  H_QP   quasi-periodic intrinsic modulation: K = A^2 exp(-dt^2 / (2 l^2)) cos(2 pi dt / P_q), coherence c = l / P_q
         (c -> infinity: a stationary sinusoid of random phase and Gaussian amplitude)
  H_RN   smooth red noise: K = A^2 exp(-dt^2 / (2 l^2))
  H_LTTE Keplerian light-travel-time orbit (Irwin 1952) + white jitter (the orbit enters as a design column at the best
         nonlinear parameters: amplitude, trend and offsets are linear)
Times and delays in days.
"""
from __future__ import annotations

import numpy as np

from .kepler_fit import _unpack, fit_keplerian
from .ltte import _gls_lnl, ltte_delay
from .oc import _augment, design, orbit_search, period_grid

DAY = 86400.0
COHERENCE = (0.5, 1.0, 2.0, 4.0, np.inf)        # l / P_q grid of H_QP


def qp_kernel(t1, t2, A, Pq, ell):
    """Quasi-periodic covariance A^2 exp(-dt^2/(2 l^2)) cos(2 pi dt / Pq); ell = inf gives A^2 cos(2 pi dt / Pq)."""
    dt = np.asarray(t1)[:, None] - np.asarray(t2)[None, :]
    env = 1.0 if np.isinf(ell) else np.exp(-0.5 * dt ** 2 / ell ** 2)
    return A ** 2 * env * np.cos(2 * np.pi * dt / Pq)


def se_kernel(t1, t2, A, ell):
    dt = np.asarray(t1)[:, None] - np.asarray(t2)[None, :]
    return A ** 2 * np.exp(-0.5 * dt ** 2 / ell ** 2)


def ml_lnl(t, y, err, band, K=None, s=0.0, priors=None, extra=None):
    """Profiled Gaussian log-likelihood of the delays: trend (+ offsets, + `extra` design columns) by GLS, covariance
    K + diag(err^2 + s^2); offset priors appended as pseudo-observations (same convention for every model)."""
    t, y, err = (np.asarray(x, float) for x in (t, y, err))
    X, names = design(t, np.asarray(band), np.average(t, weights=err ** -2))
    if extra is not None:
        X = np.column_stack([X] + list(np.atleast_2d(extra)))
        names = names + [f"x{i}" for i in range(np.atleast_2d(extra).shape[0])]
    C = np.diag(err ** 2 + s ** 2) + (0.0 if K is None else K)
    Xp, yp, sp = _augment(X, y, names, priors)
    if len(yp):
        n, m = C.shape[0], len(yp)
        C = np.block([[C, np.zeros((n, m))], [np.zeros((m, n)), np.diag(sp ** 2)]])
        X, y = np.vstack([X, Xp]), np.r_[y, yp]
    return _gls_lnl(X, y, C)[0]


def _amp_grid(err, n):
    return np.r_[0.0, np.geomspace(0.3, 30, n) * np.median(err)]


def fit_qp(t, y, err, band, priors=None, periods=None, coherence=COHERENCE, n_amp=12, n_s=6):
    """ML fit of H_QP on grids of (P_q, coherence c = l/P_q, A, s). Returns dict: lnl (best over all), and lnl_by_c
    (best for each coherence value, so that a bound c <= c_max can be applied afterwards), with the best parameters."""
    t = np.asarray(t, float)
    if periods is None:
        periods = period_grid(np.ptp(t), p_max_factor=0.5, oversample=2)
    A_grid, s_grid = _amp_grid(err, n_amp), _amp_grid(err, n_s - 1)
    by_c = {}
    for c in coherence:
        best = (-np.inf, None)
        for Pq in periods:
            base = qp_kernel(t, t, 1.0, Pq, c * Pq)
            for A in A_grid[1:]:
                for s in s_grid:
                    try:
                        l = ml_lnl(t, y, err, band, A ** 2 * base, s, priors)
                    except np.linalg.LinAlgError:
                        continue
                    if l > best[0]:
                        best = (l, dict(Pq=float(Pq), A=float(A), s=float(s), c=float(c)))
        by_c[c] = best
    l0 = max(ml_lnl(t, y, err, band, None, s, priors) for s in s_grid)     # A = 0 (no modulation) for reference
    c_best = max(by_c, key=lambda c: by_c[c][0])
    return dict(lnl=max(by_c[c_best][0], l0), lnl_by_c={c: v[0] for c, v in by_c.items()}, best=by_c[c_best][1], lnl_white=l0)


def fit_ltte(t, y, err, band, priors=None, p_max_factor=0.5, s_factors=(0.0, 0.3, 0.6, 1.0, 1.5, 2.0, 3.0, 5.0)):
    """ML fit of H_LTTE (Keplerian + white jitter). P0 and the jitter s1 from the sinusoidal O-C search over
    400 d ... p_max_factor x baseline; one multi-start Keplerian fit at s1; then the likelihood is maximized over the jitter
    grid with the orbit SHAPE (P, e, omega, t_p) fixed and amplitude, trend and offsets re-solved linearly (ml_lnl with the
    orbit as a design column: the same likelihood convention as the other models). Returns dict(lnl, P, A_s, e, omega, t_p, s)."""
    t, y, err = (np.asarray(x, float) for x in (t, y, err))
    band = np.asarray(band).astype(int)
    periods = period_grid(np.ptp(t), p_max_factor=p_max_factor)
    srch = orbit_search(t, y, err, band, periods, priors=priors)
    k = fit_keplerian(t, y, err, band, srch["P_best"], s_jit=srch["s1"], priors=priors, n_boot=0, p_factors=(0.9, 1.0, 1.12))
    P, e, w, ph = _unpack(k["x"])
    g = ltte_delay(t, P, 1.0, e, w, ph * P)
    e0 = np.median(err)
    lnl, s_best = max((ml_lnl(t, y, err, band, None, f * e0, priors, extra=g), f * e0) for f in s_factors)
    return dict(lnl=lnl, P=k["P"], A_s=k["A_s"], e=k["e"], omega=k["omega"], t_p=k["t_p"], s=s_best)
