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


class _Lik:
    """Precomputed pieces of ml_lnl for one star (design with prior pseudo-rows), and a BATCHED evaluation of the profiled
    Gaussian log-likelihood over grids of correlated-noise amplitude A and white jitter s for a fixed kernel shape K0:
    C(A, s) = A^2 K0 + diag(err^2 + s^2) (data block) and diag(sp^2) (prior block). Identical to ml_lnl element by element
    (tests/test_models.py), ~20x faster (one batched solve per kernel shape instead of one Python call per grid point)."""

    def __init__(self, t, y, err, band, priors=None):
        t, y, err = (np.asarray(x, float) for x in (t, y, err))
        X, names = design(t, np.asarray(band), np.average(t, weights=err ** -2))
        Xp, yp, sp = _augment(X, y, names, priors)
        self.n, self.m = t.size, len(yp)
        self.X = np.vstack([X, Xp]) if self.m else X
        self.y = np.r_[y, yp] if self.m else y
        self.e2 = err ** 2
        self.sp2 = np.asarray(sp, float) ** 2 if self.m else np.empty(0)

    def lnl_grid(self, K0, A_grid, s_grid):
        """lnL for every (A, s) on the grid -> array (len(A_grid), len(s_grid)). K0: (n, n) or None (no correlated noise)."""
        n, m = self.n, self.m
        A2 = np.asarray(A_grid, float)[:, None, None, None] ** 2
        S2 = np.asarray(s_grid, float)[None, :, None, None] ** 2
        N = n + m
        C = np.zeros((len(A_grid), len(s_grid), N, N))
        if K0 is not None:
            C[..., :n, :n] = A2 * K0
        idx = np.arange(n)
        C[..., idx, idx] += self.e2 + S2[..., 0]
        if m:
            j = np.arange(n, N)
            C[..., j, j] = self.sp2
        Z = np.concatenate([self.X, self.y[:, None]], axis=1)              # (N, p+1)
        L = np.linalg.cholesky(C)
        W = np.linalg.solve(L, np.broadcast_to(Z, C.shape[:-2] + Z.shape))  # L^-1 [X, y]
        Xw, yw = W[..., :-1], W[..., -1]
        F = np.einsum("...ki,...kj->...ij", Xw, Xw)
        g = np.einsum("...ki,...k->...i", Xw, yw)
        beta = np.linalg.solve(F, g[..., None])[..., 0]
        q = np.einsum("...k,...k->...", yw, yw) - np.einsum("...i,...i->...", beta, g)
        logdet = np.log(np.diagonal(L, axis1=-2, axis2=-1)).sum(-1)
        return -0.5 * q - logdet - 0.5 * N * np.log(2 * np.pi)


def _amp_grid(err, n):
    return np.r_[0.0, np.geomspace(0.3, 30, n) * np.median(err)]


def s_grid(err, n=8):
    """White-jitter grid shared by ALL models (so that nested likelihoods are comparable): 0 and 0.1-30 x median error."""
    return np.r_[0.0, np.geomspace(0.1, 30, n) * np.median(err)]


def fit_white(t, y, err, band, priors=None, n_s=8, lik=None):
    """H_W: trend + offsets + white jitter only (the common floor of every model). Returns dict(lnl, s)."""
    lik = lik or _Lik(t, y, err, band, priors)
    sg = s_grid(err, n_s)
    l = lik.lnl_grid(None, [0.0], sg)[0]
    k = int(np.argmax(l))
    return dict(lnl=float(l[k]), s=float(sg[k]))


def fit_qp(t, y, err, band, priors=None, periods=None, coherence=COHERENCE, n_amp=12, n_s=8):
    """ML fit of H_QP on grids of (P_q, coherence c = l/P_q, A, s). Returns dict: lnl (best over all), and lnl_by_c
    (best for each coherence value, so that a bound c <= c_max can be applied afterwards), with the best parameters."""
    t = np.asarray(t, float)
    if periods is None:
        periods = period_grid(np.ptp(t), p_max_factor=0.5, oversample=2)
    A_grid, sg = _amp_grid(err, n_amp)[1:], s_grid(err, n_s)
    lik = _Lik(t, y, err, band, priors)
    by_c = {}
    for c in coherence:
        best = (-np.inf, None)
        for Pq in periods:
            try:
                l = lik.lnl_grid(qp_kernel(t, t, 1.0, Pq, c * Pq), A_grid, sg)
            except np.linalg.LinAlgError:
                continue
            i, j = np.unravel_index(int(np.argmax(l)), l.shape)
            if l[i, j] > best[0]:
                best = (float(l[i, j]), dict(Pq=float(Pq), A=float(A_grid[i]), s=float(sg[j]), c=float(c)))
        by_c[c] = best
    l0 = fit_white(t, y, err, band, priors, n_s, lik=lik)["lnl"]     # A = 0 (no modulation): the floor
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


RN_ELLS = (350.0, 700.0, 1500.0, 3000.0, 6000.0)


def fit_rn(t, y, err, band, priors=None, ells=RN_ELLS, n_amp=12, n_s=8):
    """ML fit of H_RN (smooth red noise, squared-exponential kernel) on grids of (l, A, s). Returns dict(lnl, A, ell, s)."""
    t = np.asarray(t, float)
    A_grid, sg = _amp_grid(err, n_amp), s_grid(err, n_s)
    lik = _Lik(t, y, err, band, priors)
    best = (-np.inf, None)
    for ell in ells:
        try:
            l = lik.lnl_grid(se_kernel(t, t, 1.0, ell), A_grid, sg)
        except np.linalg.LinAlgError:
            continue
        i, j = np.unravel_index(int(np.argmax(l)), l.shape)
        if l[i, j] > best[0]:
            best = (float(l[i, j]), dict(A=float(A_grid[i]), ell=float(ell), s=float(sg[j])))
    return dict(lnl=best[0], **best[1])


def fit_amp_mod(t, alpha, alpha_err, band, P, n_s=12):
    """H_BL amplitude channel: is the season amplitude scale alpha_j modulated at the timing period P?

    Model (a): alpha_j = a_b (per-band constant) + white jitter s_a; model (b): (a) + a sinusoid of period P with free amplitude
    and phase (common to all bands; alpha is normalized per band, so the modulation is fractional). An LTTE orbit predicts (a);
    Blazhko-like modulation predicts (b) with the timing and amplitude modulation at the same period. Returns dict with
    lnl0, lnl1, d2lnl = 2 (lnl1 - lnl0) (~ chi^2_2 under (a) if errors are right; calibrate by simulation) and the fitted
    fractional modulation amplitude amp_mod with its error (from the WLS covariance at the best jitter)."""
    t, a, ae = (np.asarray(x, float) for x in (t, alpha, alpha_err))
    band = np.asarray(band).astype(int)
    bands = np.unique(band)
    X0 = np.column_stack([(band == b).astype(float) for b in bands])
    w = 2 * np.pi * t / P
    X1 = np.column_stack([X0, np.sin(w), np.cos(w)])
    sg = np.r_[0.0, np.geomspace(0.1, 30, n_s) * np.median(ae)]

    def best(X):
        out = (-np.inf, None, None)
        for s in sg:
            v = ae ** 2 + s ** 2
            sw = 1 / np.sqrt(v)
            beta, *_ = np.linalg.lstsq(X * sw[:, None], a * sw, rcond=None)
            r = a - X @ beta
            l = -0.5 * np.sum(r ** 2 / v + np.log(2 * np.pi * v))
            if l > out[0]:
                cov = np.linalg.inv((X * (1 / v)[:, None]).T @ X)
                out = (l, beta, cov)
        return out
    l0, _, _ = best(X0)
    l1, b1, c1 = best(X1)
    amp = float(np.hypot(b1[-2], b1[-1]))
    g = np.array([b1[-2], b1[-1]]) / max(amp, 1e-12)
    amp_err = float(np.sqrt(g @ c1[-2:, -2:] @ g))
    return dict(lnl0=l0, lnl1=l1, d2lnl=2 * (l1 - l0), amp_mod=amp, amp_mod_err=amp_err)


def fit_alpha_var(t, alpha, alpha_err, band, ells=RN_ELLS, n_amp=10, s_floor=0.0):
    """Amplitude channel of H_BL in its general form: does the season amplitude scale alpha_j vary beyond its errors?

    H0 (what an LTTE orbit predicts): alpha_j = per-band constant, variance alpha_err^2 + s_floor^2 (s_floor: population
    calibration floor, fixed). H1: + smooth variability, SE GP on any of the time scales `ells` (+ the same floor), amplitude on
    a grid. Returns dict(lnl0, lnl1, d2lnl = 2 (lnl1 - lnl0), A [fractional], ell). Variability at the timing period is a special
    case; stage 1b showed coherent intrinsic modulators vary in amplitude but mostly not as a synchronous sinusoid."""
    t, a, ae = (np.asarray(x, float) for x in (t, alpha, alpha_err))
    band = np.asarray(band).astype(int)
    X = np.column_stack([(band == b).astype(float) for b in np.unique(band)])
    D = np.diag(ae ** 2 + s_floor ** 2)
    l0 = _gls_lnl(X, a, D)[0]
    best = (l0, 0.0, np.nan)
    for ell in ells:
        base = se_kernel(t, t, 1.0, ell)
        for A in np.geomspace(0.003, 0.3, n_amp):
            try:
                l = _gls_lnl(X, a, D + A ** 2 * base)[0]
            except np.linalg.LinAlgError:
                continue
            if l > best[0]:
                best = (l, float(A), float(ell))
    return dict(lnl0=l0, lnl1=best[0], d2lnl=2 * (best[0] - l0), A=best[1], ell=best[2])
