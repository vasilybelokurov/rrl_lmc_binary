"""Predictive distributions of future O-C delays under H1 (orbit) and H0 (no orbit), with the measured red-noise model.

Delays: y = X beta [+ orbit(t)] + w(t) + eps,  w ~ GP(0, A^2 exp(-dt^2 / (2 l^2))) (smooth wander), eps ~ N(0, err^2 + s^2).
X: quadratic + MACHO band offsets (with Gaussian priors on the offsets as pseudo-rows). Prediction (universal kriging): the
trend coefficients are estimated by GLS and their uncertainty is propagated:
    mean* = X* b + K*^T C^-1 (y - X b)
    cov*  = K** - K*^T C^-1 K* + R (X^T C^-1 X)^-1 R^T,  R = X* - K*^T C^-1 X
Prediction epochs are OGLE I band (no offsets) and include the white jitter of a future SEASON delay (s^2 + sigma_new^2 added
by the user when comparing with data).
"""
from __future__ import annotations

import numpy as np

from .ltte import ltte_delay
from .oc import _augment, design

DAY = 86400.0


def gp_kernel(t1, t2, A, ell):
    return A ** 2 * np.exp(-0.5 * (np.asarray(t1)[:, None] - np.asarray(t2)[None, :]) ** 2 / ell ** 2)


def krige(t, y, err, band, t_pred, A, ell, s, priors=None):
    """Universal-kriging predictive mean and covariance of the OGLE-band delay at t_pred (trend: quadratic + offsets)."""
    t, y, err, band, t_pred = map(np.asarray, (t, y, err, band, t_pred))
    t_ref = np.average(t, weights=err ** -2)
    X, names = design(t, band, t_ref)
    Xp_rows, yp, sp = _augment(X, y, names, priors)
    C = gp_kernel(t, t, A, ell) + np.diag(err ** 2 + s ** 2)
    # priors on offsets: append as independent pseudo-observations of the coefficients
    if len(yp):
        n, m = C.shape[0], len(yp)
        C = np.block([[C, np.zeros((n, m))], [np.zeros((m, n)), np.diag(sp ** 2)]])
        Xa, ya = np.vstack([X, Xp_rows]), np.r_[y, yp]
    else:
        Xa, ya = X, y
    Ci = np.linalg.inv(C)
    F = np.linalg.inv(Xa.T @ Ci @ Xa)
    b = F @ Xa.T @ Ci @ ya
    Xs, _ = design(t_pred, np.zeros(t_pred.size, int), t_ref)
    Xs_full = np.zeros((t_pred.size, Xa.shape[1]))
    Xs_full[:, :Xs.shape[1]] = Xs          # OGLE band: offset columns are zero
    Ks = np.zeros((Xa.shape[0], t_pred.size))
    Ks[: t.size] = gp_kernel(t, t_pred, A, ell)
    Kss = gp_kernel(t_pred, t_pred, A, ell)
    resid = ya - Xa @ b
    mean = Xs_full @ b + Ks.T @ Ci @ resid
    R = Xs_full - Ks.T @ Ci @ Xa
    cov = Kss - Ks.T @ Ci @ Ks + R @ F @ R.T
    return mean, cov


def predict_h0_h1(t, y, err, band, t_pred, kep, gp_h0, gp_h1, priors=None):
    """Predictive (mean, cov) under H0 (no orbit; GP parameters gp_h0 = (s, A, l) fitted to y) and H1 (orbit `kep`
    = dict(P, A_s, e, omega, t_p) subtracted, GP parameters gp_h1 fitted to the orbit-subtracted delays)."""
    m0, c0 = krige(t, y, err, band, t_pred, gp_h0[1], gp_h0[2], gp_h0[0], priors)
    orb = lambda x: ltte_delay(np.asarray(x), kep["P"], kep["A_s"] / DAY, kep["e"], kep["omega"], kep["t_p"])
    m1, c1 = krige(t, y - orb(t), err, band, t_pred, gp_h1[1], gp_h1[2], gp_h1[0], priors)
    return (m0, c0), (m1 + orb(t_pred), c1)


def interp_prediction(t_grid, mean, cov, t_new):
    """Linear interpolation of a predictive mean (n,) and covariance (n, n) tabulated on t_grid to epochs t_new (m,):
    mean' = W mean, cov' = W cov W^T with W the (m, n) linear-interpolation weights (exact for the frozen grid; the
    error is negligible for a 30-d grid and >= 700-d correlation lengths). t_new must lie within t_grid."""
    t_grid, t_new = np.asarray(t_grid, float), np.atleast_1d(np.asarray(t_new, float))
    if t_new.min() < t_grid[0] or t_new.max() > t_grid[-1]:
        raise ValueError("t_new outside the prediction grid")
    W = np.zeros((t_new.size, t_grid.size))
    i = np.clip(np.searchsorted(t_grid, t_new) - 1, 0, t_grid.size - 2)
    f = (t_new - t_grid[i]) / (t_grid[i + 1] - t_grid[i])
    W[np.arange(t_new.size), i], W[np.arange(t_new.size), i + 1] = 1 - f, f
    return W @ mean, W @ cov @ W.T


def gauss_score(d, mean, cov):
    """ln N(d | mean, cov), chi^2 = r^T C^-1 r and the number of points, for a Gaussian predictive distribution."""
    r = np.asarray(d) - mean
    L = np.linalg.cholesky(cov)
    z = np.linalg.solve(L, r)
    chi2 = float(z @ z)
    lnl = -0.5 * chi2 - np.log(np.diag(L)).sum() - 0.5 * r.size * np.log(2 * np.pi)
    return lnl, chi2, r.size


def split_test(t, y, err, band, P, t_split, priors=None, min_train=8, min_test=3):
    """Out-of-sample orbit test on one O-C series (real or simulated; delays [d], common mode already applied).

    Training = all seasons with t < t_split; test = OGLE I (band 0) seasons with t >= t_split. Everything fitted on the
    training seasons only: robust unwrapping and band alignment, the O-C period search (P0), the H0 red-noise fit (GP, REML),
    the Keplerian orbit (white jitter = the H0 white term) and the H1 red-noise fit after removing the orbit. The test seasons
    are unwrapped sequentially from the last training season (whole cycles only) and scored against the H0/H1 kriging
    predictions (+ err^2 + white^2): ln BF = ln L1 - ln L0 (a conditional predictive log-score difference), chi^2 p-values.
    Returns a dict (empty 'ok': False if too few seasons)."""
    from scipy.stats import chi2 as chi2_dist

    from .kepler_fit import fit_keplerian
    from .oc import align_bands, gp_null, orbit_search, period_grid, robust_unwrap

    t, y, err, band = (np.asarray(x) for x in (t, y, err, band))
    band = band.astype(int)
    o = np.argsort(t)
    t, y, err, band = t[o], y[o], err[o], band[o]
    tr = t < t_split
    te = (~tr) & (band == 0)
    if tr.sum() < min_train or te.sum() < min_test or not np.any(band[tr] == 0):
        return dict(ok=False)
    tt, yt, et, bt = t[tr], y[tr], err[tr], band[tr]
    yt = align_bands(tt, robust_unwrap(tt, yt, bt, P), et, bt, P, priors)
    # test seasons: whole-cycle continuation from the last training OGLE season
    i_last = np.flatnonzero(bt == 0)[-1]
    prev, yn = yt[i_last], y[te].copy()
    for j in range(yn.size):
        yn[j] += P * np.round((prev - yn[j]) / P)
        prev = yn[j]
    periods = period_grid(tt.max() - tt.min())
    P0 = orbit_search(tt, yt, et, bt, periods, priors=priors)["P_best"]
    s0, A0, l0 = gp_null(tt, yt, et, bt, priors)
    kf = fit_keplerian(tt, yt, et, bt, P0, s_jit=s0, priors=priors, n_boot=0)
    kep = dict(P=kf["P"], A_s=kf["A_s"], e=kf["e"], omega=kf["omega"], t_p=kf["t_p"])
    orb = ltte_delay(tt, kep["P"], kep["A_s"] / DAY, kep["e"], kep["omega"], kep["t_p"])
    s1, A1, l1 = gp_null(tt, yt - orb, et, bt, priors)
    tp = t[te]
    (m0, C0), (m1, C1) = predict_h0_h1(tt, yt, et, bt, tp, kep, (s0, A0, l0), (s1, A1, l1), priors)
    e2 = err[te] ** 2
    l1_, c1_, n = gauss_score(yn, m1, C1 + np.diag(e2 + s1 ** 2))
    l0_, c0_, _ = gauss_score(yn, m0, C0 + np.diag(e2 + s0 ** 2))
    return dict(ok=True, n_train=int(tr.sum()), n_test=int(n), P0_train=P0, P_kep_train=kep["P"], A_kep_train_s=kep["A_s"],
                e_kep_train=kep["e"], A_gp_h0_s=A0 * DAY, A_gp_h1_s=A1 * DAY, lnBF=l1_ - l0_, chi2_h1=c1_, p_h1=float(chi2_dist.sf(c1_, n)),
                chi2_h0=c0_, p_h0=float(chi2_dist.sf(c0_, n)), med_err_test_s=float(np.median(err[te]) * DAY))
