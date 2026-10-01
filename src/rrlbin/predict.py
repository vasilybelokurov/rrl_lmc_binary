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
