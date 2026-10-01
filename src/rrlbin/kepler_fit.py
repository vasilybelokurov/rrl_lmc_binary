"""Keplerian light-travel-time orbit fits to season delays (Part C3).

Model (days): tau(t) = c0 + c1 x + c2 x^2 + sum_b Delta_b [band == b] + A * g(t; P, e, omega, t_p),
g = (1 - e^2) / (1 + e cos nu) * sin(nu + omega) (Irwin 1952; A = a1 sin i / c), Gaussian errors err^2 + s^2 (s fixed, e.g.
the white jitter of the circular fit). Linear parameters (A, quadratic, offsets; offset priors as pseudo-rows) are solved
exactly for each nonlinear set (variable projection); the nonlinear set (log P, sqrt(e) cos w, sqrt(e) sin w, phase) is
optimized by least squares from a grid of starts. Uncertainties: residual bootstrap.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import least_squares

from .ltte import ltte_delay
from .oc import _augment, design

DAY = 86400.0


def _unpack(q):
    P = np.exp(q[0])
    e = min(q[1] ** 2 + q[2] ** 2, 0.95)
    w = np.arctan2(q[2], q[1])
    return P, e, w, np.mod(q[3], 1.0)


def _linear(t, tau, sig, band, q, priors):
    P, e, w, ph = _unpack(q)
    g = ltte_delay(t, P, 1.0, e, w, ph * P)
    X, names = design(t, band, np.average(t, weights=sig ** -2))
    X = np.column_stack([X, g])
    names = names + ["A"]
    Xp, yp, sp = _augment(X, tau, names, priors)
    A = np.vstack([X / sig[:, None], Xp / sp[:, None]]) if len(yp) else X / sig[:, None]
    y = np.r_[tau / sig, yp / sp] if len(yp) else tau / sig
    beta, *_ = np.linalg.lstsq(A, y, rcond=None)
    return beta, names, y - A @ beta, X


def fit_keplerian(t, tau, err, band, P0, s_jit=0.0, priors=None, n_boot=100, seed=0, circular=False, x0=None):
    """Best Keplerian fit near P0 (period searched over 0.6-1.6 P0). Returns a dict with P, A [s], e, omega, t_p, chi2,
    the linear coefficients, the model, and bootstrap percentiles (16, 50, 84) of P, A, e, f(M), K1."""
    t, tau, err, band = map(np.asarray, (t, tau, err, band))
    band = band.astype(int)
    sig = np.sqrt(err ** 2 + s_jit ** 2)

    def resid(q, y):
        if circular:
            q = np.array([q[0], 0.0, 0.0, q[1]])
        return _linear(t, y, sig, band, q, priors)[2]

    def best_fit(y, starts):
        best = None
        for q0 in starts:
            q0 = np.array([q0[0], q0[3]]) if circular else np.asarray(q0)
            try:
                r = least_squares(resid, q0, args=(y,), method="trf", x_scale="jac", max_nfev=400)
            except Exception:
                continue
            if best is None or r.cost < best.cost:
                best = r
        return best

    starts = [[np.log(P0 * f), np.sqrt(e) * np.cos(w), np.sqrt(e) * np.sin(w), ph]
              for f in (0.75, 0.9, 1.0, 1.15, 1.35) for e in (0.0, 0.3, 0.6) for w in ((0.0,) if e == 0 else (0.0, 1.6, 3.1, 4.7))
              for ph in (0.0, 0.25, 0.5, 0.75)] if x0 is None else [list(x0)]
    if not np.isfinite(P0) or P0 <= 0:
        raise ValueError(f"invalid starting period {P0}")
    r = best_fit(tau, starts)
    if r is None:
        raise RuntimeError("no Keplerian start converged")
    q = np.array([r.x[0], 0.0, 0.0, r.x[1]]) if circular else r.x
    beta, names, res, X = _linear(t, tau, sig, band, q, priors)
    P, e, w, ph = _unpack(q)
    out = dict(P=P, A_s=beta[names.index("A")] * DAY, e=e, omega=w, t_p=ph * P, chi2=float(np.sum(res ** 2)), x=np.asarray(r.x),
               n=int(t.size), k=len(beta) + (2 if circular else 4), model=X @ beta, names=names, beta=beta)
    if out["A_s"] < 0:   # sign degeneracy A -> -A with omega -> omega + pi
        out["A_s"], out["omega"] = -out["A_s"], np.mod(w + np.pi, 2 * np.pi)
    # residual bootstrap (refit from the best solution)
    rng = np.random.default_rng(seed)
    data_res = (tau - out["model"]) / sig
    boots = []
    for _ in range(n_boot):
        yb = out["model"] + rng.choice(data_res, data_res.size) * sig
        rb = best_fit(yb, [r.x])
        if rb is None:
            continue
        qb = np.array([rb.x[0], 0.0, 0.0, rb.x[1]]) if circular else rb.x
        bb, nb, _, _ = _linear(t, yb, sig, band, qb, priors)
        Pb, eb, _, _ = _unpack(qb)
        Ab = abs(bb[nb.index("A")]) * DAY
        boots.append((Pb, Ab, eb))
    if boots:
        B = np.array(boots)
        fM = 4 * np.pi ** 2 * (B[:, 1] * 299792458.0) ** 3 / (6.6743e-11 * (B[:, 0] * DAY) ** 2) / 1.98841e30
        K1 = 2 * np.pi * B[:, 1] * 299792.458 / (B[:, 0] * DAY * np.sqrt(1 - B[:, 2] ** 2))
        for name, v in [("P", B[:, 0]), ("A_s", B[:, 1]), ("e", B[:, 2]), ("fM", fM), ("K1", K1)]:
            out[f"{name}_p16"], out[f"{name}_p50"], out[f"{name}_p84"] = np.percentile(v, [16, 50, 84])
    return out


def predict_ogle(fit, t, band, sig, t_pred):
    """OGLE-I-band (band 0: no MACHO offsets) prediction of the full delay model (quadratic + orbit) at t_pred [d],
    from a fit_keplerian result obtained on (t, band, sig)."""
    from .oc import design
    t_ref = np.average(t, weights=sig ** -2)
    Xp, names = design(np.asarray(t_pred, float), np.zeros(len(t_pred), int), t_ref)
    b = np.asarray(fit["beta"])
    nb = fit["names"]
    quad = sum(b[nb.index(k)] * Xp[:, names.index(k)] for k in ("q2", "q1", "q0"))
    orb = ltte_delay(np.asarray(t_pred, float), fit["P"], fit["A_s"] / DAY, fit["e"], fit["omega"], fit["t_p"])
    return quad + orb


def bootstrap_predictions(t, tau, err, band, P0, t_pred, s_jit=0.0, priors=None, n_boot=200, seed=0):
    """Residual-bootstrap envelope of the OGLE-band predicted delay curve at t_pred [d]: refits the Keplerian model to
    model + resampled (error-normalized) residuals. Returns (best-fit prediction, array (n_boot, len(t_pred)))."""
    t, tau, err, band = map(np.asarray, (t, tau, err, band))
    sig = np.sqrt(err ** 2 + s_jit ** 2)
    best = fit_keplerian(t, tau, err, band, P0, s_jit=s_jit, priors=priors, n_boot=0, seed=seed)
    p_best = predict_ogle(best, t, band, sig, t_pred)
    rng = np.random.default_rng(seed)
    res = (tau - best["model"]) / sig
    curves = []
    for _ in range(n_boot):
        yb = best["model"] + rng.choice(res, res.size) * sig
        try:
            fb = fit_keplerian(t, yb, err, band, best["P"], s_jit=s_jit, priors=priors, n_boot=0, x0=best["x"])
        except Exception:
            continue
        curves.append(predict_ogle(fb, t, band, sig, t_pred))
    return p_best, np.array(curves)
