"""Light-travel-time effect (LTTE): orbital delays and an O-C search with jitter.

Delay of the pulsating primary (Irwin 1952, ApJ 116, 211; without the constant e sin(omega) term):

    tau(t) = (a1 sin i / c) * (1 - e^2) / (1 + e cos nu) * sin(nu + omega)

Units: time in days, a1 sin i / c in days unless stated otherwise.
"""
from __future__ import annotations

import numpy as np

G = 6.67430e-11
MSUN = 1.98841e30
C = 299792458.0
DAY = 86400.0


def a1sini_over_c(P_orb, M1, M2, sini=1.0):
    """Projected light-travel semi-amplitude a1 sin(i)/c [s] of the primary; P_orb [d], masses [Msun]."""
    P = np.asarray(P_orb) * DAY
    a = (G * (np.asarray(M1) + M2) * MSUN * P ** 2 / (4 * np.pi ** 2)) ** (1 / 3)
    return a * M2 / (M1 + M2) * sini / C


def mass_function(P_orb, a1sini_s):
    """f(M2) [Msun] from P_orb [d] and a1 sin i / c [s]."""
    return 4 * np.pi ** 2 * (np.asarray(a1sini_s) * C) ** 3 / (G * (np.asarray(P_orb) * DAY) ** 2) / MSUN


def kepler_E(M, e, n_iter=30):
    """Solve Kepler's equation E - e sin E = M (Newton; e < 0.99)."""
    E = M + e * np.sin(M)
    for _ in range(n_iter):
        E = E - (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
    return E


def ltte_delay(t, P_orb, amp, e=0.0, omega=0.0, t_peri=0.0):
    """LTTE delay [same units as amp] at times t [d]. amp = a1 sin i / c."""
    M = 2 * np.pi * (np.asarray(t) - t_peri) / P_orb
    E = kepler_E(np.mod(M, 2 * np.pi), e)
    nu = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2), np.sqrt(1 - e) * np.cos(E / 2))
    return amp * (1 - e ** 2) / (1 + e * np.cos(nu)) * np.sin(nu + omega)


# ------------------------------------------------------------------ O-C search
def _profile_lnl(X, y, var, s2_grid):
    """max over jitter s^2 of the Gaussian lnL for a WLS fit y ~ X beta with variance var + s^2.
    X: (n, p); returns (lnL_max, s2_best, beta_best)."""
    best = (-np.inf, 0.0, None)
    for s2 in s2_grid:
        v = var + s2
        sw = 1 / np.sqrt(v)
        beta, *_ = np.linalg.lstsq(X * sw[:, None], y * sw, rcond=None)
        r = y - X @ beta
        lnl = -0.5 * np.sum(r ** 2 / v + np.log(2 * np.pi * v))
        if lnl > best[0]:
            best = (lnl, s2, beta)
    return best


def oc_search(t, tau, err, periods, s_grid=None, t_ref=None, X_extra=None):
    """Circular-orbit LTTE search on season delays with a quadratic ephemeris term and free jitter.

    H0: tau = c0 + c1 x + c2 x^2;  H1: H0 + A sin(2 pi t/P) + B cos(2 pi t/P); both with variance err^2 + s^2,
    s profiled on s_grid. x = (t - t_ref)/1000 d. X_extra: optional extra nuisance columns (n, q) in both
    hypotheses, e.g. an indicator of MACHO seasons for a free MACHO-OGLE delay offset.

    Returns dict with D = 2 max_P [lnL1(P) - lnL0], best period, semi-amplitude, jitter under H0/H1, and D(P).
    """
    t, tau, err = map(np.asarray, (t, tau, err))
    if t_ref is None:
        t_ref = np.average(t, weights=err ** -2)
    if s_grid is None:
        s_grid = np.r_[0.0, np.geomspace(0.1, 30, 30) * np.median(err)]
    s2 = s_grid ** 2
    x = (t - t_ref) / 1000.0
    X0 = np.vander(x, 3)
    if X_extra is not None:
        X0 = np.column_stack([X0, X_extra])
    nq = X0.shape[1]
    var = err ** 2
    l0, s20, _ = _profile_lnl(X0, tau, var, s2)
    Dp = np.empty(len(periods))
    fits = []
    for k, P in enumerate(periods):
        w = 2 * np.pi * t / P
        X1 = np.column_stack([X0, np.sin(w), np.cos(w)])
        l1, s21, b = _profile_lnl(X1, tau, var, s2)
        Dp[k] = 2 * (l1 - l0)
        fits.append((s21, b))
    k = int(np.argmax(Dp))
    s21, b = fits[k]
    return dict(D=float(Dp[k]), P_best=float(periods[k]), amp=float(np.hypot(b[nq], b[nq + 1])),
                jit0=float(np.sqrt(s20)), jit1=float(np.sqrt(s21)), Dp=Dp)


def period_grid(baseline, p_min=300.0, p_max_factor=2.0, oversample=5):
    """Periods uniform in frequency from 1/(p_max_factor*baseline) to 1/p_min."""
    f_min, f_max = 1 / (p_max_factor * baseline), 1 / p_min
    n = int(np.ceil((f_max - f_min) * baseline * oversample)) + 1
    return 1 / np.linspace(f_max, f_min, n)


# ------------------------------------------------------------------ red-noise (random-walk phase) null
def _gls_lnl(X, y, C):
    """Gaussian lnL of a GLS fit y ~ X beta with covariance C (beta at its GLS optimum)."""
    L = np.linalg.cholesky(C)
    Xw = np.linalg.solve(L, X)
    yw = np.linalg.solve(L, y)
    beta, *_ = np.linalg.lstsq(Xw, yw, rcond=None)
    r = yw - Xw @ beta
    return -0.5 * (r @ r) - np.sum(np.log(np.diag(L))) - 0.5 * y.size * np.log(2 * np.pi), beta


def oc_search_red(t, tau, err, periods, s_grid=None, q_grid=None, t_ref=None):
    """As oc_search, but H0 and H1 both include a random-walk (Wiener) phase term:

        C_jk = (err_j^2 + s^2) delta_jk + q * min(t_j - t0, t_k - t0),   t0 = min(t),

    with s and q profiled on grids (q in d^2/d: variance of the delay growing linearly with time).
    The quadratic ephemeris absorbs the random walk's mean and slope. D_red = 2 max_P [lnL1 - lnL0].
    """
    t, tau, err = map(np.asarray, (t, tau, err))
    if t_ref is None:
        t_ref = np.average(t, weights=err ** -2)
    e0 = np.median(err)
    T = np.ptp(t)
    if s_grid is None:
        s_grid = np.r_[0.0, np.geomspace(0.1, 30, 12) * e0]
    if q_grid is None:   # rms random-walk excursion over the baseline from 0.1 to 30 x the median error
        q_grid = np.r_[0.0, (np.geomspace(0.1, 30, 12) * e0) ** 2 / T]
    x = (t - t_ref) / 1000.0
    X0 = np.vander(x, 3)
    W = np.minimum.outer(t - t.min(), t - t.min())
    Cs = [np.diag(err ** 2 + s ** 2) + q * W for s in s_grid for q in q_grid]
    l0 = max(_gls_lnl(X0, tau, C)[0] for C in Cs)
    Dp = np.empty(len(periods))
    best = (-np.inf, None)
    for k, P in enumerate(periods):
        w = 2 * np.pi * t / P
        X1 = np.column_stack([X0, np.sin(w), np.cos(w)])
        l1, b1 = -np.inf, None
        for C in Cs:
            l, b = _gls_lnl(X1, tau, C)
            if l > l1:
                l1, b1 = l, b
        Dp[k] = 2 * (l1 - l0)
        if l1 > best[0]:
            best = (l1, b1)
    k = int(np.argmax(Dp))
    b = best[1]
    return dict(D=float(Dp[k]), P_best=float(periods[k]), amp=float(np.hypot(b[3], b[4])), Dp=Dp)


# ------------------------------------------------------------------ predictive (out-of-sample) test
def _design(t, t_ref, P_orb, X_extra):
    x = (t - t_ref) / 1000.0
    X = np.vander(x, 3)
    if X_extra is not None:
        X = np.column_stack([X, X_extra])
    if P_orb is not None:
        w = 2 * np.pi * t / P_orb
        X = np.column_stack([X, np.sin(w), np.cos(w)])
    return X


def predictive_lnl(t_tr, y_tr, e_tr, t_te, y_te, e_te, P_orb=None, X_tr=None, X_te=None, s_grid=None):
    """Gaussian predictive log-likelihood of test delays given a linear model fitted to training delays.

    Model: quadratic (+ extra nuisance columns) (+ circular orbit at fixed P_orb), white jitter s (profiled on the
    training set). The predictive covariance includes the parameter uncertainty: C = diag(e_te^2 + s^2) + X Cov_beta X^T.
    """
    t_ref = np.average(t_tr, weights=e_tr ** -2)
    A = _design(t_tr, t_ref, P_orb, X_tr)
    B = _design(t_te, t_ref, P_orb, X_te)
    if s_grid is None:
        s_grid = np.r_[0.0, np.geomspace(0.1, 30, 30) * np.median(e_tr)]
    _, s2, beta = _profile_lnl(A, y_tr, e_tr ** 2, s_grid ** 2)
    v = e_tr ** 2 + s2
    cov_b = np.linalg.pinv(A.T @ (A / v[:, None]))
    C = np.diag(e_te ** 2 + s2) + B @ cov_b @ B.T
    r = y_te - B @ beta
    L = np.linalg.cholesky(C)
    z = np.linalg.solve(L, r)
    return float(-0.5 * z @ z - np.sum(np.log(np.diag(L))) - 0.5 * r.size * np.log(2 * np.pi))


def predictive_score(t, y, e, test, X_extra=None, periods=None):
    """Out-of-sample test of an orbit: search the training seasons (~test) for the best circular orbit, then
    score = lnL_pred(test | orbit) - lnL_pred(test | quadratic only). Positive = the orbit predicts the held-out seasons.
    Returns dict(score, P_train, D_train)."""
    tr = ~test
    Xtr = None if X_extra is None else X_extra[tr]
    Xte = None if X_extra is None else X_extra[test]
    if periods is None:
        periods = period_grid(np.ptp(t[tr]))
    r = oc_search(t[tr], y[tr], e[tr], periods, X_extra=Xtr)
    l1 = predictive_lnl(t[tr], y[tr], e[tr], t[test], y[test], e[test], r["P_best"], Xtr, Xte)
    l0 = predictive_lnl(t[tr], y[tr], e[tr], t[test], y[test], e[test], None, Xtr, Xte)
    return dict(score=l1 - l0, P_train=r["P_best"], D_train=r["D"])


def fit_red_null(t, tau, err, X_extra=None, s_grid=None, q_grid=None):
    """ML white jitter s [d] and random-walk strength q [d^2/d] under H0 (quadratic + nuisance columns), on grids.
    Returns dict(s, q, lnL). The rms random-walk excursion over the baseline T is sqrt(q T)."""
    t, tau, err = map(np.asarray, (t, tau, err))
    e0, T = np.median(err), np.ptp(t)
    if s_grid is None:
        s_grid = np.r_[0.0, np.geomspace(0.1, 30, 16) * e0]
    if q_grid is None:
        q_grid = np.r_[0.0, (np.geomspace(0.1, 30, 16) * e0) ** 2 / T]
    X = np.vander((t - np.average(t, weights=err ** -2)) / 1000.0, 3)
    if X_extra is not None:
        X = np.column_stack([X, X_extra])
    W = np.minimum.outer(t - t.min(), t - t.min())
    best = (-np.inf, 0.0, 0.0)
    for s in s_grid:
        for q in q_grid:
            l, _ = _gls_lnl(X, tau, np.diag(err ** 2 + s ** 2) + q * W)
            if l > best[0]:
                best = (l, s, q)
    return dict(lnL=float(best[0]), s=float(best[1]), q=float(best[2]))
