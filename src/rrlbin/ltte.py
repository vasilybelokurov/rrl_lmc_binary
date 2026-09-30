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


def oc_search(t, tau, err, periods, s_grid=None, t_ref=None):
    """Circular-orbit LTTE search on season delays with a quadratic ephemeris term and free jitter.

    H0: tau = c0 + c1 x + c2 x^2;  H1: H0 + A sin(2 pi t/P) + B cos(2 pi t/P); both with variance err^2 + s^2,
    s profiled on s_grid. x = (t - t_ref)/1000 d.

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
    return dict(D=float(Dp[k]), P_best=float(periods[k]), amp=float(np.hypot(b[3], b[4])),
                jit0=float(np.sqrt(s20)), jit1=float(np.sqrt(s21)), Dp=Dp)


def period_grid(baseline, p_min=300.0, p_max_factor=2.0, oversample=5):
    """Periods uniform in frequency from 1/(p_max_factor*baseline) to 1/p_min."""
    f_min, f_max = 1 / (p_max_factor * baseline), 1 / p_min
    n = int(np.ceil((f_max - f_min) * baseline * oversample)) + 1
    return 1 / np.linspace(f_max, f_min, n)
