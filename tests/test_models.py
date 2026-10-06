"""Tests of the explicit timing models (src/rrlbin/models.py)."""
import numpy as np
from scipy.stats import multivariate_normal

from rrlbin.ltte import ltte_delay
from rrlbin.models import fit_ltte, fit_qp, ml_lnl, qp_kernel

DAY = 86400.0


def seasons(rng, n=26):
    t = np.r_[np.arange(-1000.0, 1700, 365.25), np.arange(2200.0, 11100, 365.25)][:n] + rng.normal(0, 30, n)
    return np.sort(t), np.zeros(n, int), np.full(n, 150.0 / DAY)


def test_qp_kernel_psd_and_limits():
    t = np.linspace(0, 8000, 30)
    for ell in (500.0, 3000.0, np.inf):
        K = qp_kernel(t, t, 1.0, 2500.0, ell)
        assert np.all(np.linalg.eigvalsh(K) > -1e-9)
    assert np.allclose(np.diag(qp_kernel(t, t, 2.0, 2500.0, 1000.0)), 4.0)


def test_ml_lnl_matches_brute_force():
    """Profiled likelihood: equals the Gaussian density at the GLS trend for a pure-trend model (no priors)."""
    rng = np.random.default_rng(1)
    t, band, err = seasons(rng)
    K = qp_kernel(t, t, 300 / DAY, 3000.0, 6000.0)
    y = 1e-3 + 2e-4 * (t - 5000) / 1000 + rng.multivariate_normal(np.zeros(t.size), K + np.diag(err ** 2))
    l = ml_lnl(t, y, err, band, K, 0.0)
    from rrlbin.oc import design
    X, _ = design(t, band, np.average(t, weights=err ** -2))
    C = K + np.diag(err ** 2)
    Ci = np.linalg.inv(C)
    beta = np.linalg.solve(X.T @ Ci @ X, X.T @ Ci @ y)
    assert np.isclose(l, multivariate_normal(X @ beta, C).logpdf(y))


def test_fit_ltte_recovers_orbit_and_beats_white():
    rng = np.random.default_rng(2)
    t, band, err = seasons(rng)
    y = ltte_delay(t, 3000.0, 1200 / DAY, 0.4, 1.0, 500.0) + rng.normal(0, 1, t.size) * err
    r = fit_ltte(t, y, err, band)
    assert abs(r["P"] / 3000 - 1) < 0.05 and abs(r["A_s"] / 1200 - 1) < 0.15 and abs(r["e"] - 0.4) < 0.15
    assert r["lnl"] - max(ml_lnl(t, y, err, band, None, s) for s in (0.0, 1e-3, 3e-3)) > 20


def test_fit_qp_recovers_period_of_coherent_modulation():
    rng = np.random.default_rng(3)
    t, band, err = seasons(rng)
    K = qp_kernel(t, t, 1000 / DAY, 3000.0, 4 * 3000.0)
    y = rng.multivariate_normal(np.zeros(t.size), K + np.diag(err ** 2))
    r = fit_qp(t, y, err, band, coherence=(1.0, 4.0, np.inf), n_amp=8, n_s=4)
    assert abs(r["best"]["Pq"] / 3000 - 1) < 0.15
    assert r["lnl"] - r["lnl_white"] > 10
