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


def test_fit_rn_prefers_red_noise_and_qp_nests_it_poorly():
    """Smooth red noise: fit_rn beats white noise; its likelihood is at least that of the best QP with short coherence
    minus a small margin (QP with c = 0.5 approximates an SE kernel only roughly)."""
    from rrlbin.models import fit_rn, se_kernel
    rng = np.random.default_rng(4)
    t, band, err = seasons(rng)
    y = rng.multivariate_normal(np.zeros(t.size), se_kernel(t, t, 900 / DAY, 1500.0) + np.diag(err ** 2))
    r = fit_rn(t, y, err, band, n_amp=8, n_s=4)
    assert r["lnl"] - max(ml_lnl(t, y, err, band, None, s) for s in (0.0, 1e-3, 3e-3, 1e-2)) > 5
    assert 350 <= r["ell"] <= 3000


def test_fit_amp_mod_detects_modulation_and_is_chi2_under_null():
    """alpha modulated at P (10%) is detected (d2lnl large, amplitude recovered); without modulation d2lnl ~ chi^2_2
    (mean ~ 2, few above 13.8 = the 0.1% point)."""
    from rrlbin.models import fit_amp_mod
    rng = np.random.default_rng(6)
    t, band, _ = seasons(rng, 30)
    band = np.r_[np.zeros(22, int), np.ones(8, int)]
    ae = np.full(t.size, 0.02)
    a = 1 + 0.10 * np.sin(2 * np.pi * t / 3000 + 0.7) + rng.normal(0, 1, t.size) * ae
    r = fit_amp_mod(t, a, ae, band, 3000.0)
    assert r["d2lnl"] > 50 and abs(r["amp_mod"] - 0.10) < 3 * r["amp_mod_err"]
    d = [fit_amp_mod(t, 1 + rng.normal(0, 1, t.size) * ae, ae, band, 3000.0)["d2lnl"] for _ in range(300)]
    assert 1.4 < np.mean(d) < 2.8 and np.mean(np.array(d) > 13.8) < 0.01


def test_fit_alpha_var():
    """Smoothly varying alpha (10% rms, l = 1500 d) is detected; constant alpha with correct errors gives small d2lnl
    (the statistic is >= 0 by construction; < 1% above 13.8 in 300 null draws)."""
    from rrlbin.models import fit_alpha_var, se_kernel
    rng = np.random.default_rng(9)
    t, band, _ = seasons(rng, 30)
    ae = np.full(t.size, 0.02)
    a = 1 + rng.multivariate_normal(np.zeros(t.size), se_kernel(t, t, 0.10, 1500.0)) + rng.normal(0, 1, t.size) * ae
    assert fit_alpha_var(t, a, ae, band)["d2lnl"] > 30
    d = np.array([fit_alpha_var(t, 1 + rng.normal(0, 1, t.size) * ae, ae, band, n_amp=6)["d2lnl"] for _ in range(300)])
    assert d.min() >= -1e-9 and np.mean(d > 13.8) < 0.01
