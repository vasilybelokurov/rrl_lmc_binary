"""Tests of LTTE delays, the mass function and the O-C search (src/rrlbin/ltte.py)."""
import numpy as np

from rrlbin.ltte import a1sini_over_c, kepler_E, ltte_delay, mass_function, oc_search, period_grid


def test_amplitudes_match_journal_table():
    """M1 = 0.65, P = 3000 d, edge-on: 1050, 453, 177 s for M2 = 0.6, 0.2, 0.07 (JOURNAL.md)."""
    got = a1sini_over_c(3000, 0.65, np.array([0.6, 0.2, 0.07]))
    assert np.allclose(got, [1050, 453, 177], atol=1.5)


def test_mass_function_roundtrip():
    M1, M2, sini = 0.65, 0.3, 0.8
    f = mass_function(2000, a1sini_over_c(2000, M1, M2, sini))
    assert np.isclose(f, (M2 * sini) ** 3 / (M1 + M2) ** 2, rtol=1e-10)


def test_kepler_and_circular_limit():
    M = np.linspace(0, 2 * np.pi, 50)
    for e in [0.0, 0.3, 0.9]:
        E = kepler_E(M, e)
        assert np.max(np.abs(E - e * np.sin(E) - M)) < 1e-12
    t = np.linspace(0, 5000, 300)
    assert np.allclose(ltte_delay(t, 1500, 0.01, 0.0, 0.7, 100), 0.01 * np.sin(2 * np.pi * (t - 100) / 1500 + 0.7))


def test_eccentric_peak_to_peak():
    """Peak-to-peak of the Irwin delay is 2 amp sqrt(1 - e^2 cos^2 omega)."""
    t = np.linspace(0, 1000, 200001)
    for e, w in [(0.5, 0.3), (0.8, 2.0)]:
        d = ltte_delay(t, 1000, 1.0, e, w, 0)
        assert np.isclose(np.ptp(d), 2 * np.sqrt(1 - e ** 2 * np.cos(w) ** 2), rtol=1e-4)


def test_oc_search_recovers_signal_and_null_calibration():
    rng = np.random.default_rng(0)
    t = 2200 + 365.25 * np.arange(14) + rng.uniform(0, 100, 14)
    err = np.full(t.size, 250 / 86400)
    periods = period_grid(np.ptp(t))
    # signal: 600 s, P = 2500 d, plus a quadratic
    tau = 600 / 86400 * np.sin(2 * np.pi * t / 2500 + 1) + 1e-9 * (t - 5000) ** 2 + rng.normal(0, 1, t.size) * err
    r = oc_search(t, tau, err, periods)
    assert abs(r["P_best"] - 2500) / 2500 < 0.15 and r["D"] > 15
    # null: D has a sensible distribution (max over periods of ~chi2_2 variables; median of order a few)
    D0 = [oc_search(t, rng.normal(0, 1, t.size) * err, err, periods)["D"] for _ in range(200)]
    assert 1 < np.median(D0) < 12 and np.percentile(D0, 99) < 30


def test_red_search_reduces_to_white_without_random_walk():
    """With q = 0 only, oc_search_red equals oc_search on the same jitter grid."""
    from rrlbin.ltte import oc_search_red
    rng = np.random.default_rng(5)
    t = 2200 + 365.25 * np.arange(12) + rng.uniform(0, 100, 12)
    err = np.full(t.size, 200 / 86400)
    y = 500 / 86400 * np.sin(2 * np.pi * t / 2000) + rng.normal(0, 1, t.size) * err
    P = period_grid(np.ptp(t))
    s = np.r_[0, np.geomspace(0.1, 30, 12)] * 200 / 86400
    a, b = oc_search(t, y, err, P, s_grid=s), oc_search_red(t, y, err, P, s_grid=s, q_grid=np.array([0.0]))
    assert np.allclose(a["Dp"], b["Dp"], atol=1e-6)


def test_oc_search_free_offset_column():
    """A constant offset on a subset of seasons (e.g. MACHO vs OGLE) is absorbed by an indicator column, and the
    signal is still recovered; without the column, the offset biases the search."""
    rng = np.random.default_rng(11)
    t = 450 + 365.25 * np.arange(20) - 5 * 365.25 + rng.uniform(0, 100, 20)   # 1992-2011-like
    macho = t < 1600
    err = np.full(t.size, 150 / 86400)
    y = 600 / 86400 * np.sin(2 * np.pi * t / 2200 + 0.3) + np.where(macho, 2000 / 86400, 0) + rng.normal(0, 1, t.size) * err
    P = period_grid(np.ptp(t))
    r = oc_search(t, y, err, P, X_extra=macho.astype(float)[:, None])
    assert abs(r["P_best"] / 2200 - 1) < 0.1 and abs(r["amp"] * 86400 - 600) < 200


def test_predictive_score_orbit_vs_random_walk():
    """Held-out early seasons: a real orbit gives a positive predictive score on average; a random walk does not."""
    from rrlbin.ltte import predictive_score
    from rrlbin.simulate import delay_random_walk
    rng = np.random.default_rng(21)
    t = -1300 + 365.25 * np.arange(24) + rng.uniform(0, 60, 24)     # 1992-2015-like, one point per year
    err = np.full(t.size, 150 / 86400)
    test = t < 450
    so, sr = [], []
    for _ in range(60):
        y = 800 / 86400 * np.sin(2 * np.pi * t / rng.uniform(1500, 3500) + rng.uniform(0, 6)) + rng.normal(0, 1, t.size) * err
        so.append(predictive_score(t, y, err, test)["score"])
        y = delay_random_walk(t, 800 / 86400, rng) + rng.normal(0, 1, t.size) * err
        sr.append(predictive_score(t, y, err, test)["score"])
    assert np.median(so) > 2 and np.median(so) > np.median(sr) + 2


def test_period_grid_respects_annual_nyquist():
    """Season delays sample ~once per year: the default grid must start above the 2-yr Nyquist period, and a
    long-period signal must not be assigned its annual alias."""
    assert period_grid(8000).min() >= 730
    rng = np.random.default_rng(4)
    t = 2300 + 365.25 * np.arange(16) + rng.uniform(0, 30, 16)
    err = np.full(t.size, 300 / 86400)
    y = 1050 / 86400 * np.sin(2 * np.pi * t / 3000) + rng.normal(0, 1, t.size) * err
    r = oc_search(t, y, err, period_grid(np.ptp(t)))
    assert abs(r["P_best"] / 3000 - 1) < 0.2
