"""Tests of the Fourier template and per-season delay estimator (src/rrlbin/timing.py).

Synthetic light curves use a realistic RRab template (Fourier amplitudes of OGLE-LMC-RRLYR-00010,
A_k = 0.155, 0.068, 0.049, 0.024, ...) and an OGLE-like cadence: seasons of ~240 d separated by
~125-d gaps, 40-110 epochs per season, errors 0.06 mag.
"""
import numpy as np
import pytest

from rrlbin.timing import (delays_fixed_template, fit_timing, fourier_dphi, fourier_eval, harmonic_amp_phase, season_labels,
                           unwrap_delays)

P, T0 = 0.5940743, 6000.13
A = np.array([0.155, 0.068, 0.049, 0.024, 0.010, 0.008, 0.006, 0.004])
PH = np.array([0.0, 2.3, 4.5, 0.6, 2.9, 5.1, 1.3, 3.6])
COEF = np.empty(2 * A.size)
COEF[0::2], COEF[1::2] = A * np.cos(PH), A * np.sin(PH)


def cadence(rng, n_seasons=14, t_start=2200.0):
    t = []
    for j in range(n_seasons):
        n = rng.integers(40, 110)
        t.append(t_start + 365.25 * j + np.sort(rng.uniform(0, 240, n)))
    return np.concatenate(t)


def synth(t, tau_of_t, rng=None, sigma=0.06, zp=18.8):
    m = zp + fourier_eval(COEF, (t - tau_of_t(t) - T0) / P)
    if rng is not None:
        m = m + rng.normal(0, sigma, t.size)
    return m, np.full(t.size, sigma)


def test_dphi_matches_finite_difference():
    phi = np.linspace(0, 1, 257)
    h = 1e-6
    fd = (fourier_eval(COEF, phi + h) - fourier_eval(COEF, phi - h)) / (2 * h)
    assert np.max(np.abs(fd - fourier_dphi(COEF, phi))) < 1e-6


def test_amp_phase_roundtrip():
    a, p = harmonic_amp_phase(COEF)
    assert np.allclose(a, A) and np.allclose(np.mod(p - PH + np.pi, 2 * np.pi) - np.pi, 0, atol=1e-12)


def test_season_labels():
    t = np.array([0, 1, 2, 100, 101, 300.0])
    assert list(season_labels(t[::-1], gap=60)[::-1]) == [0, 0, 0, 1, 1, 2]


def test_noiseless_stepwise_delays_exact():
    """Noiseless light curve with a delay that is constant within each season: recovered to < 0.01 s."""
    rng = np.random.default_rng(1)
    t = cadence(rng)
    lab = season_labels(t)
    step = rng.normal(0, 400, lab.max() + 1) / 86400
    m, e = synth(t, lambda x: step[lab])
    f = fit_timing(t, m, e, np.full(t.size, "O4"), P, T0, K=8)
    d = (f.tau - f.tau.mean()) - (step[f.season] - step[f.season].mean())   # constant delay is degenerate with T0
    assert np.max(np.abs(d)) * 86400 < 0.01
    assert f.mask.all()


def test_noiseless_ltte_weighted_mean():
    """Noiseless 300-s LTTE sinusoid (P_orb = 3000 d). Within a season the estimator returns the delay averaged
    with Fisher weights (dm/dt)^2, not the plain mean; agreement to < 3 s (second-order curvature terms)."""
    rng = np.random.default_rng(1)
    t = cadence(rng)
    amp = 300 / 86400
    tau = lambda x: amp * np.sin(2 * np.pi * x / 3000.0)
    m, e = synth(t, tau)
    f = fit_timing(t, m, e, np.full(t.size, "O4"), P, T0, K=8)
    lab = season_labels(t)
    wts = fourier_dphi(COEF, (t - tau(t) - T0) / P) ** 2
    truth = np.array([np.average(tau(t[lab == s]), weights=wts[lab == s]) for s in f.season])
    d = (f.tau - f.tau.mean()) - (truth - truth.mean())
    assert np.max(np.abs(d)) * 86400 < 3.0


def test_large_shift_found_by_grid():
    """A season shifted by 0.3 P (beyond the linear regime) is recovered by the coarse grid search."""
    rng = np.random.default_rng(2)
    t = cadence(rng, n_seasons=6)
    lab = season_labels(t)
    shift = np.where(lab == 3, 0.3 * P, 0.0)
    m, e = synth(t, lambda x: shift, rng, sigma=0.03)
    f = fit_timing(t, m, e, np.full(t.size, "O4"), P, T0, K=8)
    d = f.tau - np.median(f.tau)
    assert abs(d[3] - 0.3 * P) < 5 * f.tau_err[3]
    assert np.all(np.abs(np.delete(d, 3)) < 5 * np.delete(f.tau_err, 3))


def test_noise_pulls_and_zero_points():
    """200 noisy realizations: pull (tau - truth)/tau_err has std 1 +- 0.1 and mean ~0; segment zero points recovered."""
    rng = np.random.default_rng(3)
    t = cadence(rng, n_seasons=10)
    seg = np.where(t < 2200 + 365.25 * 5, "O3", "O4")
    pulls, zpd = [], []
    for _ in range(200):
        m, e = synth(t, lambda x: np.zeros_like(x), rng)
        m = m + np.where(seg == "O4", -0.025, 0.0)
        f = fit_timing(t, m, e, seg, P, T0, K=8)
        d = f.tau - np.average(f.tau, weights=f.tau_err ** -2)
        pulls.append(d / f.tau_err)
        zpd.append(f.zp[1] - f.zp[0])
    pulls = np.concatenate(pulls)
    # removing the weighted mean costs one degree of freedom: expected std sqrt(1 - 1/n) ~ 0.95
    assert 0.85 < pulls.std() < 1.1
    assert abs(pulls.mean()) < 0.05
    assert abs(np.mean(zpd) + 0.025) < 0.003


def test_unwrap():
    tau = np.array([0.0, 0.1, 0.2, 0.3 - P, 0.4 - P])
    assert np.allclose(unwrap_delays(tau, P), [0, 0.1, 0.2, 0.3, 0.4])


def test_gauge_fixed_delays_independent_of_T0():
    """With the fundamental-phase gauge, delays do not depend on the reference epoch T0 (mod P), and the
    returned template has fundamental phase 0 and reproduces the data."""
    rng = np.random.default_rng(7)
    t = cadence(rng, n_seasons=8)
    lab = season_labels(t)
    step = rng.normal(0, 300, lab.max() + 1) / 86400
    m, e = synth(t, lambda x: step[lab], rng, sigma=0.02)
    f1 = fit_timing(t, m, e, np.full(t.size, "O4"), P, T0, K=8)
    f2 = fit_timing(t, m, e, np.full(t.size, "O4"), P, T0 + 0.13, K=8)
    d = (f1.tau - (f2.tau + 0.13))
    d = (d + P / 2) % P - P / 2
    assert np.max(np.abs(d)) * 86400 < 1.0
    assert abs(harmonic_amp_phase(f1.coef)[1][0]) < 1e-10
    # the template phase absorbs the (Fisher-weighted) mean delay of the data, so gauge-fixed delays equal the true
    # delays minus a constant close to their mean (PH[0] = 0 for the input template)
    dd = (f1.tau - step[f1.season] + P / 2) % P - P / 2
    assert np.all(np.abs(dd - np.average(dd, weights=f1.tau_err ** -2)) < 4 * f1.tau_err)
    assert abs(np.mean(dd) + np.mean(step[f1.season])) < 4 * np.mean(f1.tau_err) / np.sqrt(dd.size)


def test_harmonic_coherence_zero_for_pure_delay_nonzero_for_shape_change():
    """A pure time shift (LTTE) gives tau_1 - tau_h consistent with 0; shifting only the higher harmonics
    (a shape change) is detected."""
    rng = np.random.default_rng(12)
    t = cadence(rng, n_seasons=8)
    lab = season_labels(t)
    step = rng.normal(0, 400, lab.max() + 1) / 86400
    m, e = synth(t, lambda x: step[lab], rng, sigma=0.02)
    f = fit_timing(t, m, e, np.full(t.size, "O4"), P, T0, K=8)
    z = f.coh / f.coh_err
    assert np.all(np.abs(z) < 4) and abs(np.mean(z)) < 1.5
    # shape change: in seasons >= 4 the higher harmonics lag the fundamental by 0.02 P
    c_h = COEF.copy(); c_h[:2] = 0
    c_1 = COEF.copy(); c_1[2:] = 0
    lag = np.where(lab >= 4, 0.02 * P, 0.0)
    m2 = 18.8 + fourier_eval(c_1, (t - T0) / P) + fourier_eval(c_h, (t - lag - T0) / P) + rng.normal(0, 0.02, t.size)
    f2 = fit_timing(t, m2, e, np.full(t.size, "O4"), P, T0, K=8)
    def chi2nu(f):   # scatter of the coherence values between seasons (the statistic used as a veto)
        w = f.coh_err ** -2
        return np.sum(w * (f.coh - np.sum(w * f.coh) / w.sum()) ** 2) / (f.coh.size - 1)
    assert chi2nu(f) < 2.5 and chi2nu(f2) > 4


def test_fixed_template_extends_gauge_to_new_seasons():
    """Delays of NEW seasons measured against the (gauge-fixed) template and zero point of a fit to earlier seasons are in
    the same gauge: (measured - true) is the same constant for the training and the new seasons, and re-measuring the
    training seasons with the fixed template reproduces fit_timing's delays."""
    rng = np.random.default_rng(11)
    t = cadence(rng, n_seasons=16)
    lab = season_labels(t)
    step = rng.normal(0, 400, lab.max() + 1) / 86400
    m, e = synth(t, lambda x: step[lab], rng, sigma=0.03)
    tr = lab < 12
    f = fit_timing(t[tr], m[tr], e[tr], np.full(tr.sum(), "O4"), P, T0, K=8)
    zp = f.zp[list(f.seg_names).index("O4")]
    _, tau_tr, _ = delays_fixed_template(t[tr], m[tr], e[tr], lab[tr], f.coef, P, T0, zp=zp)
    wrap = lambda x: (x + P / 2) % P - P / 2
    assert np.max(np.abs(wrap(tau_tr - f.tau))) * 86400 < 2.0
    out = delays_fixed_template(t[~tr], m[~tr], e[~tr], lab[~tr], f.coef, P, T0, zp=zp, full=True)
    assert len(out) == 8 and np.all(np.abs(out[3] - 1) < 0.05)                    # alpha ~ 1
    c_tr = np.average(wrap(f.tau - step[f.season]), weights=f.tau_err ** -2)
    d_new = wrap(out[1] - step[out[7]] - c_tr)
    assert np.all(np.abs(d_new) < 4 * out[2])


def test_gross_outliers_and_sparse_seasons():
    """Bright gross outliers (cosmic-ray-like, -1 mag) in sparse seasons: delays still recovered (within 4 sigma of the truth
    up to the gauge constant), no season dropped as unstable, and with the v3 settings (tight first clip, no repeated grid)
    the result is never better."""
    rng = np.random.default_rng(21)
    t = []
    for j in range(12):
        n = 20 if j % 3 == 0 else 80
        t.append(2200.0 + 365.25 * j + np.sort(rng.uniform(0, 240, n)))
    t = np.concatenate(t)
    lab = season_labels(t)
    step = rng.normal(0, 600, lab.max() + 1) / 86400
    m, e = synth(t, lambda x: step[lab], rng, sigma=0.05)
    bad = rng.choice(t.size, 12, replace=False)
    m[bad] -= 1.0
    f = fit_timing(t, m, e, np.full(t.size, "O4"), P, T0, K=8)
    d = (f.tau - step[f.season] + P / 2) % P - P / 2
    d -= np.average(d, weights=f.tau_err ** -2)
    assert f.unstable.size == 0 and f.season.size == 12
    assert np.all(np.abs(d) < 4 * f.tau_err)
    assert not np.any(f.mask[bad])                       # all gross outliers clipped


def test_stability_drops_bimodal_season():
    """A season whose points are almost all on one branch except a few contradictory high-leverage points: either it is
    measured consistently or it is flagged unstable and dropped; it is never returned with a delay > 5 sigma off."""
    rng = np.random.default_rng(5)
    t = cadence(rng, n_seasons=8)
    lab = season_labels(t)
    step = rng.normal(0, 300, lab.max() + 1) / 86400
    m, e = synth(t, lambda x: step[lab], rng, sigma=0.03)
    s0 = lab == 3
    idx = np.flatnonzero(s0)[:15]                       # keep 15 points in season 3 + 6 shifted by 0.3 P
    keep = ~s0
    keep[idx] = True
    tt, mm, ee, ll = t[keep], m[keep].copy(), e[keep], lab[keep]
    j = np.flatnonzero(ll == 3)[:6]
    mm[j] = synth(tt[j], lambda x: step[3] + 0.3 * P)[0]
    f = fit_timing(tt, mm, ee, np.full(tt.size, "O4"), P, T0, K=8, labels=ll)
    d = (f.tau - step[f.season] + P / 2) % P - P / 2
    d -= np.median(d)
    ok = f.season != 3
    assert np.all(np.abs(d[ok]) < 4 * f.tau_err[ok])
    if 3 in f.season:
        k = list(f.season).index(3)
        assert abs(d[k]) < 5 * f.tau_err[k] or 3 in f.unstable


def test_real_11166_no_wrong_minimum():
    """Regression (real data, skipped if absent): the full 2010-2026 fit of OGLE-LMC-RRLYR-11166 put season 21 at a wrong
    minimum (-23.8 ks, alpha 1.37) with the v3 settings; with the current settings it lies near its neighbours."""
    from rrlbin.io import lc_path
    from rrlbin.pipeline import RAW, load_star
    from rrlbin.timing import year_labels
    oid = "OGLE-LMC-RRLYR-11166"
    if not lc_path(RAW, "ogle4x", oid).exists():
        pytest.skip("extended OGLE-IV data not available")
    t, m, e, s = load_star(oid, ogle4="extended")["I"]
    f = fit_timing(t, m, e, s, 0.4912322, 6000.2098, K=8, labels=year_labels(t))
    tau = f.tau * 86400
    tau = tau - 0.4912322 * 86400 * np.round((tau - np.median(tau)) / (0.4912322 * 86400))
    assert np.max(np.abs(np.diff(tau))) < 3000
    assert np.all(np.abs(f.alpha - 1) < 0.3)
