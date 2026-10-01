"""Tests of the shared Level-2 O-C analysis (src/rrlbin/oc.py)."""
import numpy as np

from rrlbin.ltte import ltte_delay, oc_search
from rrlbin.ltte import period_grid as old_grid
from rrlbin.oc import (chi2nu_const_by_band, oc_stats, orbit_search, period_grid, robust_unwrap)
from rrlbin.simulate import delay_random_walk

P = 0.57
DAY = 86400.0


def seasons(rng, y0=1992, y1=2016):
    return np.array([245 + 365.25 * (y - 1998) + 180 + rng.normal(0, 40) for y in range(y0, y1 + 1)])


def test_robust_unwrap_repairs_single_slip_and_keeps_large_smooth_drift():
    rng = np.random.default_rng(0)
    t = seasons(rng)
    band = np.zeros(t.size, int)
    smooth = 1.5 * P * ((t - t.mean()) / np.ptp(t)) ** 2 * 4        # drift of ~1.5 cycles over the baseline
    wrapped = np.mod(smooth + P / 2, P) - P / 2                      # delays known only modulo P
    slip = wrapped.copy()
    slip[10] += P                                                   # one season off by a whole cycle
    for y in (wrapped, slip):
        u = robust_unwrap(t, y, band, P)
        d = (u - smooth) - np.median(u - smooth)
        assert np.max(np.abs(d)) < 1e-9


def test_orbit_search_matches_previous_implementation_for_ogle_only():
    rng = np.random.default_rng(1)
    t = seasons(rng, 2001, 2016)
    err = np.full(t.size, 250 / DAY)
    y = 700 / DAY * np.sin(2 * np.pi * t / 2500) + rng.normal(0, 1, t.size) * err
    per = old_grid(np.ptp(t))
    a = oc_search(t, y, err, per)
    b = orbit_search(t, y, err, np.zeros(t.size, int), per, 1)
    assert np.allclose(a["Dp"], b["Dp"], atol=1e-6)


def test_offset_prior_constrains_macho_without_overlap():
    """Without overlap years the MACHO offset is free; a prior on it lets MACHO pin the absolute level, raising D."""
    rng = np.random.default_rng(2)
    t = seasons(rng)
    band = np.where(t < 1500, 1, 0)
    band[(t > 1500) & (t < 2000)] = 0
    err = np.full(t.size, 200 / DAY)
    off = -1600 / DAY
    y = ltte_delay(t, 4000, 900 / DAY) + np.where(band == 1, off, 0) + rng.normal(0, 1, t.size) * err
    free = orbit_search(t, y, err, band, period_grid(np.ptp(t)), 1)
    pri = orbit_search(t, y, err, band, period_grid(np.ptp(t)), 1, priors={1: (off, 300 / DAY)})
    assert pri["D"] > free["D"]
    assert abs(pri["beta"][pri["names"].index("off1")] - off) * DAY < 400


def test_second_harmonic_recovers_eccentric_amplitude():
    rng = np.random.default_rng(3)
    t = np.sort(np.concatenate([seasons(rng), seasons(rng) + 0.3]))
    err = np.full(t.size, 100 / DAY)
    amp = 1000 / DAY
    y = ltte_delay(t, 3000, amp, e=0.6, omega=0.5, t_peri=100) + rng.normal(0, 1, t.size) * err
    s = oc_stats(dict(t=t, tau=y, err=err, band=np.zeros(t.size), P=P), with_red=False)
    true_half_ptp = 0.5 * np.ptp(ltte_delay(np.linspace(0, 3000, 5000), 3000, amp, e=0.6, omega=0.5)) * DAY
    assert abs(s["amp_2h_s"] - true_half_ptp) < abs(s["amp_circ_s"] - true_half_ptp)
    assert s["dD_harm2"] > 10


def test_alias_flag_and_predictive_score():
    rng = np.random.default_rng(4)
    t = seasons(rng)
    band = np.where(t < 1600, 1, 0)
    err = np.full(t.size, 150 / DAY)
    sc_orbit, sc_rw = [], []
    for _ in range(30):
        y = ltte_delay(t, rng.uniform(1500, 3500), 800 / DAY, t_peri=rng.uniform(0, 3000)) + rng.normal(0, 1, t.size) * err
        s = oc_stats(dict(t=t, tau=y, err=err, band=band, P=P), with_red=False)
        assert s["alias_dD"] > 0
        sc_orbit.append(s["pred_score"])
        y = delay_random_walk(t, 800 / DAY, rng) + rng.normal(0, 1, t.size) * err
        sc_rw.append(oc_stats(dict(t=t, tau=y, err=err, band=band, P=P), with_red=False)["pred_score"])
    assert np.median(sc_orbit) > np.median(sc_rw) + 1


def test_chi2nu_by_band():
    rng = np.random.default_rng(5)
    x = np.r_[rng.normal(1.0, 0.01, 10), rng.normal(1.3, 0.01, 10)]   # two bands with different means, each constant
    band = np.r_[np.zeros(10), np.ones(10)]
    assert 0.4 < chi2nu_const_by_band(x, np.full(20, 0.01), band) < 2.0


def test_grid_starts_at_400():
    assert abs(period_grid(8000).min() - 400) < 1


def test_common_mode_iterates_to_injected_pattern():
    """Inject a per-year common timing offset into 300 random-noise stars; the iterative estimate recovers it
    (up to the part a per-star quadratic + offset can absorb) and converges."""
    from rrlbin.oc import apply_common_mode, common_mode
    from rrlbin.timing import year_labels
    rng = np.random.default_rng(6)
    years = np.arange(-10, 25)
    true = {(0, int(y)): rng.normal(0, 60) / DAY for y in years}
    series = []
    for _ in range(300):
        t = seasons(rng)
        err = np.full(t.size, 150 / DAY)
        tau = rng.normal(0, 1, t.size) * err + np.array([true[(0, int(y))] for y in year_labels(t)])
        series.append(dict(t=t, tau=tau, err=err, band=np.zeros(t.size, int), P=P))
    cm, hist = common_mode(series, year_labels, n_iter=8)
    assert hist[-1] < 0.1 * hist[0] and hist[-1] * DAY < 3       # converges (degenerate directions projected out)
    # compare after removing the quadratic that the per-star H0 cannot distinguish from a common mode
    yv = np.array([k[1] for k in sorted(cm)])
    est = np.array([cm[k] for k in sorted(cm)])
    tru = np.array([true[k] for k in sorted(cm)])
    tru = tru - np.polyval(np.polyfit(yv, tru, 2), yv)       # the measurable part of the injected pattern
    d = est - tru
    assert np.std(d) * DAY < 15
    res_before = np.std([apply_common_mode(s, {}, year_labels)["tau"].std() for s in series])
    assert res_before >= 0


def test_robust_unwrap_clustered_slip_does_no_harm():
    """Two adjacent slipped seasons: the conservative rule may leave them, but must never increase the jumps."""
    rng = np.random.default_rng(8)
    t = seasons(rng)
    band = np.zeros(t.size, int)
    smooth = 0.3 * P * np.sin(2 * np.pi * t / 4000)
    y = np.mod(smooth + P / 2, P) - P / 2
    y[8:10] += P
    u = robust_unwrap(t, y, band, P)
    jump = lambda v: np.max(np.abs(np.diff(v)))
    seq = y[0] + np.unwrap(2 * np.pi * (y - y[0]) / P) * P / (2 * np.pi)
    assert jump(u) <= jump(seq) + 1e-12


def test_predictive_score_independent_of_held_out_wrapping():
    """Shifting a held-out season by a whole cycle must not change the training fit (no leakage via unwrapping)."""
    rng = np.random.default_rng(9)
    t = seasons(rng)
    band = np.where(t < 1600, 1, 0)
    band[(t > 400) & (t < 1600)] = 1
    tt = np.r_[t, t[(t > 450) & (t < 1700)] + 30]
    bb = np.r_[band, np.zeros(((t > 450) & (t < 1700)).sum(), int)]
    o = np.argsort(tt); tt, bb = tt[o], bb[o]
    err = np.full(tt.size, 150 / DAY)
    y = ltte_delay(tt, 2500, 800 / DAY) + rng.normal(0, 1, tt.size) * err
    a = oc_stats(dict(t=tt, tau=y, err=err, band=bb, P=P), with_red=False)
    y2 = y.copy()
    k = np.flatnonzero((bb == 1) & (tt < 450))[1]
    y2[k] += P
    b = oc_stats(dict(t=tt, tau=y2, err=err, band=bb, P=P), with_red=False)
    assert a["P_train"] == b["P_train"] and abs(a["D_train"] - b["D_train"]) < 1e-9


def test_orbit_amplitude_any_phase():
    """Regression test for a name collision (quadratic 'c1' vs orbit 'c1'): the fitted circular amplitude must not
    depend on the orbital phase."""
    rng = np.random.default_rng(10)
    t = seasons(rng)
    err = np.full(t.size, 30 / DAY)
    from rrlbin.oc import orbit_amplitude
    for phase in np.linspace(0, 2 * np.pi, 7):
        y = 1500 / DAY * np.sin(2 * np.pi * t / 3000 + phase) + rng.normal(0, 1, t.size) * err
        x = orbit_search(t, y, err, np.zeros(t.size, int), np.array([3000.0]), 1)
        assert abs(orbit_amplitude(x["beta"], x["names"], 3000.0, 1)[0] * DAY - 1500) < 60


def test_band_whole_cycle_shift_with_prior_is_harmless():
    """A MACHO band shifted by a whole pulsation cycle relative to its prior must give the same statistics
    (regression test: before align_bands, 13% of real MACHO stars had inflated jitter and D)."""
    rng = np.random.default_rng(13)
    t = seasons(rng)
    band = np.where(t < 1600, 1, 0)
    err = np.full(t.size, 150 / DAY)
    off = -1600 / DAY
    y = ltte_delay(t, 3000, 600 / DAY) + np.where(band == 1, off, 0) + rng.normal(0, 1, t.size) * err
    pri = {1: (off, 300 / DAY)}
    a = oc_stats(dict(t=t, tau=y, err=err, band=band, P=P), priors=pri, with_red=False)
    y2 = y + np.where(band == 1, P, 0.0)
    b = oc_stats(dict(t=t, tau=y2, err=err, band=band, P=P), priors=pri, with_red=False)
    assert abs(a["D"] - b["D"]) < 1e-6 and abs(a["jit0_s"] - b["jit0_s"]) < 1e-6


def test_gp_null_recovers_smooth_wander():
    """gp_null finds a large amplitude and a long correlation length for injected smooth wander, ~0 for white noise."""
    from rrlbin.oc import gp_null
    from rrlbin.simulate import delay_gp
    rng = np.random.default_rng(14)
    t = np.sort(np.concatenate([seasons(rng), seasons(rng) + 0.5]))
    err = np.full(t.size, 150 / DAY)
    band = np.zeros(t.size, int)
    A, ells = [], []
    for _ in range(10):
        y = delay_gp(t, 900 / DAY, 1500.0, rng) + rng.normal(0, 1, t.size) * err
        s, a, l = gp_null(t, y, err, band)
        A.append(a * DAY); ells.append(l)
    assert 500 < np.median(A) < 1600 and np.median(ells) >= 1000   # REML; amplitude on a coarse grid
    s, a, l = gp_null(t, rng.normal(0, 1, t.size) * err, err, band)
    assert a * DAY < 150
