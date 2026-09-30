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
