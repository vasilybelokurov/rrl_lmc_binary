"""Part C tests."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))


def test_cut_flags_toy():
    from plot_summary_stats import cut_flags
    df = pd.DataFrame(dict(D=[50, 50, 30, 50], alpha_chi2nu=[1, 3, 1, 1], amp_s=[800, 800, 800, 80000], err_med_s=[100] * 4,
                           P_best=[3000.0] * 4))
    f = cut_flags(df, "amp_s", 9000.0)
    assert list(f["all"]) == [True, False, False, False]


def test_v2_cuts_reproduce_v2_list():
    """Applying the same cut function to the v2 statistics reproduces the frozen v2 candidate list."""
    from plot_summary_stats import cut_flags
    v2 = pd.read_parquet("results/real/oc_all_cm.parquet").drop(columns=["t", "tau", "err", "flag", "alpha", "alpha_err"], errors="ignore")
    v2 = v2[v2.ok]
    sel = v2[cut_flags(v2, "amp_s", v2.baseline)["all"]]
    frozen = pd.read_csv("results/candidates/candidates_v2_corrgrid.csv")
    assert set(sel.ogle_id) == set(frozen.ogle_id)


def test_keplerian_fit_recovers_eccentric_orbits():
    """Injected Irwin orbits (e = 0, 0.4, 0.7) on a 25-season MACHO+OGLE-like cadence with a MACHO offset are recovered
    within ~2.5 bootstrap sigma (P, A, e)."""
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
    from rrlbin.kepler_fit import fit_keplerian
    from rrlbin.ltte import ltte_delay
    DAY = 86400.0
    rng = np.random.default_rng(3)
    t = np.array([245 + 365.25 * (y - 1998) + 180 + rng.normal(0, 40) for y in range(1992, 2017)])
    band = np.where(t < 1600, 1, 0)
    err = np.full(t.size, 120 / DAY)
    for e_true in (0.0, 0.4, 0.7):
        y = ltte_delay(t, 3200.0, 1200 / DAY, e_true, 1.0, 500.0) + 2e-10 * (t - 3000) ** 2 + np.where(band == 1, -1500 / DAY, 0) \
            + rng.normal(0, 1, t.size) * err
        f = fit_keplerian(t, y, err, band, 3000.0, priors={1: (-1500 / DAY, 300 / DAY)}, n_boot=60, seed=1)
        for k, truth in (("P", 3200.0), ("A_s", 1200.0), ("e", e_true)):
            lo, hi = f[f"{k}_p16"], f[f"{k}_p84"]
            half = max((hi - lo) / 2, {"P": 30, "A_s": 30, "e": 0.03}[k])
            assert abs(f[k] - truth) < 2.5 * half, (e_true, k, f[k], lo, hi)


def test_rv_from_delay_matches_keplerian_velocity():
    """c * d tau/dt equals the analytic RV curve K [cos(nu + w) + e cos w] with K = 2 pi c A / (P sqrt(1 - e^2))."""
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
    from partC_followup import rv_from_delay, C_KMS
    from rrlbin.ltte import kepler_E
    P, A, e, w, tp = 3000.0, 1000 / 86400.0, 0.5, 0.8, 200.0
    t = np.linspace(0, 6000, 500)
    M = 2 * np.pi * (t - tp) / P
    E = kepler_E(np.mod(M, 2 * np.pi), e)
    nu = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2), np.sqrt(1 - e) * np.cos(E / 2))
    K = 2 * np.pi * C_KMS * A / (P * np.sqrt(1 - e ** 2))
    analytic = K * (np.cos(nu + w) + e * np.cos(w))
    assert np.max(np.abs(rv_from_delay(t, P, A, e, w, tp, h=0.05) - analytic)) < 0.01 * K


def test_bootstrap_prediction_covers_truth():
    """On a simulated orbit, the 95% bootstrap envelope of the extrapolated (out-of-sample) O-C covers the true curve
    for most of the prediction window."""
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
    from rrlbin.kepler_fit import bootstrap_predictions
    from rrlbin.ltte import ltte_delay
    DAY = 86400.0
    rng = np.random.default_rng(5)
    t = np.array([245 + 365.25 * (y - 1998) + 180 + rng.normal(0, 40) for y in range(1992, 2016)])
    band = np.where(t < 1600, 1, 0)
    err = np.full(t.size, 150 / DAY)
    truth = lambda x: ltte_delay(x, 3000.0, 1000 / DAY, 0.3, 1.0, 400.0) + 1e-10 * (x - 3000) ** 2
    y = truth(t) + np.where(band == 1, -1500 / DAY, 0) + rng.normal(0, 1, t.size) * err
    tp = np.linspace(7600, 11000, 40)
    best, cur = bootstrap_predictions(t, y, err, band, 3000.0, tp, priors={1: (-1500 / DAY, 300 / DAY)}, n_boot=80)
    lo, hi = np.percentile(cur, [2.5, 97.5], axis=0)
    # compare shapes up to the free constant (the model's c0 is absorbed differently from truth's): remove the mean in-sample offset
    off = np.median((best - truth(tp))[:3])
    cover = np.mean((truth(tp) + off >= lo) & (truth(tp) + off <= hi))
    assert cover > 0.8


def test_krige_interpolates_and_widens():
    """Universal kriging: near the data the predictive sd is ~ the noise level; it grows with distance beyond the data;
    with A = 0 it reduces to the GLS quadratic prediction."""
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
    from rrlbin.predict import krige
    DAY = 86400.0
    rng = np.random.default_rng(2)
    t = np.sort(rng.uniform(0, 7000, 25))
    err = np.full(t.size, 100 / DAY)
    y = 1e-10 * (t - 3500) ** 2 + rng.normal(0, 1, t.size) * err
    tp = np.array([3500.0, 7500.0, 9000.0, 11000.0])
    m, c = krige(t, y, err, np.zeros(t.size, int), tp, 500 / DAY, 1500.0, 0.0)
    sd = np.sqrt(np.diag(c)) * DAY
    assert sd[0] < sd[1] < sd[2] < sd[3]
    m0, c0 = krige(t, y, err, np.zeros(t.size, int), tp, 0.0, 1500.0, 0.0)
    X = np.vander((t - t.mean()) / 1000, 3)
    b = np.linalg.lstsq(X, y, rcond=None)[0]
    assert np.allclose(m0, np.vander((tp - t.mean()) / 1000, 3) @ b, atol=1e-7)


def test_interp_prediction_and_gauss_score():
    """Interpolation of a tabulated predictive distribution is exact for a linear mean and reproduces the covariance of
    a smooth GP at grid points; gauss_score equals scipy's multivariate normal log-density, and data drawn from H1
    favour H1 over a distinct H0 (positive mean ln BF)."""
    import pytest
    from scipy.stats import multivariate_normal
    from rrlbin.predict import gauss_score, gp_kernel, interp_prediction
    tg = np.linspace(0, 3000, 101)
    mean, cov = 2.0 * tg + 5, gp_kernel(tg, tg, 300.0, 1500.0) + 1e-6 * np.eye(tg.size)
    tn = np.array([0.0, 17.0, 1234.5, 3000.0])
    m, C = interp_prediction(tg, mean, cov, tn)
    assert np.allclose(m, 2.0 * tn + 5)
    assert np.allclose(C, gp_kernel(tn, tn, 300.0, 1500.0), rtol=1e-3, atol=0)     # 1e-4 relative between grid points
    with pytest.raises(ValueError):
        interp_prediction(tg, mean, cov, [3001.0])
    C = C + np.eye(4) * 100.0 ** 2
    d = m + 50.0
    lnl, chi2, n = gauss_score(d, m, C)
    assert np.isclose(lnl, multivariate_normal(m, C).logpdf(d)) and n == 4
    rng = np.random.default_rng(3)
    m0 = m + 600.0
    bf = [gauss_score(x, m, C)[0] - gauss_score(x, m0, C)[0] for x in rng.multivariate_normal(m, C, 200)]
    assert np.mean(bf) > 3


def test_split_test_orbit_vs_red_noise():
    """split_test: an LTTE orbit (1000 s, 4000 d) plus white noise on a 1992-2026-like season grid is favoured out of sample
    (ln BF > 3, p(chi2|H1) > 0.01); smooth GP wander without an orbit is not favoured on average."""
    from rrlbin.ltte import ltte_delay
    from rrlbin.predict import gp_kernel, split_test
    rng = np.random.default_rng(8)
    t = np.r_[np.arange(-1000.0, 1700, 365.25), np.arange(2200.0, 7600, 365.25), np.arange(7800.0, 8900, 365.25),
              np.arange(9950.0, 11100, 365.25)]
    band = np.zeros(t.size, int)
    err = np.full(t.size, 120.0) / 86400
    P = 0.55
    y_orb = ltte_delay(t, 4000.0, 1000.0 / 86400, 0.2, 1.0, 1500.0) + 0.3e-3 * ((t - 5000) / 1000) ** 2 + rng.normal(0, 1, t.size) * err
    r = split_test(t, y_orb, err, band, P, 7600.0)
    assert r["ok"] and r["n_test"] == 8 and r["lnBF"] > 3 and r["p_h1"] > 0.01
    bf = []
    for k in range(6):
        K = gp_kernel(t, t, 900.0 / 86400, 1500.0) + np.eye(t.size) * err ** 2
        y = rng.multivariate_normal(np.zeros(t.size), K)
        bf.append(split_test(t, y, err, band, P, 7600.0)["lnBF"])
    assert np.median(bf) < 3
