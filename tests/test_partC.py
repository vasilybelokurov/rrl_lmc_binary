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
