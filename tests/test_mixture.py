"""Mixture-fit tests on synthetic Gaussian classes (known densities)."""
import numpy as np

from rrlbin.mixture import ClassDensities, fit_fractions


def test_em_recovers_fractions_with_known_densities():
    rng = np.random.default_rng(0)
    means = {"a": [0, 0], "b": [2, 0], "c": [0, 3]}
    sims = {k: rng.normal(m, 1, (4000, 2)) for k, m in means.items()}
    dens = ClassDensities(sims, jitter=0.0)
    f_true = np.array([0.6, 0.3, 0.1])
    n = rng.multinomial(5000, f_true)
    X = np.vstack([rng.normal(means[k], 1, (m, 2)) for k, m in zip(means, n)])
    f = fit_fractions(dens.logpdf(X))
    assert np.allclose(f, f_true, atol=0.03)
    assert np.isclose(f.sum(), 1)
