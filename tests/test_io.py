"""I/O tests: time conversion."""
import numpy as np

from rrlbin.io import mjd_to_hjd


def test_mjd_to_hjd_lmc():
    """Toward the LMC (near the south ecliptic pole, beta ~ -85 deg) the heliocentric correction is small,
    |dt| < 499 s * cos(beta) ~ 45 s, and has a period of 1 yr."""
    mjd = 50500.0 + np.linspace(0, 365.25, 200)
    h = mjd_to_hjd(mjd, 80.9, -69.75)
    dt = (h - (mjd - 49999.5)) * 86400
    assert np.max(np.abs(dt)) < 60 and np.ptp(dt) > 40
    assert np.allclose(mjd_to_hjd([50500.0], 80.9, -69.75), mjd_to_hjd([50500.0 + 365.25], 80.9, -69.75) - 365.25,
                       atol=2e-5)
