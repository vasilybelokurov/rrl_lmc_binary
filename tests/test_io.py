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


def test_load_star_extended_ogle4():
    """The extended OGLE-IV files (2010-2026) contain the public epochs plus later ones; in the high-cadence fields
    (LMC502/503/509/510/511/516; 10449 is in one) the photometry is re-reduced at the few-mmag level, elsewhere identical.
    The default loader is unchanged (public). Skipped if the data are not on disk."""
    import pytest
    from rrlbin.io import lc_path
    from rrlbin.pipeline import RAW, load_star

    oid = "OGLE-LMC-RRLYR-10449"
    if not lc_path(RAW, "ogle4x", oid).exists():
        pytest.skip("extended OGLE-IV data not available")
    pub, ext = load_star(oid)["I"], load_star(oid, ogle4="extended")["I"]
    t_pub = pub[0][pub[3] == "O4"]
    t_ext, m_ext = ext[0][ext[3] == "O4"], ext[1][ext[3] == "O4"]
    assert np.array_equal(pub[0][pub[3] != "O4"], ext[0][ext[3] != "O4"])         # OGLE-II/III untouched
    common, i_pub, i_ext = np.intersect1d(np.round(t_pub, 5), np.round(t_ext, 5), return_indices=True)
    assert common.size >= 0.99 * t_pub.size
    dm = pub[1][pub[3] == "O4"][i_pub] - m_ext[i_ext]
    assert np.median(np.abs(dm)) < 0.005 and np.abs(dm).max() < 0.05
    assert t_ext.max() > 11000 > t_pub.max()                                       # extends to 2026
    with pytest.raises(ValueError):
        load_star(oid, ogle4="bogus")
