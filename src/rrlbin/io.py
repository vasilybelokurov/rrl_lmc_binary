"""I/O for OGLE-III / OGLE-IV LMC RR Lyrae catalogues and light curves.

Column layouts follow the README files of the two collections:
OGLE-IV (Soszynski+2016, Acta Astron. 66, 131) and OGLE-III (Soszynski+2009, Acta Astron. 59, 1).
Times in the light-curve files are HJD - 2450000.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


def _sexa_to_deg(ra: pd.Series, dec: pd.Series) -> tuple[np.ndarray, np.ndarray]:
    """'hh:mm:ss.ss', '+dd:mm:ss.s' -> degrees."""
    r = ra.str.split(":", expand=True).astype(float)
    d = dec.str.split(":", expand=True)
    sign = np.where(d[0].str.strip().str.startswith("-"), -1.0, 1.0)
    dd = d[0].str.replace("-", "").str.replace("+", "").astype(float)
    ra_deg = 15.0 * (r[0] + r[1] / 60 + r[2] / 3600)
    dec_deg = sign * (dd + d[1].astype(float) / 60 + d[2].astype(float) / 3600)
    return ra_deg.to_numpy(), np.asarray(dec_deg)


def _read_fwf(path, specs, names) -> pd.DataFrame:
    df = pd.read_fwf(path, colspecs=specs, names=names, header=None, dtype=str, keep_default_na=False)
    return df.apply(lambda c: c.str.strip())


def read_ogle4_ident(path) -> pd.DataFrame:
    """OGLE-IV ident.dat -> ogle_id, subtype, ra, dec [deg], ogle4_field_id, ogle3_id, ogle2_id, other_id."""
    specs = [(0, 20), (22, 26), (28, 39), (40, 51), (53, 69), (70, 85), (86, 101), (102, 300)]
    names = ["ogle_id", "subtype", "ra_s", "dec_s", "ogle4_field_id", "ogle3_id", "ogle2_id", "other_id"]
    df = _read_fwf(path, specs, names)
    df["ra"], df["dec"] = _sexa_to_deg(df.ra_s, df.dec_s)
    return df.drop(columns=["ra_s", "dec_s"])


def read_ogle3_ident(path) -> pd.DataFrame:
    """OGLE-III ident.dat -> ogle_id, field, dbnum, subtype, ra, dec [deg], ogle2_id, macho_id, gcvs_id, other_id."""
    specs = [(0, 20), (22, 30), (31, 37), (38, 43), (44, 55), (56, 67), (68, 83), (84, 98), (99, 104), (105, 300)]
    names = ["ogle_id", "field", "dbnum", "subtype", "ra_s", "dec_s", "ogle2_id", "macho_id", "gcvs_id", "other_id"]
    df = _read_fwf(path, specs, names)
    df["ra"], df["dec"] = _sexa_to_deg(df.ra_s, df.dec_s)
    return df.drop(columns=["ra_s", "dec_s"])


def read_lc(path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read an OGLE light-curve file: (t [HJD-2450000], mag, mag_err). Empty arrays for an empty file."""
    a = np.loadtxt(path, ndmin=2)
    if a.size == 0:
        return np.empty(0), np.empty(0), np.empty(0)
    return a[:, 0], a[:, 1], a[:, 2]


def lc_path(raw: str | Path, survey: str, ogle_id: str, band: str = "I") -> Path:
    """Path of a light curve in data/raw; survey in {'ogle3', 'ogle4'}."""
    return Path(raw) / f"{survey}_lmc_rrlyr" / "phot" / band / f"{ogle_id}.dat"


# ------------------------------------------------------------------ MACHO
def mjd_to_hjd(mjd, ra_deg: float, dec_deg: float) -> np.ndarray:
    """MJD (UTC) -> HJD - 2450000 (UTC), observer at Mount Stromlo (the MACHO site; topocentric term < 0.03 s).

    Note: whether MACHO `dateobs` is the exposure start or mid-point is not documented in the TAP
    metadata; a constant offset is calibrated empirically against OGLE-II (overlap 1997-1999).
    """
    from astropy import units as u
    from astropy.coordinates import EarthLocation, SkyCoord
    from astropy.time import Time

    stromlo = EarthLocation.from_geodetic(149.0089 * u.deg, -35.3207 * u.deg, 770 * u.m)
    tm = Time(np.asarray(mjd, float), format="mjd", scale="utc", location=stromlo)
    ltt = tm.light_travel_time(SkyCoord(ra_deg * u.deg, dec_deg * u.deg), kind="heliocentric")
    return (tm.mjd - 49999.5) + ltt.to_value(u.day)


def read_macho(tile_file, seqn: int, ra_deg: float, dec_deg: float, band: str = "b",
               max_err: float = 0.5) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """MACHO light curve of one star from a fetched tile file: (HJD-2450000, instrumental mag, err).

    band 'b' (MACHO blue) or 'r' (MACHO red). Rows with the -99 sentinels or err <= 0 or err > max_err are dropped.
    """
    d = pd.read_parquet(tile_file)
    d = d[d.seqn == seqn]
    m, e = d[f"{band}mag"].to_numpy(float), d[f"{band}err"].to_numpy(float)
    ok = (m > -50) & (e > 0) & (e < max_err)
    t = mjd_to_hjd(d.dateobs.to_numpy(float)[ok], ra_deg, dec_deg)
    o = np.argsort(t)
    return t[o], m[ok][o], e[ok][o]
