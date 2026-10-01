"""Level-1 pipeline shared by real data and simulations: light curves -> per-band timing fits -> season series.

Bands: 0 = OGLE I (segments O2/O3/O4 with separate zero points, one template), 1 = MACHO B, 2 = MACHO R (own templates).
Seasons are observing years (timing.year_labels). A band needs >= MIN_EPOCHS epochs to be used.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from .io import lc_path, read_lc, read_macho
from .timing import fit_timing, year_labels

MIN_EPOCHS = 100
RAW = Path("data/raw")


def load_star(ogle_id: str, macho_id=None, ra=None, dec=None, raw=RAW) -> dict:
    """Light curves: {'I': (t, m, e, seg), 'MB': (t, m, e, seg), 'MR': (t, m, e, seg)} (MACHO bands only if available)."""
    out = {}
    p3, p4 = lc_path(raw, "ogle3", ogle_id), lc_path(raw, "ogle4", ogle_id)
    t3, m3, e3 = read_lc(p3) if p3.exists() else (np.empty(0),) * 3
    t4, m4, e4 = read_lc(p4) if p4.exists() else (np.empty(0),) * 3
    seg = np.r_[np.where(t3 < 2000, "O2", "O3"), np.full(t4.size, "O4")]
    out["I"] = (np.r_[t3, t4], np.r_[m3, m4], np.r_[e3, e4], seg)
    if isinstance(macho_id, str) and macho_id.count(".") == 2:
        f_, t_, s_ = map(int, macho_id.split("."))
        tile = Path(raw) / "macho" / f"{f_}.{t_}.parquet"
        if tile.exists():
            for band in ("b", "r"):
                t, m, e = read_macho(tile, s_, ra, dec, band=band)
                if t.size >= MIN_EPOCHS:
                    out["M" + band.upper()] = (t, m, e, np.full(t.size, "M" + band.upper()))
    return out


def fit_star(lcs: dict, P: float, T0: float, K: int = 8, failures: dict | None = None) -> dict:
    """Timing fit per band (year seasons). Bands with too few epochs or failing fits are skipped; the reason is recorded
    in `failures` (band -> message) if a dict is given (same rule for real stars and simulations)."""
    fits = {}
    for b, (t, m, e, seg) in lcs.items():
        if t.size < MIN_EPOCHS:
            if failures is not None:
                failures[b] = f"only {t.size} epochs"
            continue
        try:
            fits[b] = fit_timing(t, m, e, seg, P, T0, K=K, labels=year_labels(t))
        except Exception as ex:
            if failures is not None:
                failures[b] = str(ex)[:100]
    return fits


BAND_CODE = {"I": 0, "MB": 1, "MR": 2}


def series_from_fits(fits: dict, P: float) -> dict:
    """Concatenate the per-band season results into one series dict (see rrlbin.oc)."""
    keys = ["t", "tau", "err", "band", "alpha", "alpha_err", "coh", "coh_err"]
    s = {k: [] for k in keys}
    for b, f in fits.items():
        n = f.season.size
        s["t"].append(f.t_season)
        s["tau"].append(f.tau)
        s["err"].append(f.tau_err)
        s["band"].append(np.full(n, BAND_CODE[b]))
        s["alpha"].append(f.alpha)
        s["alpha_err"].append(f.alpha_err)
        s["coh"].append(f.coh)
        s["coh_err"].append(f.coh_err)
    out = {k: (np.concatenate(v) if v else np.empty(0)) for k, v in s.items()}
    out["P"] = P
    out["chi2nu_I"] = fits["I"].chi2nu if "I" in fits else np.nan
    return out


def star_series(ogle_id, P, T0, macho_id=None, ra=None, dec=None, raw=RAW):
    """Real star: (series, fits, light curves)."""
    lcs = load_star(ogle_id, macho_id, ra, dec, raw)
    failures = {}
    fits = fit_star(lcs, P, T0, failures=failures)
    s = series_from_fits(fits, P)
    s["failures"] = failures
    return s, fits, lcs
