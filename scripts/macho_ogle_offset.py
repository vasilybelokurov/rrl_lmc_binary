"""Calibrate the MACHO-vs-OGLE delay offset on stars observed by both MACHO (1992-1999) and OGLE-II (1997-2000).

Separate band templates (MACHO B or R vs OGLE I) leave an arbitrary constant delay offset per star. In seasons
observed by both surveys the true delay is common, so Delta = tau_MACHO - tau_OGLE measures that offset
(band-dependent template phase difference + any time-system offset, e.g. exposure start vs mid-point).
Per star: weighted mean Delta and chi2 of the season differences about it (consistency check).

Usage
-----
    python scripts/macho_ogle_offset.py --band b --out results/macho/offset_b.parquet
"""
from __future__ import annotations

import argparse
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from rrlbin.io import read_macho, read_ogle4_ident  # noqa: E402
from rrlbin.timing import fit_timing, harmonic_amp_phase, year_labels  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from timing_sample import ogle_I  # noqa: E402

DAY = 86400.0
BAND = "b"


def one(row):
    oid, mid, P, T0, ra, dec = row
    f_, t_, s_ = map(int, mid.split("."))
    tile = Path(f"data/raw/macho/{f_}.{t_}.parquet")
    if not tile.exists():
        return None
    try:
        tM, mM, eM = read_macho(tile, s_, ra, dec, band=BAND)
        if tM.size < 100:
            return dict(ogle_id=oid, ok=False, why="few MACHO epochs")
        fM = fit_timing(tM, mM, eM, np.full(tM.size, "M"), P, T0, K=8, labels=year_labels(tM))
        t, m, e, seg = ogle_I(oid)
        fO = fit_timing(t, m, e, seg, P, T0, K=8, labels=year_labels(t))
        # the template phase is arbitrary per band: wrap both delay series to within P/2 of their median
        dm = []
        for j, tj in enumerate(fM.t_season):
            k = np.flatnonzero(fO.season == fM.season[j])
            if k.size:
                k = k[0]
                d = fM.tau[j] - fO.tau[k]
                dm.append((tj, d, np.hypot(fM.tau_err[j], fO.tau_err[k])))
        if len(dm) == 0:
            return dict(ogle_id=oid, ok=False, why="no overlap")
        tj, d, s = map(np.array, zip(*dm))
        d = d[0] + ((d - d[0] + P / 2) % P - P / 2)          # remove whole-cycle ambiguities between seasons
        d = d - P * np.round(np.median(d) / P)               # and bring the offset into (-P/2, P/2]
        w = s ** -2
        mean = np.sum(w * d) / w.sum()
        A, ph = harmonic_amp_phase(fO.coef)
        AM, phM = harmonic_amp_phase(fM.coef)
        return dict(ogle_id=oid, ok=True, n_overlap=len(d), delta_s=mean * DAY, delta_err_s=w.sum() ** -0.5 * DAY,
                    chi2nu=float(np.sum(w * (d - mean) ** 2) / max(len(d) - 1, 1)), P=P,
                    nM=int(fM.mask.sum()), errM_med_s=float(np.median(fM.tau_err) * DAY),
                    errO_med_s=float(np.median(fO.tau_err) * DAY), A1_I=A[0], A1_M=AM[0],
                    dphi1=float(np.mod(phM[0] - ph[0] + np.pi, 2 * np.pi) - np.pi))
    except Exception as ex:
        return dict(ogle_id=oid, ok=False, why=str(ex)[:120])


def main() -> None:
    global BAND
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--band", default="b")
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--max-stars", type=int, default=None)
    ap.add_argument("--out", default="results/macho/offset_b.parquet")
    a = ap.parse_args()
    BAND = a.band
    inv = pd.read_parquet("data/lc_inventory.parquet")
    ident = read_ogle4_ident("data/raw/ogle4_lmc_rrlyr/ident.dat")[["ogle_id", "ra", "dec"]]
    par = pd.read_fwf("data/raw/ogle4_lmc_rrlyr/RRab.dat", colspecs=[(0, 20), (37, 46), (58, 68)],
                      names=["ogle_id", "P", "T0"], header=None)
    s = inv[(inv.subtype == "RRab") & (inv.t3_first < 2000) & (inv.macho_id.fillna("").str.count(r"\.") == 2)]
    s = s.merge(par, on="ogle_id").merge(ident, on="ogle_id")
    if a.max_stars:
        s = s.head(a.max_stars)
    jobs = list(zip(s.ogle_id, s.macho_id, s.P.astype(float), s.T0.astype(float), s.ra, s.dec))
    with Pool(a.workers, initializer=_set_band, initargs=(BAND,)) as pool:
        rows = [r for r in pool.map(one, jobs, chunksize=4) if r is not None]
    df = pd.DataFrame(rows)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(a.out)
    ok = df[df.ok]
    print(f"candidates {len(s)}, with tile {len(df)}, ok {len(ok)}; failures: {df[~df.ok].why.value_counts().head().to_dict()}")
    if len(ok):
        q = lambda c: np.percentile(ok[c], [10, 50, 90]).round(2)
        for c in ["n_overlap", "nM", "errM_med_s", "errO_med_s", "delta_s", "delta_err_s", "chi2nu", "dphi1"]:
            print(f"{c:12s} p10/50/90 {q(c)}")
        good = ok[ok.delta_err_s < 150]
        print(f"Delta (err < 150 s, N={len(good)}): weighted mean {np.average(good.delta_s, weights=good.delta_err_s ** -2):.1f} s,"
              f" robust scatter {1.4826 * np.median(np.abs(good.delta_s - good.delta_s.median())):.1f} s,"
              f" median err {good.delta_err_s.median():.1f} s")


def _set_band(b):
    global BAND
    BAND = b


if __name__ == "__main__":
    main()
