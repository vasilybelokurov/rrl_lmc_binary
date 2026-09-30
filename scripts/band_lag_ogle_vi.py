"""OGLE V vs I delay offset (same time system): the pure band-lag effect vs period, as a control for the
MACHO B - OGLE I offset. Delta_VI = weighted mean over common years of tau_V - tau_I (fundamental-phase gauge).

Usage
-----
    python scripts/band_lag_ogle_vi.py --n 600 --out results/macho/lag_ogle_vi.parquet
"""
from __future__ import annotations

import argparse
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from rrlbin.io import lc_path, read_lc  # noqa: E402
from rrlbin.timing import fit_timing, year_labels  # noqa: E402

RAW = Path("data/raw")
DAY = 86400.0


def band(oid, b):
    t3, m3, e3 = read_lc(lc_path(RAW, "ogle3", oid, b)) if lc_path(RAW, "ogle3", oid, b).exists() else (np.empty(0),) * 3
    t4, m4, e4 = read_lc(lc_path(RAW, "ogle4", oid, b))
    seg = np.r_[np.where(t3 < 2000, "O2", "O3"), np.full(t4.size, "O4")]
    return np.r_[t3, t4], np.r_[m3, m4], np.r_[e3, e4], seg


def one(row):
    oid, P, T0 = row
    try:
        tV, mV, eV, sV = band(oid, "V")
        if tV.size < 60:
            return None
        tI, mI, eI, sI = band(oid, "I")
        fV = fit_timing(tV, mV, eV, sV, P, T0, K=6, labels=year_labels(tV), min_season=8)
        fI = fit_timing(tI, mI, eI, sI, P, T0, K=8, labels=year_labels(tI))
        d, w = [], []
        for j, s in enumerate(fV.season):
            k = np.flatnonzero(fI.season == s)
            if k.size:
                d.append(fV.tau[j] - fI.tau[k[0]])
                w.append(1 / (fV.tau_err[j] ** 2 + fI.tau_err[k[0]] ** 2))
        if not d:
            return None
        d, w = np.array(d), np.array(w)
        d = d[0] + ((d - d[0] + P / 2) % P - P / 2)
        mean = np.sum(w * d) / w.sum()
        mean -= P * np.round(mean / P)
        return dict(ogle_id=oid, P=P, delta_s=mean * DAY, delta_err_s=w.sum() ** -0.5 * DAY, n=len(d))
    except Exception:
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=600)
    ap.add_argument("--out", default="results/macho/lag_ogle_vi.parquet")
    a = ap.parse_args()
    inv = pd.read_parquet("data/lc_inventory.parquet")
    par = pd.read_fwf(RAW / "ogle4_lmc_rrlyr" / "RRab.dat", colspecs=[(0, 20), (37, 46), (58, 68)],
                      names=["ogle_id", "P", "T0"], header=None)
    s = inv[(inv.subtype == "RRab") & (inv.n4_I > 0)].merge(par, on="ogle_id").sample(a.n, random_state=3)
    with Pool(6) as pool:
        rows = [r for r in pool.map(one, list(zip(s.ogle_id, s.P.astype(float), s.T0.astype(float))), chunksize=5) if r]
    df = pd.DataFrame(rows)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(a.out)
    g = df[df.delta_err_s < 300]
    c, cov = np.polyfit(g.P, g.delta_s, 1, w=1 / g.delta_err_s, cov=True)
    r = g.delta_s - np.polyval(c, g.P)
    x = g.delta_s / DAY / g.P
    print(f"N={len(g)}; V-I lag median {x.median():.4f} cycles (robust sd {1.4826 * np.median(np.abs(x - x.median())):.4f})")
    print(f"weighted fit Delta_VI = ({c[0]:.0f} +- {cov[0, 0] ** .5:.0f}) P + ({c[1]:.0f} +- {cov[1, 1] ** .5:.0f}) s; "
          f"robust resid sd {1.4826 * np.median(np.abs(r - np.median(r))):.0f} s, median err {g.delta_err_s.median():.0f} s")


if __name__ == "__main__":
    main()
