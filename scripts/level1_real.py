"""Level 1, real stars: light curves (OGLE I; MACHO B and R where available) -> timing fits -> season series.

Output (chunked, resumable): <out>/part_NNNN.parquet with one row per star: ogle_id, P, T0, chi2nu_I and the season-series
arrays t, tau, err, band, alpha, alpha_err, coh, coh_err (see rrlbin.oc). No O-C statistics here (that is Level 2).

Usage
-----
    python scripts/level1_real.py --out results/real/series_v3 [--limit 150 --require-o2] [--workers 6] [--ogle4 extended]
"""
from __future__ import annotations

import argparse
import sys

sys.path.insert(0, __import__("os").path.dirname(__file__))
import chunked  # noqa: E402,F401  (single-threaded BLAS)
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from rrlbin.io import read_ogle4_ident  # noqa: E402
from rrlbin.pipeline import star_series  # noqa: E402

KEYS = ["t", "tau", "err", "band", "alpha", "alpha_err", "coh", "coh_err"]


def sample_table():
    inv = pd.read_parquet("data/lc_inventory.parquet")
    ident = read_ogle4_ident("data/raw/ogle4_lmc_rrlyr/ident.dat")[["ogle_id", "ra", "dec"]]
    par = pd.read_fwf("data/raw/ogle4_lmc_rrlyr/RRab.dat", colspecs=[(0, 20), (37, 46), (58, 68)],
                      names=["ogle_id", "P", "T0"], header=None)
    s = inv[(inv.subtype == "RRab") & (inv.n3_I > 0) & (inv.n4_I > 0)].merge(par, on="ogle_id").merge(ident, on="ogle_id")
    return s.sort_values("ogle_id").reset_index(drop=True)


def one(row):
    oid, mid, P, T0, ra, dec, ogle4 = row
    try:
        s, fits, _ = star_series(oid, P, T0, mid, ra, dec, ogle4=ogle4)
        # usable only with an OGLE I fit (the same requirement as in the simulations)
        out = dict(ogle_id=oid, ok="I" in fits, P=P, T0=T0, chi2nu_I=s["chi2nu_I"], bands="".join(sorted(fits)),
                   failures="; ".join(f"{k}: {v}" for k, v in s["failures"].items()))
        out.update({k: s[k] for k in KEYS})
        return out
    except Exception as ex:
        return dict(ogle_id=oid, ok=False, msg=str(ex)[:200])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="results/real/series_v3")
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--require-o2", action="store_true")
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--ogle4", choices=["public", "extended"], default="public",
                    help="OGLE-IV light curves: public OCVS (to 2016) or the extended 2010-2026 files")
    a = ap.parse_args()
    s = sample_table()
    if a.require_o2:
        s = s[(s.t3_first < 2000) & s.macho_id.fillna("").str.count(r"\.").eq(2)]
    if a.limit:
        s = s.sample(a.limit, random_state=0).sort_values("ogle_id")
    jobs = [j + (a.ogle4,) for j in zip(s.ogle_id, s.macho_id, s.P.astype(float), s.T0.astype(float), s.ra, s.dec)]
    print(f"stars: {len(jobs)}", flush=True)
    chunked.run_chunked(one, jobs, a.out, workers=a.workers, chunk=250,
                        manifest=dict(script="level1_real", ids=[j[0] for j in jobs], require_o2=a.require_o2, ogle4=a.ogle4))
    d = chunked.merge_parts(a.out)
    print(f"ok {d.ok.sum()} / {len(d)}; bands: {d[d.ok].bands.value_counts().to_dict()}")
    if "failures" in d:
        f = d.failures.fillna("")
        print("band-fit failures:", f[f != ""].str.split(";").explode().str.split(":").str[0].str.strip().value_counts().to_dict())


if __name__ == "__main__":
    main()
