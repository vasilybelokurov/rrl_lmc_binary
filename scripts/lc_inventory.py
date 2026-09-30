"""Inventory of OGLE-III and OGLE-IV LMC RR Lyrae I-band light curves.

For every star in the OGLE-IV collection: number of I-band epochs and time span in OGLE-III and OGLE-IV,
and the combined baseline. OGLE-III light curves are linked through the OGLE-III ID in the OGLE-IV
ident.dat (the OGLE-III and OGLE-IV collections use the same OGLE-LMC-RRLYR-NNNNN numbering for
stars already in OGLE-III; this is checked below by comparing coordinates).

Usage
-----
    python scripts/lc_inventory.py --raw data/raw --out data/lc_inventory.parquet
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from rrlbin.io import read_lc, read_ogle3_ident, read_ogle4_ident  # noqa: E402


def summarize(path: Path) -> tuple[int, float, float]:
    """Return (N epochs, first HJD-2450000, last HJD-2450000); (0, nan, nan) if the file is missing."""
    if not path.exists():
        return 0, np.nan, np.nan
    t = read_lc(path)[0]
    if t.size == 0:
        return 0, np.nan, np.nan
    return t.size, float(t.min()), float(t.max())


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--raw", default="data/raw")
    ap.add_argument("--out", default="data/lc_inventory.parquet")
    a = ap.parse_args()
    raw = Path(a.raw)
    o4 = read_ogle4_ident(raw / "ogle4_lmc_rrlyr" / "ident.dat")
    o3 = read_ogle3_ident(raw / "ogle3_lmc_rrlyr" / "ident.dat")

    # Consistency check of the shared numbering: same ID -> same position.
    m = o4.merge(o3[["ogle_id", "ra", "dec"]], on="ogle_id", suffixes=("", "_o3"))
    sep = 3600 * np.hypot((m.ra - m.ra_o3) * np.cos(np.radians(m.dec)), m.dec - m.dec_o3)
    print(f"IDs in both collections: {len(m)}; separation median {np.median(sep):.3f}\", "
          f"max {sep.max():.2f}\", >1\": {(sep > 1).sum()}")

    rows = []
    for oid in o4.ogle_id:
        n4, a4, b4 = summarize(raw / "ogle4_lmc_rrlyr" / "phot" / "I" / f"{oid}.dat")
        n3, a3, b3 = summarize(raw / "ogle3_lmc_rrlyr" / "phot" / "I" / f"{oid}.dat")
        rows.append((oid, n3, a3, b3, n4, a4, b4))
    inv = pd.DataFrame(rows, columns=["ogle_id", "n3_I", "t3_first", "t3_last", "n4_I", "t4_first", "t4_last"])
    inv = o4[["ogle_id", "subtype", "ogle3_id", "ogle2_id"]].merge(inv, on="ogle_id")
    inv = inv.merge(o3[["ogle_id", "macho_id"]], on="ogle_id", how="left")
    inv["t_first"] = inv[["t3_first", "t4_first"]].min(axis=1)
    inv["t_last"] = inv[["t3_last", "t4_last"]].max(axis=1)
    inv["baseline_d"] = inv.t_last - inv.t_first
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    inv.to_parquet(a.out)

    has3, has4 = inv.n3_I > 0, inv.n4_I > 0
    print(f"OGLE-IV stars: {len(inv)}; with O-IV I LC: {has4.sum()}; with O-III I LC: {has3.sum()}; both: {(has3 & has4).sum()}")
    print(f"with MACHO ID (via OGLE-III ident): {(inv.macho_id.fillna('').str.len() > 0).sum()}; "
          f"with OGLE-II ID: {(inv.ogle2_id.fillna('').str.len() > 0).sum()}")
    for st in ["RRab", "RRc", "RRd"]:
        s = inv[(inv.subtype == st) & has3 & has4]
        print(f"{st} both: N={len(s)}; median n3={s.n3_I.median():.0f}, n4={s.n4_I.median():.0f}; "
              f"baseline median {s.baseline_d.median():.0f} d, p10 {s.baseline_d.quantile(.1):.0f}, p90 {s.baseline_d.quantile(.9):.0f}")
    s = inv[has4 & ~has3]
    print(f"O-IV only: N={len(s)}; median n4={s.n4_I.median():.0f}; baseline median {s.baseline_d.median():.0f} d")
    print(f"O-III: t range {inv.t3_first.min():.1f}-{inv.t3_last.max():.1f}; O-IV: {inv.t4_first.min():.1f}-{inv.t4_last.max():.1f}")


if __name__ == "__main__":
    main()
