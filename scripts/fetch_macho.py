"""Fetch MACHO light curves of OGLE LMC RR Lyrae through the NCI MACHO TAP service.

Stars are those with a MACHO ID (field.tile.seqn) in the OGLE-III ident.dat. One TAP query per (field, tile);
the result is saved as data/raw/macho/<field>.<tile>.parquet (all requested stars of that tile). Resumable:
tiles already on disk are skipped. Raw columns are kept (including -99 / 999 sentinels); cleaning and the
MJD -> HJD conversion are done downstream (rrlbin.io.read_macho).

Service: https://machotap.asvo.nci.org.au/ncitap/tap (table public.photometry_view; dateobs in MJD).

Usage
-----
    python scripts/fetch_macho.py --workers 4            # all tiles
    python scripts/fetch_macho.py --limit 5              # first 5 tiles (smoke test)
"""
from __future__ import annotations

import argparse
import io
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
import requests

TAP_URL = "https://machotap.asvo.nci.org.au/ncitap/tap"
COLS = ["fieldid", "tileid", "seqn", "obsid", "dateobs", "sideofpier", "exposure", "airmass",
        "rmag", "rerr", "rcrowd", "rtype", "rchi2", "rcosmicrf", "rfwhm", "rcut",
        "bmag", "berr", "bcrowd", "btype", "bchi2", "bcosmicrf", "bfwhm", "bcut"]
SELECT = ", ".join(COLS)


def tap_csv(query: str, timeout: float = 600) -> pd.DataFrame:
    """Synchronous TAP query returning CSV. (The VOTable route fails: astropy rejects the service's
    'character(1)' datatype of `sideofpier`, and the service's ADQL does not support CAST.)"""
    r = requests.post(f"{TAP_URL}/sync", data={"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv",
                                               "MAXREC": 10_000_000, "QUERY": query}, timeout=timeout)
    r.raise_for_status()
    if not r.text.startswith(COLS[0]):
        raise RuntimeError(r.text[:300])
    return pd.read_csv(io.StringIO(r.text))


def fetch_tile(field: int, tile: int, seqns: list[int], out: Path, retries: int = 4) -> int:
    q = (f"SELECT {SELECT} FROM public.photometry_view WHERE fieldid = {field} AND tileid = {tile} "
         f"AND seqn IN ({', '.join(map(str, seqns))})")
    for k in range(retries):
        try:
            df = tap_csv(q)
            tmp = out.with_suffix(".tmp")
            df.to_parquet(tmp)
            tmp.rename(out)
            return len(df)
        except Exception as e:  # network / service errors: back off and retry
            if k == retries - 1:
                raise RuntimeError(f"{field}.{tile}: {e}") from e
            time.sleep(5 * 2 ** k)
    return 0


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--inventory", default="data/lc_inventory.parquet")
    ap.add_argument("--out", default="data/raw/macho")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--limit", type=int, default=None)
    a = ap.parse_args()

    inv = pd.read_parquet(a.inventory)
    m = inv.macho_id.dropna()
    m = m[m.str.count(r"\.") == 2]
    ids = m.str.split(".", expand=True).astype(int)
    ids.columns = ["field", "tile", "seqn"]
    groups = sorted(ids.groupby(["field", "tile"]).seqn.apply(lambda s: sorted(set(s))).items())
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    todo = [(f, t, s) for (f, t), s in groups if not (out / f"{f}.{t}.parquet").exists()]
    if a.limit:
        todo = todo[:a.limit]
    print(f"stars {len(ids)}, tiles {len(groups)}, to fetch {len(todo)}", flush=True)

    t0, done, nrow, failed = time.time(), 0, 0, []
    with ThreadPoolExecutor(a.workers) as ex:
        futs = {ex.submit(fetch_tile, f, t, s, out / f"{f}.{t}.parquet"): (f, t) for f, t, s in todo}
        for fu in as_completed(futs):
            try:
                nrow += fu.result()
            except Exception as e:
                failed.append(futs[fu])
                print("FAILED", e, flush=True)
            done += 1
            if done % 50 == 0 or done == len(todo):
                dt = time.time() - t0
                print(f"{done}/{len(todo)} tiles, {nrow} rows, {dt:.0f} s, ETA {dt / done * (len(todo) - done):.0f} s",
                      flush=True)
    print(f"failed tiles: {len(failed)} {failed[:20]}")


if __name__ == "__main__":
    main()
