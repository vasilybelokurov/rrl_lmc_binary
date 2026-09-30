"""Chunked, resumable parallel map for long runs: results are written every `chunk` jobs to <out_dir>/part_NNNN.parquet
(existing parts are skipped), with a progress line per chunk. Use merge_parts() to combine.

Set single-threaded BLAS before numpy is imported in the workers (avoids oversubscription with many processes):
OMP_NUM_THREADS / OPENBLAS_NUM_THREADS / VECLIB_MAXIMUM_THREADS = 1.
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import time  # noqa: E402
from multiprocessing import Pool  # noqa: E402
from pathlib import Path  # noqa: E402

import pandas as pd  # noqa: E402


def run_chunked(func, jobs, out_dir, workers=6, chunk=500, flatten=False):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    n_chunks = (len(jobs) + chunk - 1) // chunk
    t0 = time.time()
    with Pool(workers) as pool:
        for c in range(n_chunks):
            f = out / f"part_{c:04d}.parquet"
            if f.exists():
                continue
            res = pool.map(func, jobs[c * chunk:(c + 1) * chunk], chunksize=4)
            rows = [r for rs in res for r in rs] if flatten else [r for r in res if r is not None]
            pd.DataFrame(rows).to_parquet(f)
            dt = time.time() - t0
            print(f"chunk {c + 1}/{n_chunks} done, {dt:.0f} s elapsed", flush=True)


def merge_parts(out_dir) -> pd.DataFrame:
    return pd.concat([pd.read_parquet(f) for f in sorted(Path(out_dir).glob("part_*.parquet"))], ignore_index=True)
