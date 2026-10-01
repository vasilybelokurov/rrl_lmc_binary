"""Level 2: O-C statistics (rrlbin.oc.oc_stats) for saved season series (real or simulated), identical code for both.

Steps
  1. (real only, --common-mode out.json) iterative per-(band, year) common-mode estimate from well-behaved stars
     (>= 12 seasons), saved; or --apply-cm file.json to apply a saved one (sims: none, by construction).
  2. (optional, --band-lag json) Gaussian priors on the MACHO B/R offsets: mean = slope P + intercept, sd = scatter
     (from scripts/calibrate_band_lag.py).
  3. oc_stats per series, in parallel.

Usage
-----
    python scripts/level2.py --series results/real/series_v3 --common-mode results/calib/common_mode_v3.json \
        --band-lag results/calib/band_lag_v3.json --out results/real/stats_v3.parquet
    python scripts/level2.py --series results/inject/series_v3 --band-lag results/calib/band_lag_v3.json \
        --out results/inject/stats_v3.parquet
"""
from __future__ import annotations

import argparse
import json
import sys

sys.path.insert(0, __import__("os").path.dirname(__file__))
import chunked  # noqa: E402,F401  (single-threaded BLAS)
from multiprocessing import Pool  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from rrlbin.oc import apply_common_mode, common_mode, oc_stats  # noqa: E402
from rrlbin.timing import year_labels  # noqa: E402

DAY = 86400.0
KEYS = ["t", "tau", "err", "band", "alpha", "alpha_err", "coh", "coh_err"]
G = {}


def load(src):
    p = Path(src)
    return chunked.merge_parts(p) if p.is_dir() else pd.read_parquet(p)


def to_series(r):
    s = {k: np.asarray(r[k], float) for k in KEYS}
    s["P"] = float(r["P"])
    return s


def priors_for(P, lag):
    if not lag:
        return None
    return {int(b): ((v["slope"] * P + v["intercept"]) / DAY, v["sd"] / DAY) for b, v in lag.items()}


def work(r):
    try:
        s = to_series(r)
        if G["cm"]:
            s = apply_common_mode(s, G["cm"], year_labels)
        st = oc_stats(s, priors=priors_for(s["P"], G["lag"]), with_red=G["red"])
        st.pop("tau_unwrapped", None)
        return st
    except Exception as ex:
        return dict(err_msg=str(ex)[:200])


def _init(cm, lag, red):
    G.update(cm=cm, lag=lag, red=red)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--series", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--common-mode", default=None, help="estimate the common mode from this population and save it here")
    ap.add_argument("--apply-cm", default=None, help="apply a saved common mode")
    ap.add_argument("--band-lag", default=None)
    ap.add_argument("--no-red", action="store_true")
    ap.add_argument("--cm-exclude", default=None, help="CSV with ogle_id: stars excluded from the common-mode estimate")
    ap.add_argument("--workers", type=int, default=6)
    a = ap.parse_args()
    d = load(a.series)
    if "ok" in d:
        d = d[d.ok.fillna(True).astype(bool)]
    if "kind" in d:
        d = d[d.kind != "error"]
    d = d[d.t.apply(len) >= 6].reset_index(drop=True)
    cm = {}
    if a.common_mode:
        series = [to_series(r) for r in d.to_dict("records")]
        ex = set(pd.read_csv(a.cm_exclude).ogle_id) if a.cm_exclude else set()
        sel = [(len(s["t"]) >= 12) and (oid not in ex) for s, oid in zip(series, d.ogle_id)]
        cm, hist = common_mode(series, year_labels, n_iter=6, select=sel)
        Path(a.common_mode).parent.mkdir(parents=True, exist_ok=True)
        Path(a.common_mode).write_text(json.dumps({f"{b},{y}": v * DAY for (b, y), v in sorted(cm.items())}, indent=1))
        print("common mode [s]:", {k: round(v * DAY, 1) for k, v in sorted(cm.items())}, "; max update per iteration [s]:",
              [round(h * DAY, 2) for h in hist], flush=True)
    elif a.apply_cm:
        cm = {tuple(map(int, k.split(","))): v / DAY for k, v in json.loads(Path(a.apply_cm).read_text()).items()}
    lag = json.loads(Path(a.band_lag).read_text()) if a.band_lag else None
    recs = d.to_dict("records")
    with Pool(a.workers, initializer=_init, initargs=(cm, lag, not a.no_red)) as pool:
        rows = pool.map(work, recs, chunksize=20)
    st = pd.DataFrame(rows)
    meta = d.drop(columns=[k for k in KEYS if k in d])
    out = pd.concat([meta.reset_index(drop=True), st], axis=1)
    out = out.loc[:, ~out.columns.duplicated()]
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_parquet(a.out)
    print(f"written {a.out}: {len(out)} rows; failures {out.err_msg.notna().sum() if 'err_msg' in out else 0}")


if __name__ == "__main__":
    main()
