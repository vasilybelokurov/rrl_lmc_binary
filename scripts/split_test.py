"""Out-of-sample orbit test (rrlbin.predict.split_test) for candidates, applied identically to real stars and simulations.

Train on all seasons before --t-split (default HJD' 7600, i.e. 1992-2016.5), predict the OGLE I seasons after it (2016/17-2025/26):
orbit (Keplerian + red noise refitted after the orbit) vs no orbit (red noise), all fitted on the training seasons only.
Selection: the v4 candidate cuts (plot_summary_stats.cut_flags on the Level-2 stats) and/or an explicit --ids list. For the
simulations the same cuts select the rows, so the pass rate of the test among SELECTED nuisance sims (null, Blazhko, jumps,
random walks, empirical red noise) is the false-pass rate including selection.

Usage
-----
    python scripts/split_test.py --series results/real/series_v4 --stats results/real/stats_v4.parquet \
        --cm results/calib/common_mode_v4.json --lag results/calib/band_lag_v4.json --extra-ids results/predictions/test_2026-10-06.csv \
        --out results/split_v4/real.csv --workers 2
    python scripts/split_test.py --series results/inject/series_v4_macho --stats results/inject/stats_v4_macho.parquet \
        --lag results/calib/band_lag_v4.json --out results/split_v4/sims_macho.csv --workers 6
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
from level2 import KEYS, load, priors_for, to_series  # noqa: E402
from plot_summary_stats import cut_flags  # noqa: E402
from rrlbin.oc import apply_common_mode  # noqa: E402
from rrlbin.predict import split_test  # noqa: E402
from rrlbin.timing import year_labels  # noqa: E402

DAY = 86400.0
G = {}


def work(job):
    i, rec = job
    try:
        s = to_series(rec)
        if G["cm"]:
            s = apply_common_mode(s, G["cm"], year_labels)
        r = split_test(s["t"], s["tau"], s["err"], s["band"], s["P"], G["t_split"], priors=priors_for(s["P"], G["lag"]))
        return dict(row=i, **r)
    except Exception as ex:
        return dict(row=i, ok=False, err_msg=str(ex)[:200])


def _init(cm, lag, t_split):
    G.update(cm=cm, lag=lag, t_split=t_split)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--series", required=True)
    ap.add_argument("--stats", required=True)
    ap.add_argument("--cm", default=None)
    ap.add_argument("--lag", default=None)
    ap.add_argument("--t-split", type=float, default=7600.0)
    ap.add_argument("--extra-ids", default=None, help="CSV with ogle_id: also test these stars (real data), whether or not selected")
    ap.add_argument("--max-per-kind", type=int, default=None, help="simulations: cap the number of selected rows per class")
    ap.add_argument("--out", required=True)
    ap.add_argument("--workers", type=int, default=2)
    a = ap.parse_args()
    # the same row filter as level2.py, so that series row i <-> stats row i
    d = load(a.series)
    if "ok" in d:
        d = d[d.ok.fillna(True).astype(bool)]
    if "kind" in d:
        d = d[d.kind != "error"]
    d = d[d.t.apply(len) >= 6].reset_index(drop=True)
    st = pd.read_parquet(a.stats)
    if len(st) != len(d) or not (st.ogle_id.values == d.ogle_id.values).all():
        raise SystemExit("series and stats rows do not align")
    f = st.join(cut_flags(st, "amp_s", st.baseline))
    sel = f["all"].fillna(False).to_numpy(bool).copy()
    if a.extra_ids:
        sel |= d.ogle_id.isin(set(pd.read_csv(a.extra_ids).ogle_id)).to_numpy()
    idx = np.flatnonzero(sel)
    if a.max_per_kind and "kind" in d:
        rng = np.random.default_rng(1)
        keep = []
        for k, g in pd.Series(idx).groupby(d.kind.to_numpy()[idx]):
            keep += list(rng.permutation(g.to_numpy())[: a.max_per_kind])
        idx = np.sort(np.array(keep, int))
    print(f"rows {len(d)}; selected {len(idx)}" + (f"; by class {d.kind.iloc[idx].value_counts().to_dict()}" if "kind" in d else ""), flush=True)
    cm = {tuple(map(int, k.split(","))): v / DAY for k, v in json.loads(Path(a.cm).read_text()).items()} if a.cm else {}
    lag = json.loads(Path(a.lag).read_text()) if a.lag else None
    recs = d.to_dict("records")
    with Pool(a.workers, initializer=_init, initargs=(cm, lag, a.t_split)) as pool:
        rows = pool.map(work, [(i, recs[i]) for i in idx], chunksize=1)
    R = pd.DataFrame(rows).set_index("row")
    meta = d.drop(columns=[k for k in KEYS if k in d]).iloc[R.index]
    out = pd.concat([meta, f.iloc[R.index][["D", "P_best", "amp_s", "all"]].rename(columns={"all": "selected"}), R], axis=1)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(a.out, index=False)
    ok = out[out.ok.fillna(False).astype(bool)]
    pas = (ok.p_h1 > 0.01) & (ok.lnBF > 3)
    print(f"tested {len(ok)} / {len(out)}; pass (p_h1 > 0.01 and lnBF > 3): {int(pas.sum())}")
    if "kind" in ok:
        print(ok.assign(pas=pas).groupby("kind").pas.agg(["mean", "sum", "size"]).round(3).to_string())


if __name__ == "__main__":
    main()
