"""Is the excess season-to-season O-C scatter of real stars intrinsic timing noise (common to all bands) or band-specific
(e.g. correlated photometric noise)? For stars with OGLE I, MACHO B and MACHO R, take each band's residuals about a
common H0 fit (quadratic + band offsets) and correlate the normalized residuals of the SAME season:
  B vs R (same MACHO nights, different filters), and MACHO vs OGLE (independent telescopes) in 1997-2000.
Intrinsic timing noise -> positive correlation in both; photometric noise of the MACHO system -> B-R only;
independent white noise -> neither.

Usage
-----
    python scripts/noise_origin_test.py --series results/validation_v3/series_real
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from level2 import load, to_series  # noqa: E402
from rrlbin.oc import residuals_h0  # noqa: E402
from rrlbin.timing import year_labels  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--series", default="results/validation_v3/series_real")
    a = ap.parse_args()
    d = load(a.series)
    d = d[d.ok.fillna(True).astype(bool)]
    pairs = {"B-R": [], "M-O": []}
    for r in d.to_dict("records"):
        s = to_series(r)
        if not ({0, 1, 2} <= set(s["band"].astype(int))):
            continue
        res = residuals_h0(s) / s["err"]
        yl = year_labels(s["t"])
        df = pd.DataFrame(dict(y=yl, b=s["band"].astype(int), z=res))
        w = df.pivot_table(index="y", columns="b", values="z", aggfunc="mean")
        if {1, 2} <= set(w.columns):
            x = w[[1, 2]].dropna()
            pairs["B-R"] += list(map(tuple, x.to_numpy()))
        if {0, 1} <= set(w.columns):
            x = w[[0, 1]].dropna()
            pairs["M-O"] += list(map(tuple, x.to_numpy()))
    for k, v in pairs.items():
        v = np.array(v)
        if len(v) < 10:
            print(k, "too few pairs")
            continue
        # robust correlation (Spearman) and clipped Pearson
        from scipy.stats import spearmanr
        rho, p = spearmanr(v[:, 0], v[:, 1])
        m = (np.abs(v) < 5).all(axis=1)
        print(f"{k}: {len(v)} same-season pairs; Spearman rho = {rho:.3f} (p = {p:.1e}); Pearson (|z|<5) = "
              f"{np.corrcoef(v[m, 0], v[m, 1])[0, 1]:.3f} (n = {m.sum()})")


if __name__ == "__main__":
    main()
