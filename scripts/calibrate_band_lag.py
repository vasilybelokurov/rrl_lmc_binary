"""Band-lag relations for the MACHO B and R delay offsets relative to OGLE I, from stars observed by MACHO and OGLE-II
in the same years (1997-2000), where the offset is fixed by the data.

Per star: H0 fit (quadratic + free band offsets; white jitter profiled; common mode applied if given) -> offset_b.
Relation: offset_b = slope * P + intercept (robust: iteratively 4-sigma clipped weighted fit), scatter = 1.4826 MAD of
residuals, corrected in quadrature for the median measurement error. Output JSON {band: {slope, intercept, sd, n}}.

Usage
-----
    python scripts/calibrate_band_lag.py --series results/real/series_v3 --apply-cm results/calib/common_mode_v3.json \
        --out results/calib/band_lag_v3.json
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from level2 import load, to_series  # noqa: E402
from rrlbin.oc import apply_common_mode, default_s_grid, design, profile_lnl, robust_unwrap  # noqa: E402
from rrlbin.timing import year_labels  # noqa: E402

DAY = 86400.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--series", required=True)
    ap.add_argument("--apply-cm", default=None)
    ap.add_argument("--out", default="results/calib/band_lag_v3.json")
    a = ap.parse_args()
    d = load(a.series)
    d = d[d.ok.fillna(True).astype(bool)] if "ok" in d else d
    cm = {}
    if a.apply_cm:
        cm = {tuple(map(int, k.split(","))): v / DAY for k, v in json.loads(Path(a.apply_cm).read_text()).items()}
    rows = []
    for r in d.to_dict("records"):
        s = apply_common_mode(to_series(r), cm, year_labels) if cm else to_series(r)
        band = s["band"].astype(int)
        if not ((band == 0) & (s["t"] < 2000)).any():
            continue
        tau = robust_unwrap(s["t"], s["tau"], band, s["P"])
        X, names = design(s["t"], band, np.average(s["t"], weights=s["err"] ** -2))
        _, _, b, cov = profile_lnl(X, tau, s["err"], names, default_s_grid(s["err"]))
        for bb in (1, 2):
            k = f"off{bb}"
            # require overlap: MACHO band seasons in 1997-2000 while OGLE-II also observed
            if k in names and ((band == bb) & (s["t"] > 400) & (s["t"] < 1700)).sum() >= 1:
                off = b[names.index(k)]
                off -= s["P"] * np.round(off / s["P"])
                rows.append(dict(ogle_id=r["ogle_id"], band=bb, P=s["P"], off_s=off * DAY,
                                 err_s=np.sqrt(cov[names.index(k), names.index(k)]) * DAY))
    R = pd.DataFrame(rows)
    res = {}
    for bb, g in R.groupby("band"):
        g = g[g.err_s < 300]
        m = np.ones(len(g), bool)
        for _ in range(5):
            c = np.polyfit(g.P[m], g.off_s[m], 1, w=1 / np.hypot(g.err_s[m], 300))
            r = g.off_s - np.polyval(c, g.P)
            mad = 1.4826 * np.median(np.abs(r[m] - np.median(r[m])))
            m = np.abs(r) < 4 * mad
        sd = float(np.sqrt(max(mad ** 2 - np.median(g.err_s[m]) ** 2, 0)))
        res[str(bb)] = dict(slope=float(c[0]), intercept=float(c[1]), sd=max(sd, 50.0), n=int(m.sum()),
                            mad_s=float(mad), median_err_s=float(np.median(g.err_s[m])),
                            lag_cycles_median=float(np.median(g.off_s[m] / DAY / g.P[m])))
        print(f"band {bb}: N={m.sum()}, offset = {c[0]:.0f} P + {c[1]:.0f} s; residual MAD {mad:.0f} s, "
              f"median error {np.median(g.err_s[m]):.0f} s -> intrinsic sd {sd:.0f} s; median lag "
              f"{res[str(bb)]['lag_cycles_median']:.4f} cycles")
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
