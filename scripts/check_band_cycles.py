"""How often is a MACHO band's delay series a whole number of cycles away from the band-lag prior?
For real stars with MACHO: free-offset H0 fit (no priors) -> offset_b; k_b = round((offset_b - prior_mean_b)/P).
Usage: python scripts/check_band_cycles.py
"""
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
d = load("results/real/series_v3")
d = d[d.ok & d.bands.str.contains("M")].sample(1500, random_state=0)
lag = json.loads(Path("results/calib/band_lag_v3.json").read_text())
cm = {tuple(map(int, k.split(","))): v / DAY for k, v in json.loads(Path("results/calib/common_mode_v3.json").read_text()).items()}
st = pd.read_parquet("results/real/stats_v3.parquet", columns=["ogle_id", "s_red_s", "err_med_s", "D"]).set_index("ogle_id")
rows = []
for r in d.to_dict("records"):
    s = apply_common_mode(to_series(r), cm, year_labels)
    band = s["band"].astype(int)
    tau = robust_unwrap(s["t"], s["tau"], band, s["P"])
    X, names = design(s["t"], band, np.average(s["t"], weights=s["err"] ** -2))
    _, _, b, _ = profile_lnl(X, tau, s["err"], names, default_s_grid(s["err"]))
    row = dict(ogle_id=r["ogle_id"], P=s["P"])
    for bb in (1, 2):
        k = f"off{bb}"
        if k in names:
            mu = (lag[str(bb)]["slope"] * s["P"] + lag[str(bb)]["intercept"]) / DAY
            row[f"k{bb}"] = int(np.round((b[names.index(k)] - mu) / s["P"]))
    rows.append(row)
R = pd.DataFrame(rows).set_index("ogle_id").join(st)
for bb in (1, 2):
    c = f"k{bb}"
    print(f"band {bb}: fraction with whole-cycle mismatch |k|>=1: {np.mean(R[c].abs() >= 1):.3f}; values {R[c].value_counts().head(5).to_dict()}")
mis = (R.k1.abs() >= 1) | (R.k2.abs() >= 1)
q = lambda x: np.nanpercentile(x, [50, 90]).round(2)
print("s_red/err p50/p90: mismatched", q(R[mis].s_red_s / R[mis].err_med_s), " aligned", q(R[~mis].s_red_s / R[~mis].err_med_s))
print("D p50/p90: mismatched", q(R[mis].D), " aligned", q(R[~mis].D))
