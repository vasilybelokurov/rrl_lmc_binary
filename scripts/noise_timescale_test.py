"""Is the real excess timing noise 'red' (slow wander, absorbed within each survey) or Brownian?
For real MACHO+OGLE stars (common mode + band-lag priors + whole-cycle alignment applied): chi2_nu with no jitter
(a) per band, each band with its own quadratic; (b) joint H0 (one quadratic for all bands + band offsets with priors).
Brownian noise inflates (a) as well; slow wander inflates (b) much more than (a).
Usage: python scripts/noise_timescale_test.py [--n 400] [--series ...]
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from level2 import load, priors_for, to_series  # noqa: E402
from rrlbin.oc import align_bands, apply_common_mode, design, profile_lnl, robust_unwrap  # noqa: E402
from rrlbin.timing import year_labels  # noqa: E402

DAY = 86400.0


def chi2_fit(t, y, e, band, pri):
    X, names = design(t, band, np.average(t, weights=e ** -2))
    _, _, b, _ = profile_lnl(X, y, e, names, np.array([0.0]), pri)
    return float(np.sum(((y - X @ b) / e) ** 2) / max(t.size - X.shape[1], 1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=400)
    ap.add_argument("--series", default="results/real/series_v3")
    ap.add_argument("--cm", default="results/calib/common_mode_v3.json")
    ap.add_argument("--kind", default=None, help="for simulated series: restrict to this class")
    a = ap.parse_args()
    d = load(a.series)
    if "ok" in d:
        d = d[d.ok.fillna(True).astype(bool)]
    if a.kind:
        d = d[d.kind == a.kind]
    d = d[d.band.apply(lambda b: len(set(np.asarray(b).astype(int))) >= 2)]
    d = d.sample(min(a.n, len(d)), random_state=0)
    lag = json.loads(Path("results/calib/band_lag_v3.json").read_text())
    cm = {tuple(map(int, k.split(","))): v / DAY for k, v in json.loads(Path(a.cm).read_text()).items()} if a.cm != "none" else {}
    rows = []
    for r in d.to_dict("records"):
        s = apply_common_mode(to_series(r), cm, year_labels) if cm else to_series(r)
        t, y, e, band = s["t"], s["tau"], s["err"], s["band"].astype(int)
        o = np.argsort(t); t, y, e, band = t[o], y[o], e[o], band[o]
        y = robust_unwrap(t, y, band, s["P"])
        pri = priors_for(s["P"], lag)
        y = align_bands(t, y, e, band, s["P"], pri)
        per = [chi2_fit(t[band == b], y[band == b], e[band == b], np.zeros((band == b).sum(), int), None)
               for b in np.unique(band) if (band == b).sum() >= 5]
        rows.append(dict(per_band=np.median(per), joint=chi2_fit(t, y, e, band, pri)))
    R = pd.DataFrame(rows)
    print(f"N={len(R)}; chi2_nu (no jitter) p50/p75/p90: per band {np.percentile(R.per_band, [50, 75, 90]).round(1)}, "
          f"joint {np.percentile(R.joint, [50, 75, 90]).round(1)}; median joint/per-band ratio {np.median(R.joint / R.per_band):.1f}")


if __name__ == "__main__":
    main()
