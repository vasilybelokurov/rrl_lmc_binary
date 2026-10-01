"""Is the simulated photometric noise realistic at the timing level? For stars with null simulations (Level-1 series),
fit the REAL light curves of the same stars and compare, band by band: (1) median per-season delay error (photometric
noise -> timing precision); (2) the number of seasons; (3) the scatter of season delays about each star's own H0 fit,
normalized by the errors (chi2_nu with no extra jitter) - real stars include intrinsic timing noise, nulls do not.

Usage
-----
    python scripts/compare_noise_real_sim.py --sims results/validation_v3/series_ltte
"""
from __future__ import annotations

import argparse
import sys

sys.path.insert(0, __import__("os").path.dirname(__file__))
import chunked  # noqa: E402
from multiprocessing import Pool  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from level1_real import one, sample_table  # noqa: E402
from rrlbin.oc import design, profile_lnl, robust_unwrap  # noqa: E402

DAY = 86400.0


def chi2_h0(r, b):
    t, tau, err, band = (np.asarray(r[k], float) for k in ("t", "tau", "err", "band"))
    m = band == b
    if m.sum() < 5:
        return np.nan
    t, tau, err = t[m], tau[m], err[m]
    o = np.argsort(t)
    t, tau, err = t[o], tau[o], err[o]
    tau = robust_unwrap(t, tau, np.zeros(t.size, int), float(r["P"]))
    X, names = design(t, np.zeros(t.size, int), np.average(t, weights=err ** -2))
    _, _, beta, _ = profile_lnl(X, tau, err, names, np.array([0.0]))
    return float(np.sum(((tau - X @ beta) / err) ** 2) / (t.size - 3))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sims", default="results/validation_v3/series_ltte")
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--kind", default="null", help="simulated class to compare with the real stars")
    a = ap.parse_args()
    sims = chunked.merge_parts(a.sims)
    null = sims[sims.kind == a.kind]
    s = sample_table()
    s = s[s.ogle_id.isin(set(null.ogle_id))]
    jobs = list(zip(s.ogle_id, s.macho_id, s.P.astype(float), s.T0.astype(float), s.ra, s.dec))
    with Pool(a.workers) as p:
        real = pd.DataFrame(p.map(one, jobs, chunksize=4))
    real = real[real.ok].set_index("ogle_id")
    rows = []
    for oid, g in null.groupby("ogle_id"):
        if oid not in real.index:
            continue
        r = real.loc[oid]
        for b in (0, 1, 2):
            rb = np.asarray(r["band"]) == b
            if rb.sum() < 5:
                continue
            sb = [x for x in g.to_dict("records") if (np.asarray(x["band"]) == b).sum() >= 5]
            if not sb:
                continue
            rows.append(dict(ogle_id=oid, band=b,
                             err_real=np.median(np.asarray(r["err"])[rb]) * DAY,
                             err_sim=np.median([np.median(np.asarray(x["err"])[np.asarray(x["band"]) == b]) for x in sb]) * DAY,
                             n_real=int(rb.sum()), n_sim=np.median([(np.asarray(x["band"]) == b).sum() for x in sb]),
                             chi2_real=chi2_h0(r, b), chi2_sim=np.median([chi2_h0(x, b) for x in sb])))
    R = pd.DataFrame(rows)
    R.to_csv("results/validation_v3/noise_real_vs_sim.csv", index=False)
    for b, gb in R.groupby("band"):
        ratio = gb.err_sim / gb.err_real
        print(f"band {b}: {len(gb)} stars; median season error real {gb.err_real.median():.0f} s vs null-sim {gb.err_sim.median():.0f} s; "
              f"ratio sim/real p10/50/90 = {np.percentile(ratio, [10, 50, 90]).round(2)}; seasons real/sim {gb.n_real.median():.0f}/{gb.n_sim.median():.0f}; "
              f"chi2_nu about H0 (no jitter) real p50/p90 {np.nanpercentile(gb.chi2_real, [50, 90]).round(1)} vs sim {np.nanpercentile(gb.chi2_sim, [50, 90]).round(1)}")


if __name__ == "__main__":
    main()
