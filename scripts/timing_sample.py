"""Fail-fast test: per-season timing precision and O-C scatter on a random sample of OGLE RRab.

For each star (OGLE-III [+OGLE-II] + OGLE-IV I band): joint template + per-season delay fit
(rrlbin.timing.fit_timing, K=8), then a weighted quadratic fit to the unwrapped season delays
(tau = c0 + c1 t + c2 t^2, i.e. period offset + constant Pdot). Records the median per-season error,
chi2_nu of the quadratic fit, and the rms of its residuals. chi2_nu >> 1 means timing structure
beyond a constant Pdot (LTTE, Blazhko, stochastic phase jitter, or underestimated errors).

Usage
-----
    python scripts/timing_sample.py --n 1000 --seed 0 --out results/timing_sample/rrab_1000.parquet
"""
from __future__ import annotations

import argparse
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from rrlbin.io import lc_path, read_lc  # noqa: E402
from rrlbin.timing import fit_timing, harmonic_amp_phase, unwrap_delays  # noqa: E402

RAW = Path("data/raw")


def ogle_I(oid: str):
    """Concatenated OGLE-II/III/IV I-band light curve with segment labels (OGLE-II epochs are those at t < 2000)."""
    t3, m3, e3 = read_lc(lc_path(RAW, "ogle3", oid))
    t4, m4, e4 = read_lc(lc_path(RAW, "ogle4", oid))
    seg = np.r_[np.where(t3 < 2000, "O2", "O3"), np.full(t4.size, "O4")]
    return np.r_[t3, t4], np.r_[m3, m4], np.r_[e3, e4], seg


def one(row) -> dict:
    oid, P, T0 = row
    try:
        t, m, e, seg = ogle_I(oid)
        f = fit_timing(t, m, e, seg, P, T0, K=8)
        tau = unwrap_delays(f.tau, P)
        x = (f.t_season - 5000.0) / 1000.0
        W = 1 / f.tau_err ** 2
        X = np.vander(x, 3)
        c, *_ = np.linalg.lstsq(X * np.sqrt(W)[:, None], tau * np.sqrt(W), rcond=None)
        r = tau - X @ c
        A, _ = harmonic_amp_phase(f.coef)
        return dict(ogle_id=oid, ok=True, n_ep=int(f.mask.sum()), n_season=f.season.size,
                    has_O2=bool((f.seg_season == "O2").any()),
                    err_med_s=float(np.median(f.tau_err) * 86400), err_min_s=float(f.tau_err.min() * 86400),
                    chi2nu_lc=f.chi2nu, chi2nu_quad=float(np.sum(W * r ** 2) / max(x.size - 3, 1)),
                    rms_quad_s=float(np.std(r) * 86400), A1=float(A[0]),
                    c2_s=float(c[0] * 86400), baseline=float(np.ptp(f.t_season)))
    except Exception as ex:  # keep going; record the failure
        return dict(ogle_id=oid, ok=False, err=str(ex)[:200])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--out", default="results/timing_sample/rrab_1000.parquet")
    a = ap.parse_args()
    inv = pd.read_parquet("data/lc_inventory.parquet")
    par = pd.read_fwf(RAW / "ogle4_lmc_rrlyr" / "RRab.dat", colspecs=[(0, 20), (22, 28), (37, 46), (58, 68), (70, 75)],
                      names=["ogle_id", "I", "P", "T0", "amp_I"], header=None)
    for c in ["I", "P", "T0", "amp_I"]:   # missing values are '-' in the OGLE tables
        par[c] = pd.to_numeric(par[c], errors="coerce")
    s = inv[(inv.subtype == "RRab") & (inv.n3_I > 0) & (inv.n4_I > 0)].merge(par, on="ogle_id")
    s = s.sample(n=min(a.n, len(s)), random_state=a.seed)
    with Pool(a.workers) as pool:
        rows = pool.map(one, list(zip(s.ogle_id, s.P, s.T0)), chunksize=10)
    out = pd.DataFrame(rows).merge(s[["ogle_id", "I", "P", "amp_I"]], on="ogle_id")
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_parquet(a.out)
    ok = out[out.ok]
    print(f"stars {len(out)}, ok {len(ok)}")
    q = lambda c: np.nanpercentile(ok[c].astype(float), [10, 50, 90]).round(2)
    for c in ["n_season", "err_med_s", "err_min_s", "chi2nu_lc", "chi2nu_quad", "rms_quad_s", "A1", "amp_I", "I"]:
        print(f"{c:12s} p10/50/90 = {q(c)}")
    print("fraction chi2nu_quad > 3:", round((ok.chi2nu_quad > 3).mean(), 3), "; > 10:", round((ok.chi2nu_quad > 10).mean(), 3))
    for lo, hi in [(0, 0.4), (0.4, 0.6), (0.6, 0.8), (0.8, 2)]:
        b = ok[(ok.amp_I >= lo) & (ok.amp_I < hi)]
        print(f"amp_I [{lo},{hi}): N={len(b)}, median err {b.err_med_s.median():.0f} s, median chi2nu_quad {b.chi2nu_quad.median():.2f}")


if __name__ == "__main__":
    main()
