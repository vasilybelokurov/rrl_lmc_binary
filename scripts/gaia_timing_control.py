"""Control for the Gaia DR3 timing test: for ~300 random non-candidate RRab, compare the Gaia half-window delay change with
the change measured by OGLE itself over the same epochs (OGLE season delays, common-mode corrected, linearly interpolated;
extrapolated from the last two seasons when t > the last OGLE season). If Gaia timing is reliable, z = (obs - OGLE)/sigma
has unit scatter. Same Gaia measurement as scripts/gaia_timing.py.

Usage: python scripts/gaia_timing_control.py [--n 300]
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import sqlutilpy

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from level2 import load, to_series  # noqa: E402
from rrlbin.oc import apply_common_mode, robust_unwrap  # noqa: E402
from rrlbin.timing import delays_fixed_template, fit_timing, year_labels  # noqa: E402

DAY = 86400.0
SPLIT = 7470.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=300)
    ap.add_argument("--template", choices=["gaia", "ogle"], default="ogle")
    a = ap.parse_args()
    cand = set(pd.read_csv("results/partC/partC_tiers.csv").ogle_id)
    g = pd.read_parquet("../rrl_lmc_rotation/data/ogle4_lmc_rrl_gaiadr3.parquet", columns=["ogle_id", "source_id", "phot_g_mean_mag"])
    par = pd.read_fwf("data/raw/ogle4_lmc_rrlyr/RRab.dat", colspecs=[(0, 20), (37, 46), (58, 68)], names=["ogle_id", "P", "T0"], header=None)
    st = pd.read_parquet("results/real/stats_v3.parquet", columns=["ogle_id", "jit0_s", "err_med_s"])
    T = g.merge(par, on="ogle_id").merge(st, on="ogle_id")
    # control stars: non-candidates with quiet O-C (jitter < 2 sigma) and similar G to the candidates
    T = T[(T.source_id > 0) & ~T.ogle_id.isin(cand) & (T.jit0_s < 2 * T.err_med_s) & T.phot_g_mean_mag.between(18.8, 19.8)]
    T = T.sample(min(a.n, len(T)), random_state=0)
    ser = load("results/real/series_v3")
    ser = ser[ser.ogle_id.isin(set(T.ogle_id))].set_index("ogle_id")
    cm = {tuple(map(int, k.split(","))): v / DAY for k, v in json.loads(Path("results/calib/common_mode_v3.json").read_text()).items()}
    q = """select e.source_id, e.g_transit_time, e.g_transit_mag, e.g_transit_flux_over_error, e.variability_flag_g_reject
           from mytab t join gaia_dr3.epoch_photometry e on e.source_id = t.sid"""
    res = sqlutilpy.local_join(q, "mytab", (T.source_id.astype(np.int64).to_numpy(),), ("sid",), asDict=True)
    ep = {int(s): tuple(np.asarray(x) for x in v) for s, *v in zip(res["source_id"], res["g_transit_time"], res["g_transit_mag"],
                                                                   res["g_transit_flux_over_error"], res["variability_flag_g_reject"])}
    rows = []
    for r in T.itertuples():
        if int(r.source_id) not in ep or r.ogle_id not in ser.index:
            continue
        tt, mm, snr, rej = ep[int(r.source_id)]
        tt, mm, snr = tt.astype(float), mm.astype(float), snr.astype(float)
        ok = np.isfinite(tt) & np.isfinite(mm) & np.isfinite(snr) & (snr > 0) & ~rej.astype(bool)
        t = tt[ok] + 2455197.5 - 2450000.0
        lab = (t > SPLIT).astype(int)
        if min((lab == 0).sum(), (lab == 1).sum()) < 10:
            continue
        try:
            if a.template == "ogle":
                from rrlbin.pipeline import load_star
                tI, mI, eI, sI = load_star(r.ogle_id)["I"]
                fI = fit_timing(tI, mI, eI, sI, float(r.P), float(r.T0), K=8, labels=year_labels(tI))
                ts, taus, errs = delays_fixed_template(t, mm[ok], 1.0857 / snr[ok], lab, fI.coef, float(r.P), float(r.T0))
            else:
                f = fit_timing(t, mm[ok], 1.0857 / snr[ok], np.full(t.size, "G"), float(r.P), float(r.T0), K=4, labels=lab,
                               min_season=10, clip=5.0)
                ts, taus, errs = f.t_season, f.tau, f.tau_err
        except Exception:
            continue
        if ts.size < 2:
            continue
        class _F: pass
        f = _F(); f.t_season, f.tau, f.tau_err = ts, taus, errs
        s = apply_common_mode(to_series(dict(ser.loc[r.ogle_id].to_dict(), ogle_id=r.ogle_id)), cm, year_labels)
        m0 = s["band"].astype(int) == 0
        to, yo, eo = s["t"][m0], s["tau"][m0], s["err"][m0]
        o = np.argsort(to)
        to, yo, eo = to[o], robust_unwrap(to[o], yo[o], np.zeros(o.size, int), float(r.P)), eo[o]
        def oc_w(x):
            """weights w over OGLE seasons such that OGLE delay at x = w . yo (linear interpolation / extrapolation)."""
            w = np.zeros(to.size)
            j = np.searchsorted(to, x)
            j = min(max(j, 1), to.size - 1)
            f_ = (x - to[j - 1]) / (to[j] - to[j - 1])
            w[j - 1], w[j] = 1 - f_, f_
            return w
        wdiff = oc_w(f.t_season[1]) - oc_w(f.t_season[0])
        d_ogle = wdiff @ yo
        sig_ogle = np.sqrt(np.sum((wdiff * eo) ** 2))
        d_obs = (f.tau[1] - f.tau[0] + r.P / 2) % r.P - r.P / 2
        sig_g = np.hypot(*f.tau_err[:2])
        sig = np.hypot(sig_g, sig_ogle)
        rows.append(dict(ogle_id=r.ogle_id, nG=int(t.size), d_obs_s=d_obs * DAY, d_ogle_s=d_ogle * DAY, sigma_gaia_s=sig_g * DAY,
                         sigma_ogle_s=sig_ogle * DAY, sigma_s=sig * DAY, z=(d_obs - d_ogle) / sig, z_gaia_only=(d_obs - d_ogle) / sig_g))
    R = pd.DataFrame(rows)
    R.to_csv("results/partC/gaia_timing_control.csv", index=False)
    z = R.z.to_numpy()
    print(f"OGLE reference error median {np.median(R.sigma_ogle_s):.0f} s; Gaia difference error median {np.median(R.sigma_gaia_s):.0f} s")
    print(f"control stars: {len(R)}; total difference error median {np.median(R.sigma_s):.0f} s; "
          f"z: median {np.median(z):.2f}, robust sd {1.4826 * np.median(np.abs(z - np.median(z))):.2f}, sd {z.std():.2f}; "
          f"|z| > 3: {np.mean(np.abs(z) > 3):.3f}; |z| > 5: {np.mean(np.abs(z) > 5):.3f}")


if __name__ == "__main__":
    main()
