"""Freeze dated O-C predictions for 2016.4-2026.5 (OGLE I band) for the Tier-1/2 candidates, from PUBLIC data only
(OGLE-II/III/IV to 2016.3 + MACHO 1992-1999), ahead of a request for post-2016 OGLE photometry.

For each candidate: Keplerian model (orbit + quadratic + MACHO offsets with band-lag priors + white jitter), residual-bootstrap
envelope of the predicted delay curve (200 refits); and the no-orbit model (quadratic + offsets, H0) prediction with its
parameter uncertainty. Delays are in the fundamental-phase gauge of the OGLE I template with the linear ephemeris (P, T0) of
the OGLE-IV catalogue and the v3 common-mode correction (OGLE years beyond 2016 have no correction; assumed 0).

Outputs (frozen; do not edit):
  results/predictions/predictions_<date>.parquet   long table: ogle_id, t_hjdp, model, best, lo95, lo68, med, hi68, hi95
  results/predictions/predictions_<date>_meta.json  per-candidate parameters, ephemeris, inputs, sha256 of the table, git commit
  plots/predictions/<id>.png                        data + predictions
Usage: python scripts/freeze_predictions.py --workers 2
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import subprocess
import sys

sys.path.insert(0, __import__("os").path.dirname(__file__))
import chunked  # noqa: E402,F401
from multiprocessing import Pool  # noqa: E402
from pathlib import Path  # noqa: E402

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from astropy.time import Time  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from level2 import load, priors_for, to_series  # noqa: E402
from rrlbin.kepler_fit import bootstrap_predictions  # noqa: E402
from rrlbin.ltte import ltte_delay  # noqa: E402
from rrlbin.oc import gp_null  # noqa: E402
from rrlbin.predict import predict_h0_h1  # noqa: E402
from rrlbin.oc import align_bands, apply_common_mode, default_s_grid, design, profile_lnl, robust_unwrap  # noqa: E402
from rrlbin.timing import year_labels  # noqa: E402

DAY = 86400.0
G = {}


def one(args):
    r, P0, jit, P_puls, T0, gp0 = args
    s = apply_common_mode(to_series(r), G["cm"], year_labels)
    t, y, e, band = s["t"], s["tau"], s["err"], s["band"].astype(int)
    o = np.argsort(t)
    t, y, e, band = t[o], y[o], e[o], band[o]
    pri = priors_for(s["P"], G["lag"])
    y = align_bands(t, robust_unwrap(t, y, band, s["P"]), e, band, s["P"], pri)
    tp = G["tp"]
    best, cur = bootstrap_predictions(t, y, e, band, P0, tp, s_jit=jit / DAY, priors=pri, n_boot=200, seed=1)
    # no-orbit prediction (+ parameter-uncertainty band)
    sig = np.sqrt(e ** 2 + (jit / DAY) ** 2)
    t_ref = np.average(t, weights=sig ** -2)
    X, names = design(t, band, t_ref)
    _, _, b0, cov0 = profile_lnl(X, y, e, names, default_s_grid(e), pri)
    Xp, _ = design(tp, np.zeros(tp.size, int), t_ref)
    i3 = [names.index(k) for k in ("q2", "q1", "q0")]
    h0 = Xp[:, :3] @ b0[i3]
    h0_sd = np.sqrt(np.einsum("ij,jk,ik->i", Xp[:, :3], cov0[np.ix_(i3, i3)], Xp[:, :3]))
    q = np.percentile(cur, [2.5, 16, 50, 84, 97.5], axis=0) if len(cur) else np.full((5, tp.size), np.nan)
    # red-noise (GP) predictive distributions: H0 with the star's REML noise fit; H1 with the noise refitted after removing
    # the best orbit (refit the Keplerian once to get its parameters)
    from rrlbin.kepler_fit import fit_keplerian
    kf = fit_keplerian(t, y, e, band, P0, s_jit=jit / DAY, priors=pri, n_boot=0)
    kep = dict(P=kf["P"], A_s=kf["A_s"], e=kf["e"], omega=kf["omega"], t_p=kf["t_p"])
    orb_t = ltte_delay(t, kep["P"], kep["A_s"] / DAY, kep["e"], kep["omega"], kep["t_p"])
    s1, A1, l1 = gp_null(t, y - orb_t, e, band, pri)
    (mh0, ch0), (mh1, ch1) = predict_h0_h1(t, y, e, band, tp, kep, (gp0[0] / DAY, gp0[1] / DAY, gp0[2]), (s1, A1, l1), pri)
    rows = pd.DataFrame(dict(ogle_id=r["ogle_id"], t_hjdp=tp, best=best * DAY, lo95=q[0] * DAY, lo68=q[1] * DAY, med=q[2] * DAY,
                             hi68=q[3] * DAY, hi95=q[4] * DAY, noorbit=h0 * DAY, noorbit_sd=h0_sd * DAY,
                             h1_mean=mh1 * DAY, h1_sd=np.sqrt(np.diag(ch1)) * DAY, h0_mean=mh0 * DAY, h0_sd=np.sqrt(np.diag(ch0)) * DAY))
    np.savez(G["outdir"] / f"cov_{r['ogle_id']}.npz", t_pred=tp, h1_mean=mh1 * DAY, h1_cov=ch1 * DAY ** 2, h0_mean=mh0 * DAY,
             h0_cov=ch0 * DAY ** 2, white_h1_s=s1 * DAY, white_h0_s=gp0[0])
    data = dict(t=t, y=y * DAY, e=e * DAY, band=band, offsets={k: float(b0[names.index(k)] * DAY) for k in names if k.startswith("off")})
    return r["ogle_id"], rows, data, dict(P_puls=P_puls, T0_ogle=T0, n_boot=int(len(cur)), jitter_s=jit,
                                         gp_h0=dict(s=gp0[0], A=gp0[1], ell=gp0[2]), gp_h1=dict(s=s1 * DAY, A=A1 * DAY, ell=l1))


def _init(cm, lag, tp, outdir):
    G.update(cm=cm, lag=lag, tp=tp, outdir=outdir)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--tiers", default="1,2")
    a = ap.parse_args()
    today = dt.date.today().isoformat()
    out, pl = Path("results/predictions"), Path("plots/predictions")
    out.mkdir(parents=True, exist_ok=True)
    pl.mkdir(parents=True, exist_ok=True)
    T = pd.read_csv("results/partC/partC_tiers.csv")
    T = T[T.tier.isin([int(x) for x in a.tiers.split(",")]) & np.isfinite(T.P_kep) & ~T.e_at_bound]
    par = pd.read_fwf("data/raw/ogle4_lmc_rrlyr/RRab.dat", colspecs=[(0, 20), (37, 46), (58, 68)], names=["ogle_id", "P", "T0"],
                      header=None).set_index("ogle_id")
    st = pd.read_parquet("results/real/stats_v3.parquet", columns=["ogle_id", "jit1_s", "s_gp_s", "A_gp_s", "ell_gp_d"]).set_index("ogle_id")
    ser = load("results/real/series_v3").set_index("ogle_id")
    cm = {tuple(map(int, k.split(","))): v / DAY for k, v in json.loads(Path("results/calib/common_mode_v3.json").read_text()).items()}
    lag = json.loads(Path("results/calib/band_lag_v3.json").read_text())
    tp = np.linspace(Time("2016-06-01").jd, Time("2026-07-01").jd, 123) - 2450000.0
    jobs = [(dict(ser.loc[o].to_dict(), ogle_id=o), float(P0), float(st.loc[o, "jit1_s"]), float(par.loc[o, "P"]), float(par.loc[o, "T0"]),
             (float(st.loc[o, "s_gp_s"]), float(st.loc[o, "A_gp_s"]), float(st.loc[o, "ell_gp_d"])))
            for o, P0 in zip(T.ogle_id, T.P_kep)]
    covdir = out / f"cov_{today}"
    covdir.mkdir(exist_ok=True)
    with Pool(a.workers, initializer=_init, initargs=(cm, lag, tp, covdir)) as pool:
        res = pool.map(one, jobs, chunksize=1)
    tab = pd.concat([r[1] for r in res], ignore_index=True)
    fn = out / f"predictions_{today}.parquet"
    tab.to_parquet(fn)
    sha = hashlib.sha256(fn.read_bytes()).hexdigest()
    commit = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    meta = dict(created=dt.datetime.now().isoformat(timespec="seconds"), table=str(fn), sha256=sha, git_commit_of_code=commit,
                inputs=["results/real/series_v3 (public OGLE-II/III/IV to 2016.3, MACHO B/R 1992-1999)",
                        "results/calib/common_mode_v3.json", "results/calib/band_lag_v3.json", "results/partC/partC_tiers.csv"],
                gauge="OGLE I fundamental-phase gauge; linear ephemeris P, T0 of the OGLE-IV catalogue (RRab.dat); delays in s",
                candidates={})
    Tt = T.set_index("ogle_id")
    for oid, rows, data, info in res:
        k = Tt.loc[oid]
        meta["candidates"][oid] = dict(tier=int(k.tier), P_orb_d=float(k.P_kep), A_s=float(k.A_kep_s), e=float(k.e_kep),
                                       omega=float(k.omega_kep), t_peri=float(k.tp_kep), **info)
        fig, ax = plt.subplots(figsize=(10, 4.2))
        for b, col, lab in [(0, "C0", "OGLE I"), (1, "C1", "MACHO B (offset removed)"), (2, "C2", "MACHO R (offset removed)")]:
            m = data["band"] == b
            if m.any():
                yy = data["y"][m] - data["offsets"].get(f"off{b}", 0.0)
                ax.errorbar(data["t"][m], yy, data["e"][m], fmt="o", ms=3, color=col, label=lab)
        ax.fill_between(rows.t_hjdp, rows.h1_mean - 2 * rows.h1_sd, rows.h1_mean + 2 * rows.h1_sd, color="0.6", alpha=0.5,
                        label="H1 orbit + red noise, ±2σ")
        ax.plot(rows.t_hjdp, rows.h1_mean, "k-", lw=1)
        ax.fill_between(rows.t_hjdp, rows.h0_mean - 2 * rows.h0_sd, rows.h0_mean + 2 * rows.h0_sd, color="C3", alpha=0.15,
                        label="H0 no orbit + red noise, ±2σ")
        ax.plot(rows.t_hjdp, rows.h0_mean, "C3--", lw=1)
        ax.axvline(7516, c="k", lw=0.6, ls=":")
        ax.set(xlabel="HJD − 2450000", ylabel="O−C [s] (OGLE I gauge)",
               title=f"{oid} (Tier {int(k.tier)}): frozen prediction {today}, P_orb={k.P_kep:.0f} d, e={k.e_kep:.2f}")
        ax.legend(fontsize=7, ncol=3)
        fig.tight_layout()
        fig.savefig(pl / f"{oid}.png", dpi=100)
        plt.close(fig)
    (out / f"predictions_{today}_meta.json").write_text(json.dumps(meta, indent=1))
    print(f"frozen {len(res)} candidates -> {fn} (sha256 {sha[:16]}...); plots in {pl}")


if __name__ == "__main__":
    main()
