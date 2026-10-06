"""Test the frozen O-C predictions (tag predictions-2026-10-01) with the post-2016 OGLE-IV seasons (OGLE-team files, 2026-10-06).

Nothing frozen is refitted. Per candidate:
  1. Frozen frame: the v3 season series (public data) with the v3 common mode, robust unwrap and band alignment, exactly as in
     scripts/freeze_predictions.py -> OGLE I delays y_I (the frame of the predictions).
  2. Public OGLE I light curve refitted (fit_timing; deterministic) -> template + O4 zero point; check that its season delays
     equal the frozen ones up to whole cycles (|residual| < 1 s), which fixes the cycle offset k.
  3. New seasons (year labels after the last public OGLE season) from the extended light curve, measured against that FIXED
     template and zero point (rrlbin.timing.delays_fixed_template; same gauge, tests/test_timing.py) -> + k P, then
     unwrapped sequentially from the last public season (largest jump reported). Common mode for the new years: 0 (as frozen).
  4. Score against the frozen predictive distributions (cov_<date>/cov_<id>.npz, interpolated to the season epochs):
     covariance + diag(err^2 + white^2) per hypothesis; ln BF = ln L(H1) - ln L(H0); chi^2 and its p-value per hypothesis;
     for all new seasons and for 2017-2020 / 2022-2026 separately. Also the coverage of the white-noise bootstrap band (95%).
  5. Calibration (Monte Carlo, N_MC draws at the observed epochs and errors): P(ln BF >= observed | H0) and
     P(ln BF <= observed | H1), each hypothesis' own predictive distribution.
  6. POST-HOC H1 variants (sensitivity only; the frozen H1 is the test): 'boot' = frozen H1 covariance + the variance of the
     white-noise bootstrap envelope (orbit-parameter uncertainty); 'h0noise' = best orbit + the star's H0 red-noise GP
     (kriging re-run on the frozen public series with the frozen orbit and gp_h0) - a deliberately generous H1.
  Outliers (the 2026 files are raw database extractions, Soszynski, pers. comm.): per-season 4-sigma clipping (robust scale);
     a season is flagged UNSTABLE if its delay with only gross outliers removed (8 sigma) differs by > 3 sigma; scores are also
     given for the stable seasons only (suffix _stable).
  7. Check of the measurement: new-season delays re-measured with a template refitted to the whole 2010-2026 light curve
     (gauge tied to the frozen frame by the median offset over the public seasons) -> max |difference| [s].

Outputs: results/predictions/test_<today>.csv (per candidate), results/predictions/test_<today>_seasons.csv (per season),
plots/predictions_test/<id>.png.
Usage: python scripts/test_frozen_predictions.py [--date 2026-10-01] [--workers 4] [--cm-shift 0]
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import sys

sys.path.insert(0, __import__("os").path.dirname(__file__))
import chunked  # noqa: E402,F401  (single-threaded BLAS)
from multiprocessing import Pool  # noqa: E402
from pathlib import Path  # noqa: E402

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import chi2 as chi2_dist  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from level2 import load, priors_for, to_series  # noqa: E402
from rrlbin.oc import align_bands, apply_common_mode, default_s_grid, design, profile_lnl, robust_unwrap  # noqa: E402
from rrlbin.pipeline import load_star  # noqa: E402
from rrlbin.predict import gauss_score, interp_prediction, krige  # noqa: E402
from rrlbin.ltte import ltte_delay  # noqa: E402
from rrlbin.timing import delays_fixed_template, fit_timing, year_labels  # noqa: E402

DAY = 86400.0
N_MC = 4000
SPLIT = 9500.0          # HJD' between the 2017-2020 and 2022-2026 blocks (COVID gap)
G = {}


def frozen_frame(r):
    """Series of the frozen predictions (as in freeze_predictions.one): t, y, e, band [d], sorted; plus MACHO offsets [d]."""
    s = apply_common_mode(to_series(r), G["cm"], year_labels)
    t, y, e, band = s["t"], s["tau"], s["err"], s["band"].astype(int)
    o = np.argsort(t)
    t, y, e, band = t[o], y[o], e[o], band[o]
    pri = priors_for(s["P"], G["lag"])
    y = align_bands(t, robust_unwrap(t, y, band, s["P"]), e, band, s["P"], pri)
    X, names = design(t, band, np.average(t, weights=e ** -2))
    _, _, b0, _ = profile_lnl(X, y, e, names, default_s_grid(e), pri)
    offs = {int(k[3:]): float(b0[names.index(k)]) for k in names if k.startswith("off")}
    return t, y, e, band, offs


def score_block(d, err, m1, C1, w1, m0, C0, w0):
    n = d.size
    if n == 0:
        return dict(n=0)
    l1, c1, _ = gauss_score(d, m1, C1 + np.diag(err ** 2 + w1 ** 2))
    l0, c0, _ = gauss_score(d, m0, C0 + np.diag(err ** 2 + w0 ** 2))
    return dict(n=n, lnBF=l1 - l0, chi2_h1=c1, p_h1=float(chi2_dist.sf(c1, n)), chi2_h0=c0, p_h0=float(chi2_dist.sf(c0, n)))


def calib(err, m1, C1, w1, m0, C0, w0, lnbf_obs, rng):
    """Monte Carlo: P(lnBF >= obs | H0), P(lnBF <= obs | H1) at the observed epochs/errors."""
    if err.size == 0:
        return np.nan, np.nan
    S1, S0 = C1 + np.diag(err ** 2 + w1 ** 2), C0 + np.diag(err ** 2 + w0 ** 2)
    L1, L0 = np.linalg.cholesky(S1), np.linalg.cholesky(S0)
    ld1, ld0 = np.log(np.diag(L1)).sum(), np.log(np.diag(L0)).sum()

    def lnbf(D):
        z1, z0 = np.linalg.solve(L1, (D - m1).T), np.linalg.solve(L0, (D - m0).T)
        return -0.5 * (z1 ** 2).sum(0) - ld1 + 0.5 * (z0 ** 2).sum(0) + ld0
    d0 = rng.multivariate_normal(m0, S0, N_MC)
    d1 = rng.multivariate_normal(m1, S1, N_MC)
    return float(np.mean(lnbf(d0) >= lnbf_obs)), float(np.mean(lnbf(d1) <= lnbf_obs))


def one(oid):
    try:
        return _one(oid)
    except Exception as ex:  # recorded, not fatal
        return dict(ogle_id=oid, ok=False, msg=str(ex)[:200]), None


def _one(oid):
    r = G["ser"].loc[oid]
    P, T0 = float(r["P"]), float(r["T0"])
    t, y, e, band, offs = frozen_frame(dict(r.to_dict(), ogle_id=oid))
    I = band == 0
    tI, yI = t[I], y[I]
    # 2. public refit -> template, zero point, cycle offset
    lc = load_star(oid)["I"]
    f = fit_timing(*lc, P, T0, K=8, labels=year_labels(lc[0]))
    lab_pub = year_labels(f.t_season)
    cm_I = np.array([G["cm"].get((0, int(L)), 0.0) for L in lab_pub])
    # match seasons by label (the frozen series has the same seasons)
    lab_fr = year_labels(tI)
    common, i_f, i_y = np.intersect1d(lab_pub, lab_fr, return_indices=True)
    off = yI[i_y] - (f.tau[i_f] - cm_I[i_f])
    k = np.round(off / P)
    resid = np.abs(off - k * P) * DAY
    last = np.argmax(f.t_season[i_f])
    k_last, y_last = k[last], yI[i_y][last]
    # 3. new seasons, fixed template and O4 zero point
    tx, mx, ex, sx = load_star(oid, ogle4="extended")["I"]
    labx = year_labels(tx)
    new = (sx == "O4") & (labx > lab_pub.max())
    zp = f.zp[list(f.seg_names).index("O4")]
    ts, tau, terr, al, al_err, nuse, chi2nu, labs = delays_fixed_template(tx[new], mx[new], ex[new], labx[new], f.coef, P, T0,
                                                                          min_season=15, clip=G["clip"], zp=zp, full=True)
    o = np.argsort(ts)
    ts, tau, terr, al, al_err, nuse, chi2nu, labs = (a[o] for a in (ts, tau, terr, al, al_err, nuse, chi2nu, labs))
    # stability flag: the same seasons with only gross outliers removed (8 sigma); in sparse seasons a 4-sigma clip can remove
    # rising-branch points (most timing information) and move the delay to another minimum (found for 17610 season 23)
    _, tau8, terr8, *_, labs8 = delays_fixed_template(tx[new], mx[new], ex[new], labx[new], f.coef, P, T0, min_season=15, clip=8.0,
                                                      zp=zp, full=True)
    d8 = dict(zip(labs8, zip(tau8, terr8)))
    dtau8 = np.array([(tau[j] - d8[L][0] + P / 2) % P - P / 2 if L in d8 else np.nan for j, L in enumerate(labs)])
    unstable = np.abs(dtau8) > 3 * np.array([max(terr[j], d8[L][1]) if L in d8 else np.inf for j, L in enumerate(labs)])
    yn = tau + k_last * P + G["cm_shift"] / DAY
    prev, jumps = y_last, []
    for j in range(yn.size):
        kk = np.round((prev - yn[j]) / P)
        yn[j] += kk * P
        jumps.append(abs(yn[j] - prev) / P)
        prev = yn[j]
    # 4. score against the frozen predictive distributions
    z = np.load(G["covdir"] / f"cov_{oid}.npz")
    tg = z["t_pred"]
    inside = (ts >= tg[0]) & (ts <= tg[-1])
    d, err_s, tt = yn[inside] * DAY, terr[inside] * DAY, ts[inside]
    m1, C1 = interp_prediction(tg, z["h1_mean"], z["h1_cov"], tt)
    m0, C0 = interp_prediction(tg, z["h0_mean"], z["h0_cov"], tt)
    w1, w0 = float(z["white_h1_s"]), float(z["white_h0_s"])
    out = dict(ogle_id=oid, ok=True, tier=G["tier"][oid], P_orb=G["meta"][oid]["P_orb_d"], A_s=G["meta"][oid]["A_s"],
               n_new=int(inside.sum()), t_first=float(tt.min()) if tt.size else np.nan, t_last=float(tt.max()) if tt.size else np.nan,
               med_err_s=float(np.median(err_s)) if tt.size else np.nan, frame_resid_max_s=float(resid.max()),
               k_const=bool(np.all(k == k_last)), max_jump_cycles=float(max(jumps)) if jumps else np.nan,
               alpha_dev_max=float(np.max(np.abs(al[inside] - 1))) if tt.size else np.nan)
    stab = ~unstable[inside]
    out["n_unstable"] = int((~stab).sum())
    for lab_, sel in (("all", np.ones(tt.size, bool)), ("1720", tt < SPLIT), ("2226", tt >= SPLIT), ("stable", stab)):
        ix = np.flatnonzero(sel)
        sc = score_block(d[ix], err_s[ix], m1[ix], C1[np.ix_(ix, ix)], w1, m0[ix], C0[np.ix_(ix, ix)], w0)
        out.update({f"{k_}_{lab_}": v for k_, v in sc.items()})
    rng = np.random.default_rng(int(oid[-5:]))
    out["pfa_h0"], out["pmiss_h1"] = calib(err_s, m1, C1, w1, m0, C0, w0, out.get("lnBF_all", np.nan), rng)
    # 6. post-hoc H1 variants
    pr = G["pred"][G["pred"].ogle_id == oid].sort_values("t_hjdp")
    vb = ((np.interp(tt, pr.t_hjdp, pr.hi95) - np.interp(tt, pr.t_hjdp, pr.lo95)) / 3.92) ** 2
    sc = score_block(d, err_s, m1, C1 + np.diag(vb), w1, m0, C0, w0)
    out.update(lnBF_boot=sc.get("lnBF", np.nan), p_h1_boot=sc.get("p_h1", np.nan))
    km = G["meta"][oid]
    orb = lambda x: ltte_delay(np.asarray(x), km["P_orb_d"], km["A_s"] / DAY, km["e"], km["omega"], km["t_peri"])
    g0 = km["gp_h0"]
    pri = priors_for(P, G["lag"])
    mh, Ch = krige(t, y - orb(t), e, band, tg, g0["A"] / DAY, g0["ell"], g0["s"] / DAY, pri)
    mh, Ch = (mh + orb(tg)) * DAY, Ch * DAY ** 2
    m1n, C1n = interp_prediction(tg, mh, Ch, tt)
    sc = score_block(d, err_s, m1n, C1n, g0["s"], m0, C0, w0)
    out.update(lnBF_h0noise=sc.get("lnBF", np.nan), p_h1_h0noise=sc.get("p_h1", np.nan))
    # 7. re-measure with a template refitted to the full extended light curve
    fx = fit_timing(tx, mx, ex, sx, P, T0, K=8, labels=labx)
    labfx = year_labels(fx.t_season)
    cmx = np.array([G["cm"].get((0, int(L)), 0.0) for L in labfx])
    cpub, a_, b_ = np.intersect1d(labfx, lab_fr, return_indices=True)
    dlt = yI[b_] - (fx.tau[a_] - cmx[a_])
    dlt = dlt - P * np.round((dlt - dlt[-1]) / P)
    shift = np.median(dlt)
    cn, a2, b2 = np.intersect1d(labfx, labs[inside], return_indices=True)
    if cn.size:
        yx = fx.tau[a2] + shift + G["cm_shift"] / DAY
        yx = yx + P * np.round((yn[inside][b2] - yx) / P)
        out["remeasure_maxdiff_s"] = float(np.max(np.abs(yx - yn[inside][b2])) * DAY)
        out["remeasure_public_rms_s"] = float(np.std(dlt - shift) * DAY)
    # bootstrap (white-noise) 95% band coverage, for reference
    lo, hi = np.interp(tt, pr.t_hjdp, pr.lo95), np.interp(tt, pr.t_hjdp, pr.hi95)
    out["boot95_cover"] = float(np.mean((d >= lo) & (d <= hi))) if tt.size else np.nan
    seasons = pd.DataFrame(dict(ogle_id=oid, label=labs[inside], t=tt, y_s=d, err_s=err_s, alpha=al[inside], alpha_err=al_err[inside],
                                n_epochs=nuse[inside], chi2nu=chi2nu[inside], h1_mean=m1, h1_sd=np.sqrt(np.diag(C1) + w1 ** 2),
                                h0_mean=m0, h0_sd=np.sqrt(np.diag(C0) + w0 ** 2), dtau_clip8_s=dtau8[inside] * DAY,
                                unstable=unstable[inside]))
    seasons["z_h1"] = (seasons.y_s - seasons.h1_mean) / np.hypot(seasons.h1_sd, seasons.err_s)
    seasons["z_h0"] = (seasons.y_s - seasons.h0_mean) / np.hypot(seasons.h0_sd, seasons.err_s)
    plot(oid, out, t, y, e, band, offs, z, seasons)
    return out, seasons


def plot(oid, out, t, y, e, band, offs, z, S):
    fig, ax = plt.subplots(2, 1, figsize=(10, 6.5), sharex=True, gridspec_kw=dict(height_ratios=[3, 1.3]))
    a = ax[0]
    for b, col, lab in [(0, "C0", "OGLE I (public)"), (1, "C1", "MACHO B"), (2, "C2", "MACHO R")]:
        m = band == b
        if m.any():
            a.errorbar(t[m], (y[m] - offs.get(b, 0.0)) * DAY, e[m] * DAY, fmt="o", ms=3, color=col, label=lab)
    tg = z["t_pred"]
    s1 = np.sqrt(np.diag(z["h1_cov"]) + float(z["white_h1_s"]) ** 2)
    s0 = np.sqrt(np.diag(z["h0_cov"]) + float(z["white_h0_s"]) ** 2)
    a.fill_between(tg, z["h1_mean"] - 2 * s1, z["h1_mean"] + 2 * s1, color="0.6", alpha=0.5, label="frozen H1 (orbit), ±2σ")
    a.plot(tg, z["h1_mean"], "k-", lw=1)
    a.fill_between(tg, z["h0_mean"] - 2 * s0, z["h0_mean"] + 2 * s0, color="C3", alpha=0.12, label="frozen H0 (no orbit), ±2σ")
    a.plot(tg, z["h0_mean"], "C3--", lw=1)
    a.errorbar(S.t, S.y_s, S.err_s, fmt="s", ms=5, mfc="gold", mec="k", ecolor="k", label="OGLE I 2017–2026 (new)")
    a.axvline(7516, c="k", lw=0.6, ls=":")
    a.set(ylabel="O−C [s] (OGLE I gauge)",
          title=f"{oid} (Tier {out['tier']}): P_orb={out['P_orb']:.0f} d; new seasons {out['n_new']}; "
                f"ln BF(H1:H0) = {out.get('lnBF_all', np.nan):.1f}; p(χ²|H1) = {out.get('p_h1_all', np.nan):.2g}")
    lo = min(np.min((y - np.vectorize(lambda b: offs.get(b, 0.0))(band)) * DAY), S.y_s.min() if len(S) else 0)
    hi = max(np.max((y - np.vectorize(lambda b: offs.get(b, 0.0))(band)) * DAY), S.y_s.max() if len(S) else 0)
    a.set_ylim(lo - 0.25 * (hi - lo), hi + 0.25 * (hi - lo))
    a.legend(fontsize=7, ncol=3)
    b_ = ax[1]
    b_.errorbar(S.t, S.z_h1, 1 / np.hypot(S.h1_sd, S.err_s) * S.err_s, fmt="ko", ms=4, label="vs H1")
    b_.plot(S.t, S.z_h0.clip(-20, 20), "C3x", label="vs H0 (clipped ±20)")
    b_.axhline(0, c="k", lw=0.6)
    b_.axhspan(-2, 2, color="0.85")
    b_.set(xlabel="HJD − 2450000", ylabel="(new − pred)/σ")
    b_.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(G["pl"] / f"{oid}.png", dpi=100)
    plt.close(fig)


def _init(state):
    G.update(state)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--date", default="2026-10-01")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--cm-shift", type=float, default=0.0, help="common-mode shift [s] added to all new seasons (sensitivity)")
    ap.add_argument("--clip", type=float, default=4.0, help="per-season sigma clipping of the new epochs (robust scale)")
    ap.add_argument("--tag", default="")
    a = ap.parse_args()
    today = dt.date.today().isoformat()
    out_dir, pl = Path("results/predictions"), Path("plots/predictions_test" + a.tag)
    pl.mkdir(parents=True, exist_ok=True)
    meta = json.loads(Path(f"results/predictions/predictions_{a.date}_meta.json").read_text())["candidates"]
    covdir = Path(f"results/predictions/cov_{a.date}")
    ids = sorted(f.stem[4:] for f in covdir.glob("cov_*.npz"))
    ser = load("results/real/series_v3").set_index("ogle_id")
    cm = {tuple(map(int, k.split(","))): v / DAY for k, v in json.loads(Path("results/calib/common_mode_v3.json").read_text()).items()}
    lag = json.loads(Path("results/calib/band_lag_v3.json").read_text())
    pred = pd.read_parquet(f"results/predictions/predictions_{a.date}.parquet")
    have = [o for o in ids if Path(f"data/raw/ogle4_lmc_rrlyr_2026/phot/I/{o}.dat").exists()]
    print(f"frozen candidates {len(ids)}; with new OGLE data {len(have)}; missing {sorted(set(ids) - set(have))}", flush=True)
    state = dict(meta=meta, covdir=covdir, ser=ser.loc[have], cm=cm, lag=lag, pred=pred[pred.ogle_id.isin(have)], pl=pl,
                 tier={o: meta[o]["tier"] for o in have}, cm_shift=a.cm_shift, clip=a.clip)
    with Pool(a.workers, initializer=_init, initargs=(state,)) as pool:
        res = pool.map(one, have, chunksize=1)
    R = pd.DataFrame([r[0] for r in res])
    S = pd.concat([r[1] for r in res if r[1] is not None], ignore_index=True)
    R.to_csv(out_dir / f"test_{today}{a.tag}.csv", index=False)
    S.to_csv(out_dir / f"test_{today}{a.tag}_seasons.csv", index=False)
    cols = ["ogle_id", "tier", "P_orb", "A_s", "n_new", "med_err_s", "frame_resid_max_s", "k_const", "max_jump_cycles",
            "alpha_dev_max", "n_unstable", "lnBF_1720", "lnBF_2226", "lnBF_all", "p_h1_all", "p_h0_all", "lnBF_stable", "p_h1_stable", "pfa_h0", "pmiss_h1", "lnBF_boot", "p_h1_boot",
            "lnBF_h0noise", "p_h1_h0noise", "remeasure_maxdiff_s", "remeasure_public_rms_s", "boot95_cover"]
    print(R.sort_values(["tier", "lnBF_all"], ascending=[True, False])[[c for c in cols if c in R]].round(3).to_string(index=False))
    if "msg" in R:
        print("failures:", R[~R.ok.astype(bool)][["ogle_id", "msg"]].to_string(index=False))


if __name__ == "__main__":
    main()
