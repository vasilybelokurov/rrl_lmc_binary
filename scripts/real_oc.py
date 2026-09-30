"""Science run: season delays (O-C) of all OGLE-III+IV RRab (+ MACHO B where fetched), with detection statistics.

Per star: OGLE I (O2/O3/O4) timing fit with year seasons; if a MACHO tile exists, MACHO B timing fit; combined delay
series with a free MACHO offset. Stores the series (t, tau, err, survey flag, alpha, alpha_err) for later modelling, plus:
D (circular LTTE vs quadratic, white jitter), P_best, amplitude, the alpha veto statistic, the H0 red-noise fit (s, q),
and, where possible, the predictive score (held-out MACHO seasons before 1997).

Usage
-----
    python scripts/real_oc.py --workers 6 --out results/real/oc_all.parquet [--limit 200]
"""
from __future__ import annotations

import argparse
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from rrlbin.io import read_macho, read_ogle4_ident  # noqa: E402
from rrlbin.ltte import fit_red_null, oc_search, period_grid, predictive_score  # noqa: E402
from rrlbin.timing import fit_timing, unwrap_delays, year_labels  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from timing_sample import ogle_I  # noqa: E402

DAY = 86400.0


def alpha_chi2(a, ae):
    w = ae ** -2
    return float(np.sum(w * (a - np.sum(w * a) / w.sum()) ** 2) / max(a.size - 1, 1))


def one(row):
    oid, mid, P, T0, ra, dec = row
    try:
        tO, mO, eO, segO = ogle_I(oid)
        fO = fit_timing(tO, mO, eO, segO, P, T0, K=8, labels=year_labels(tO))
        t, tau, err = fO.t_season, unwrap_delays(fO.tau, P), fO.tau_err
        alpha, alpha_e = fO.alpha, fO.alpha_err
        flag = np.zeros(t.size)
        has_M = False
        if isinstance(mid, str) and mid.count(".") == 2:
            f_, t_, s_ = map(int, mid.split("."))
            tile = Path(f"data/raw/macho/{f_}.{t_}.parquet")
            if tile.exists():
                tM, mM, eM = read_macho(tile, s_, ra, dec, band="b")
                if tM.size >= 100:
                    fM = fit_timing(tM, mM, eM, np.full(tM.size, "M"), P, T0, K=8, labels=year_labels(tM))
                    t = np.r_[fM.t_season, t]
                    tau = np.r_[unwrap_delays(fM.tau, P), tau]
                    err = np.r_[fM.tau_err, err]
                    alpha, alpha_e = np.r_[fM.alpha, alpha], np.r_[fM.alpha_err, alpha_e]
                    flag = np.r_[np.ones(fM.season.size), flag]
                    has_M = True
        X = flag[:, None] if has_M else None
        r = oc_search(t, tau, err, period_grid(np.ptp(t)), X_extra=X)
        rn = fit_red_null(t, tau, err, X_extra=X)
        out = dict(ogle_id=oid, ok=True, has_M=has_M, has_O2=bool((segO == "O2").any()), n_season=t.size,
                   baseline=float(np.ptp(t)), D=r["D"], P_best=r["P_best"], amp_s=r["amp"] * DAY, jit0_s=r["jit0"] * DAY,
                   err_med_s=float(np.median(err) * DAY), alpha_chi2nu=alpha_chi2(fO.alpha, fO.alpha_err),
                   s_red_s=rn["s"] * DAY, rw_rms_s=float(np.sqrt(rn["q"] * np.ptp(t)) * DAY), chi2nu_lc=fO.chi2nu,
                   t=t, tau=tau, err=err, flag=flag, alpha=alpha, alpha_err=alpha_e)
        test = (flag == 1) & (t < 450)
        if has_M and test.sum() >= 3 and ((flag == 1) & ~test).sum() >= 2 and (segO == "O2").any():
            ps = predictive_score(t, tau, err, test, X_extra=X)
            out.update(pred_score=ps["score"], P_train=ps["P_train"], D_train=ps["D_train"])
        return out
    except Exception as ex:
        return dict(ogle_id=oid, ok=False, msg=str(ex)[:200])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--out", default="results/real/oc_all.parquet")
    a = ap.parse_args()
    inv = pd.read_parquet("data/lc_inventory.parquet")
    ident = read_ogle4_ident("data/raw/ogle4_lmc_rrlyr/ident.dat")[["ogle_id", "ra", "dec"]]
    par = pd.read_fwf("data/raw/ogle4_lmc_rrlyr/RRab.dat", colspecs=[(0, 20), (37, 46), (58, 68)],
                      names=["ogle_id", "P", "T0"], header=None)
    s = inv[(inv.subtype == "RRab") & (inv.n3_I > 0) & (inv.n4_I > 0)].merge(par, on="ogle_id").merge(ident, on="ogle_id")
    if a.limit:
        s = s.sample(a.limit, random_state=0)
    jobs = list(zip(s.ogle_id, s.macho_id, s.P.astype(float), s.T0.astype(float), s.ra, s.dec))
    print(f"stars: {len(jobs)}", flush=True)
    with Pool(a.workers) as pool:
        rows = list(pool.imap_unordered(one, jobs, chunksize=8))
    df = pd.DataFrame(rows)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(a.out)
    ok = df[df.ok]
    print(f"ok {len(ok)} / {len(df)}; with MACHO {ok.has_M.sum()}; with predictive score {ok.pred_score.notna().sum() if 'pred_score' in ok else 0}")


if __name__ == "__main__":
    main()
