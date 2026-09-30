"""Injection-recovery with MACHO B (1992-1999) + OGLE I (OGLE-II/III/IV) on real cadences of the same stars.

Like inject_recover.py, but each simulation draws ONE delay/amplitude realization on the union of MACHO and OGLE
epochs, then synthesizes each survey from its own fitted template. Seasons are observing years (year_labels).
The O-C search runs on the combined season delays with a free MACHO offset (indicator column; constrained by the
years both surveys observed). For comparison, the OGLE-only search on the same simulation is also recorded.

Usage
-----
    python scripts/inject_recover_macho.py --n-stars 200 --workers 6 --out results/inject/macho_run1.parquet
"""
from __future__ import annotations

import argparse
import sys

sys.path.insert(0, __import__("os").path.dirname(__file__))
import chunked  # noqa: E402,F401  (sets single-threaded BLAS before numpy)
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from rrlbin.io import read_macho, read_ogle4_ident  # noqa: E402
from rrlbin.ltte import oc_search, period_grid, predictive_score  # noqa: E402
from rrlbin.simulate import delay_and_amplitude, simulate_lc  # noqa: E402
from rrlbin.timing import fit_timing, unwrap_delays, year_labels  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from inject_recover import N_PER_CLASS, draw  # noqa: E402
from timing_sample import ogle_I  # noqa: E402

DAY = 86400.0


def alpha_chi2(f):
    wa = f.alpha_err ** -2
    abar = np.sum(wa * f.alpha) / wa.sum()
    return float(np.sum(wa * (f.alpha - abar) ** 2) / max(f.alpha.size - 1, 1))


def analyse(tO, mO, eO, segO, tM, mM, eM, P, T0):
    fO = fit_timing(tO, mO, eO, segO, P, T0, K=8, labels=year_labels(tO))
    fM = fit_timing(tM, mM, eM, np.full(tM.size, "M"), P, T0, K=8, labels=year_labels(tM))
    tauO, tauM = unwrap_delays(fO.tau, P), unwrap_delays(fM.tau, P)
    out = {}
    rO = oc_search(fO.t_season, tauO, fO.tau_err, period_grid(np.ptp(fO.t_season)))
    out.update(D_O=rO["D"], P_O=rO["P_best"], amp_O_s=rO["amp"] * DAY)
    t = np.r_[fM.t_season, fO.t_season]
    tau = np.r_[tauM, tauO]
    err = np.r_[fM.tau_err, fO.tau_err]
    ind = np.r_[np.ones(tauM.size), np.zeros(tauO.size)][:, None]
    r = oc_search(t, tau, err, period_grid(np.ptp(t)), X_extra=ind)
    n_overlap = int(np.isin(fM.season, fO.season).sum())
    # out-of-sample test: hold out the MACHO seasons before 1997 (HJD' < 450, i.e. before OGLE-II), train on the rest;
    # the MACHO offset is then fixed by the overlap years in the training set
    test = (ind[:, 0] == 1) & (t < 450)
    if test.sum() >= 3 and ((ind[:, 0] == 1) & ~test).sum() >= 2 and (tO < 2000).any():
        ps = predictive_score(t, tau, err, test, X_extra=ind)
        out.update(pred_score=ps["score"], P_train=ps["P_train"], D_train=ps["D_train"], n_test=int(test.sum()))
    out.update(D=r["D"], P_best=r["P_best"], amp_best_s=r["amp"] * DAY, jit0_s=r["jit0"] * DAY,
               errO_med_s=float(np.median(fO.tau_err) * DAY), errM_med_s=float(np.median(fM.tau_err) * DAY),
               n_overlap=n_overlap, alpha_chi2nu=alpha_chi2(fO), alpha_chi2nu_M=alpha_chi2(fM),
               err_med_s=float(np.median(err) * DAY), baseline=float(np.ptp(t)),
               t=t, tau=tau, err=err, flag=ind[:, 0])   # season series saved: statistics can be recomputed without re-simulating
    return fO, fM, out


def one_star(args):
    oid, mid, P, T0, ra, dec, seed, per_class = args   # per_class passed explicitly: spawned workers re-import defaults
    rng = np.random.default_rng(seed)
    out = []
    try:
        f_, t_, s_ = map(int, mid.split("."))
        tM, mM, eM = read_macho(f"data/raw/macho/{f_}.{t_}.parquet", s_, ra, dec, band="b")
        tO, mO, eO, segO = ogle_I(oid)
        fO, fM, real = analyse(tO, mO, eO, segO, tM, mM, eM, P, T0)
        out.append(dict(ogle_id=oid, kind="real", has_O2=bool((segO == "O2").any()), **real))
        zpO, zpM = dict(zip(fO.seg_names, fO.zp)), dict(zip(fM.seg_names, fM.zp))
        nsO, nsM = np.sqrt(max(fO.chi2nu, 1.0)), np.sqrt(max(fM.chi2nu, 1.0))
        t_all = np.r_[tO, tM]
        for kind, n in per_class.items():
            for _ in range(n):
                p, info = draw(kind, rng, t_all)
                tau, A = delay_and_amplitude(t_all, P, rng, **p)
                mOs, _ = simulate_lc(tO, eO, segO, fO.coef, zpO, P, T0, rng, nsO, tau=tau[:tO.size], A=A[:tO.size])
                mMs, _ = simulate_lc(tM, eM, np.full(tM.size, "M"), fM.coef, zpM, P, T0, rng, nsM,
                                     tau=tau[tO.size:], A=A[tO.size:])
                _, _, res = analyse(tO, mOs, eO, segO, tM, mMs, eM, P, T0)
                out.append(dict(ogle_id=oid, has_O2=bool((segO == "O2").any()), **info, **res))
    except Exception as ex:
        out.append(dict(ogle_id=oid, kind="error", msg=str(ex)[:200]))
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n-stars", type=int, default=200)
    ap.add_argument("--seed", type=int, default=2)
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--out", default="results/inject/macho_run1.parquet")
    ap.add_argument("--per-class", default=None, help="sims per class, e.g. null=2,blazhko=2,jump=1,rwalk=1,ltte=4")
    ap.add_argument("--require-o2", action="store_true", help="only stars with OGLE-II epochs (for the predictive test)")
    a = ap.parse_args()
    per_class = dict(N_PER_CLASS)
    if a.per_class:
        per_class = {k: int(v) for k, v in (x.split("=") for x in a.per_class.split(","))}
    print("sims per class:", per_class, flush=True)
    inv = pd.read_parquet("data/lc_inventory.parquet")
    ident = read_ogle4_ident("data/raw/ogle4_lmc_rrlyr/ident.dat")[["ogle_id", "ra", "dec"]]
    par = pd.read_fwf("data/raw/ogle4_lmc_rrlyr/RRab.dat", colspecs=[(0, 20), (37, 46), (58, 68)],
                      names=["ogle_id", "P", "T0"], header=None)
    s = inv[(inv.subtype == "RRab") & (inv.n3_I > 0) & (inv.n4_I > 0) & (inv.macho_id.fillna("").str.count(r"\.") == 2)]
    s = s[[Path("data/raw/macho/" + ".".join(m.split(".")[:2]) + ".parquet").exists() for m in s.macho_id]]
    if a.require_o2:
        s = s[s.t3_first < 2000]
    s = s.merge(par, on="ogle_id").merge(ident, on="ogle_id")
    print(f"eligible stars (tile fetched): {len(s)}", flush=True)
    s = s.sample(n=min(a.n_stars, len(s)), random_state=a.seed)
    jobs = [(o, mi, float(P), float(T0), ra, de, a.seed * 100000 + k, per_class)
            for k, (o, mi, P, T0, ra, de) in enumerate(zip(s.ogle_id, s.macho_id, s.P, s.T0, s.ra, s.dec))]
    parts = Path(a.out).with_suffix("")
    chunked.run_chunked(one_star, jobs, parts, workers=a.workers, chunk=60, flatten=True)
    df = chunked.merge_parts(parts)
    df.to_parquet(a.out)
    print(df.kind.value_counts().to_string())


if __name__ == "__main__":
    main()
