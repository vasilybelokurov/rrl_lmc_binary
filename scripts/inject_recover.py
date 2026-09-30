"""Injection-recovery on real OGLE cadences: nulls, timing nuisances and LTTE orbits through the full pipeline.

For each star (random RRab with OGLE-III+IV I band):
  1. fit the real light curve (template, zero points, noise scale sqrt(max(chi2nu, 1))); compute the real O-C statistic D;
  2. simulate light curves at the real epochs and errors for:
       null    : constant Pdot ~ N(0, 0.3) d/Myr only
       blazhko : + coherent phase modulation eps_phi ~ U(0.005, 0.04) cycles and amplitude modulation eps_A ~ U(0, 0.3),
                 P_B ~ logU(20, 3000) d
       jump    : + abrupt period change dP/P ~ logU(1e-6, 2e-5) at a random time
       rwalk   : + random-walk phase, rms ~ logU(50, 1000) s
       ltte    : + orbit P_orb ~ logU(300, 10000) d, M2 ~ logU(0.05, 1.5) Msun, M1 = 0.65, cos i ~ U(0, 1),
                 e = 0 (half) or ~ U(0, 0.7) (half), random omega and periastron time
     (all non-null classes also carry the random Pdot);
  3. run fit_timing -> unwrap -> oc_search, and record D, the best period and amplitude, the jitter, and the
     amplitude-constancy diagnostic chi2nu(alpha_j).

Usage
-----
    python scripts/inject_recover.py --n-stars 300 --workers 6 --out results/inject/run1.parquet
"""
from __future__ import annotations

import argparse
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from rrlbin.ltte import a1sini_over_c, oc_search, period_grid  # noqa: E402
from rrlbin.simulate import simulate_lc  # noqa: E402
from rrlbin.timing import fit_timing, unwrap_delays  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from timing_sample import ogle_I  # noqa: E402

N_PER_CLASS = dict(null=4, blazhko=4, jump=2, rwalk=2, ltte=8)
DAY = 86400.0


def analyse(t, m, e, seg, P, T0):
    f = fit_timing(t, m, e, seg, P, T0, K=8)
    tau = unwrap_delays(f.tau, P)
    periods = period_grid(np.ptp(f.t_season))
    r = oc_search(f.t_season, tau, f.tau_err, periods)
    wa = f.alpha_err ** -2
    abar = np.sum(wa * f.alpha) / wa.sum()
    return f, dict(D=r["D"], P_best=r["P_best"], amp_best_s=r["amp"] * DAY, jit0_s=r["jit0"] * DAY,
                   jit1_s=r["jit1"] * DAY, err_med_s=float(np.median(f.tau_err) * DAY), n_season=f.season.size,
                   alpha_chi2nu=float(np.sum(wa * (f.alpha - abar) ** 2) / max(f.alpha.size - 1, 1)))


def draw(kind, rng, t):
    lo, hi = t.min(), t.max()
    p = dict(pdot=rng.normal(0, 0.3))
    info = dict(kind=kind, pdot=p["pdot"])
    if kind == "blazhko":
        b = dict(P_B=float(np.exp(rng.uniform(np.log(20), np.log(3000)))), eps_A=rng.uniform(0, 0.3),
                 eps_phi=rng.uniform(0.005, 0.04), psi_A=rng.uniform(0, 2 * np.pi), psi_phi=rng.uniform(0, 2 * np.pi))
        p["blazhko"] = b
        info.update(P_B=b["P_B"], eps_A=b["eps_A"], eps_phi=b["eps_phi"])
    elif kind == "jump":
        j = dict(dP_over_P=float(np.exp(rng.uniform(np.log(1e-6), np.log(2e-5)))) * rng.choice([-1, 1]),
                 t_break=rng.uniform(lo + 0.1 * (hi - lo), hi - 0.1 * (hi - lo)))
        p["jump"] = j
        info.update(dP_over_P=j["dP_over_P"])
    elif kind == "jump_big":   # large abrupt period changes (tail of the real O-C noise)
        j = dict(dP_over_P=float(np.exp(rng.uniform(np.log(2e-5), np.log(2e-4)))) * rng.choice([-1, 1]),
                 t_break=rng.uniform(lo + 0.1 * (hi - lo), hi - 0.1 * (hi - lo)))
        p["jump"] = j
        info.update(dP_over_P=j["dP_over_P"])
    elif kind == "rwalk_big":  # large random-walk phase wander (tail of the real O-C noise)
        p["rw_rms"] = float(np.exp(rng.uniform(np.log(1000), np.log(15000)))) / DAY
        info.update(rw_rms_s=p["rw_rms"] * DAY)
    elif kind == "rwalk":
        p["rw_rms"] = float(np.exp(rng.uniform(np.log(50), np.log(1000)))) / DAY
        info.update(rw_rms_s=p["rw_rms"] * DAY)
    elif kind == "ltte":
        P_orb = float(np.exp(rng.uniform(np.log(300), np.log(10000))))
        M2 = float(np.exp(rng.uniform(np.log(0.05), np.log(1.5))))
        cosi = rng.uniform(0, 1)
        ecc = 0.0 if rng.uniform() < 0.5 else rng.uniform(0, 0.7)
        amp = a1sini_over_c(P_orb, 0.65, M2, np.sqrt(1 - cosi ** 2)) / DAY
        p["ltte"] = dict(P_orb=P_orb, amp=amp, e=ecc, omega=rng.uniform(0, 2 * np.pi), t_peri=rng.uniform(0, P_orb))
        info.update(P_orb=P_orb, M2=M2, cosi=cosi, ecc=ecc, amp_s=amp * DAY)
    return p, info


def one_star(args):
    oid, P, T0, seed = args
    rng = np.random.default_rng(seed)
    out = []
    try:
        t, m, e, seg = ogle_I(oid)
        f, real = analyse(t, m, e, seg, P, T0)
        out.append(dict(ogle_id=oid, kind="real", **real))
        zp = dict(zip(f.seg_names, f.zp))
        ns = np.sqrt(max(f.chi2nu, 1.0))
        for kind, n in N_PER_CLASS.items():
            for _ in range(n):
                p, info = draw(kind, rng, t)
                ms, _ = simulate_lc(t, e, seg, f.coef, zp, P, T0, rng, noise_scale=ns, **p)
                _, res = analyse(t, ms, e, seg, P, T0)
                out.append(dict(ogle_id=oid, **info, **res))
    except Exception as ex:
        out.append(dict(ogle_id=oid, kind="error", msg=str(ex)[:200]))
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n-stars", type=int, default=300)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--out", default="results/inject/run1.parquet")
    a = ap.parse_args()
    inv = pd.read_parquet("data/lc_inventory.parquet")
    par = pd.read_fwf("data/raw/ogle4_lmc_rrlyr/RRab.dat", colspecs=[(0, 20), (37, 46), (58, 68)],
                      names=["ogle_id", "P", "T0"], header=None)
    s = inv[(inv.subtype == "RRab") & (inv.n3_I > 0) & (inv.n4_I > 0)].merge(par, on="ogle_id")
    s = s.sample(n=a.n_stars, random_state=a.seed)
    jobs = [(o, float(P), float(T0), a.seed * 100000 + k) for k, (o, P, T0) in enumerate(zip(s.ogle_id, s.P, s.T0))]
    with Pool(a.workers) as pool:
        rows = [r for rs in pool.imap_unordered(one_star, jobs, chunksize=2) for r in rs]
    df = pd.DataFrame(rows)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(a.out)
    print(df.kind.value_counts().to_string())


if __name__ == "__main__":
    main()
