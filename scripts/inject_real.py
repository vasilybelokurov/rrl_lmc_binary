"""Detection efficiency of the binary search measured on the REAL stars: Keplerian LTTE orbits are added to each real star's
measured season delays (all bands; the light-travel delay is achromatic), keeping its real red timing noise, its real
amplitude series (alpha) and its real cadence and errors; the Level-2 statistics are recomputed with the same code and the
candidate cuts applied (rrlbin.oc.oc_stats, plot_summary_stats.cut_flags).

Approximation: the orbit is added at each season's (Fisher-weighted) epoch instead of to the light curve, i.e. the within-season
smearing of the delay is neglected (relative amplitude loss ~ (2 pi sigma_t / P)^2 / 2 with sigma_t ~ 70 d: 10% at P = 1000 d,
1% at P = 3000 d). Validated by --target sims: the same series-level injection into simulated stars WITHOUT timing noise
(class 'null' of the v4 light-curve simulations), compared with the light-curve-level injections (class 'ltte').

Orbit prior (as the light-curve simulations): P ~ logU(300, 1e4) d, M2 ~ logU(0.05, 1.5) Msun, M1 = 0.65, cos i ~ U(0, 1),
e = 0 (half) or U(0, 0.7), omega and periastron uniform.
Output: results/partB_v4/inject_real.parquet (or inject_null_sims.parquet): one row per injection with the injected
parameters, the statistics and the cut flags.
Usage:
    python scripts/inject_real.py --n-per-star 2 --workers 6
    python scripts/inject_real.py --target sims --n-per-star 2 --workers 6
"""
from __future__ import annotations

import argparse
import json
import sys

sys.path.insert(0, __import__("os").path.dirname(__file__))
import chunked  # noqa: E402  (single-threaded BLAS)
from multiprocessing import Pool  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from level2 import load, priors_for, to_series  # noqa: E402
from plot_summary_stats import cut_flags  # noqa: E402
from rrlbin.ltte import a1sini_over_c, ltte_delay  # noqa: E402
from rrlbin.oc import apply_common_mode, oc_stats  # noqa: E402
from rrlbin.timing import year_labels  # noqa: E402

DAY = 86400.0
G = {}


def draw_orbit(rng):
    P = float(np.exp(rng.uniform(np.log(300), np.log(1e4))))
    M2 = float(np.exp(rng.uniform(np.log(0.05), np.log(1.5))))
    cosi = rng.uniform(0, 1)
    e = 0.0 if rng.uniform() < 0.5 else rng.uniform(0, 0.7)
    A = float(a1sini_over_c(P, 0.65, M2, np.sqrt(1 - cosi ** 2)))
    return dict(P_orb=P, M2=M2, cosi=cosi, ecc=e, inj_amp_s=A, omega=rng.uniform(0, 2 * np.pi), t_peri=rng.uniform(0, P))


def work(job):
    rec, seed, n = job
    rng = np.random.default_rng(seed)
    out = []
    try:
        s0 = to_series(rec)
        if G["cm"]:
            s0 = apply_common_mode(s0, G["cm"], year_labels)
        pri = priors_for(s0["P"], G["lag"])
        for k in range(n):
            o = draw_orbit(rng)
            s = dict(s0)
            s["tau"] = s0["tau"] + ltte_delay(s0["t"], o["P_orb"], o["inj_amp_s"] / DAY, o["ecc"], o["omega"], o["t_peri"])
            st = oc_stats(s, priors=pri, with_red=False)
            st.pop("tau_unwrapped", None)
            out.append(dict(ogle_id=rec["ogle_id"], rep=k, **o, **st))
    except Exception as ex:
        out.append(dict(ogle_id=rec["ogle_id"], err_msg=str(ex)[:200]))
    return out


def _init(cm, lag):
    G.update(cm=cm, lag=lag)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--target", choices=["real", "sims"], default="real")
    ap.add_argument("--n-per-star", type=int, default=2)
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--seed", type=int, default=31)
    ap.add_argument("--max-stars", type=int, default=None)
    a = ap.parse_args()
    lag = json.loads(Path("results/calib/band_lag_v4.json").read_text())
    if a.target == "real":
        d = load("results/real/series_v4")
        d = d[d.ok.astype(bool) & d.bands.str.contains("M")]
        cm = {tuple(map(int, k.split(","))): v / DAY for k, v in json.loads(Path("results/calib/common_mode_v4.json").read_text()).items()}
        out = "results/partB_v4/inject_real.parquet"
    else:   # simulated stars with photometric noise only (light-curve level), no injected timing signal
        d = chunked.merge_parts("results/inject/series_v4_macho")
        d = d[d.kind == "null"].drop_duplicates("ogle_id")
        cm = {}
        out = "results/partB_v4/inject_null_sims.parquet"
    d = d[d.t.apply(len) >= 6]
    if a.max_stars:
        d = d.sample(a.max_stars, random_state=a.seed)
    recs = d.to_dict("records")
    jobs = [(r, a.seed * 100000 + i, a.n_per_star) for i, r in enumerate(recs)]
    print(f"target {a.target}: stars {len(jobs)}, injections {len(jobs) * a.n_per_star}", flush=True)
    with Pool(a.workers, initializer=_init, initargs=(cm, lag)) as pool:
        rows = []
        for i, r in enumerate(pool.imap_unordered(work, jobs, chunksize=4)):
            rows += r
            if (i + 1) % 500 == 0:
                print(f"{i + 1}/{len(jobs)} stars", flush=True)
    R = pd.DataFrame(rows)
    ok = R.err_msg.isna() if "err_msg" in R else np.ones(len(R), bool)
    R = R.join(cut_flags(R[ok], "amp_s", R[ok].baseline))
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    R.to_parquet(out)
    print(f"written {out}: {len(R)} rows; failures {int((~ok).sum())}; overall pass {R['all'].mean():.3f}")


if __name__ == "__main__":
    main()
