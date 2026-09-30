"""Level 1, simulations: for each star, simulate its light curves in every available band (OGLE I, MACHO B, MACHO R) at the
real epochs and errors, from the star's own fitted templates, with ONE delay/amplitude realization shared by all bands; then
run exactly the same timing fits as for real stars and save the season series.

Classes (inject_recover.draw): null, blazhko, jump, rwalk, jump_big, rwalk_big, ltte, and
  empirical : white per-season jitter s and random-walk phase with rms r, with (s, r) drawn jointly from the real stars'
              H0 noise fits (--noise-table, columns s_red_s, rw_rms_s), i.e. timing noise with the real population's
              distribution (it includes any real binaries' contribution to the noise fits; binaries are rare).

Usage
-----
    python scripts/level1_sims.py --n-stars 1500 --per-class null=2,blazhko=2,jump=1,rwalk=1,jump_big=2,rwalk_big=2,empirical=4,ltte=4 \
        --noise-table results/real/oc_all_cm.parquet --out results/inject/series_v3
"""
from __future__ import annotations

import argparse
import sys

sys.path.insert(0, __import__("os").path.dirname(__file__))
import chunked  # noqa: E402,F401  (single-threaded BLAS)
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from rrlbin.pipeline import fit_star, load_star, series_from_fits  # noqa: E402
from rrlbin.simulate import delay_and_amplitude, delay_random_walk, simulate_lc  # noqa: E402
from rrlbin.timing import year_labels  # noqa: E402
from inject_recover import draw  # noqa: E402
from level1_real import KEYS, sample_table  # noqa: E402

DAY = 86400.0


def one_star(args):
    oid, mid, P, T0, ra, dec, seed, per_class, noise, lag = args
    rng = np.random.default_rng(seed)
    out = []
    try:
        lcs = load_star(oid, mid, ra, dec)
        fits = fit_star(lcs, P, T0)
        if "I" not in fits:
            return [dict(ogle_id=oid, kind="error", msg="no OGLE I fit")]
        bands = list(fits)
        t_all = np.concatenate([lcs[b][0] for b in bands])
        for kind, n in per_class.items():
            for _ in range(n):
                if kind == "empirical":
                    s_s, rw_s = noise[rng.integers(len(noise))]
                    p, info = dict(pdot=rng.normal(0, 0.3)), dict(kind=kind, s_inj_s=s_s, rw_rms_s=rw_s)
                    tau, A = delay_and_amplitude(t_all, P, rng, pdot=p["pdot"])
                    if rw_s > 0:
                        tau = tau + delay_random_walk(t_all, rw_s / DAY, rng)
                    yl = year_labels(t_all)
                    jit = {y: rng.normal(0, s_s / DAY) for y in np.unique(yl)}
                    tau = tau + np.array([jit[y] for y in yl])
                else:
                    p, info = draw(kind, rng, t_all)
                    tau, A = delay_and_amplitude(t_all, P, rng, **p)
                # physical band lag: the templates are gauge-fixed (fundamental phase 0), so without this the simulated
                # MACHO-OGLE offset would be 0, whereas real stars have the measured band lag (with intrinsic scatter)
                code = {"MB": "1", "MR": "2"}
                boff = {b: rng.normal(lag[code[b]]["slope"] * P + lag[code[b]]["intercept"], lag[code[b]]["sd"]) / DAY
                        for b in bands if b in code and lag and code[b] in lag}
                info.update({f"lag{code[b]}_s": v * DAY for b, v in boff.items()})
                sim, i0 = {}, 0
                for b in bands:
                    t, m, e, seg = lcs[b]
                    f = fits[b]
                    ms, _ = simulate_lc(t, e, seg, f.coef, dict(zip(f.seg_names, f.zp)), P, T0, rng,
                                        np.sqrt(max(f.chi2nu, 1.0)), tau=tau[i0:i0 + t.size] + boff.get(b, 0.0),
                                        A=A[i0:i0 + t.size])
                    sim[b] = (t, ms, e, seg)
                    i0 += t.size
                s = series_from_fits(fit_star(sim, P, T0), P)
                row = dict(ogle_id=oid, P=P, bands="".join(sorted(bands)), **info)
                row.update({k: s[k] for k in KEYS})
                out.append(row)
    except Exception as ex:
        out.append(dict(ogle_id=oid, kind="error", msg=str(ex)[:200]))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n-stars", type=int, default=1500)
    ap.add_argument("--seed", type=int, default=21)
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--per-class", default="null=2,blazhko=2,jump=1,rwalk=1,jump_big=2,rwalk_big=2,empirical=4,ltte=4")
    ap.add_argument("--noise-table", default="results/real/oc_all_cm.parquet")
    ap.add_argument("--out", default="results/inject/series_v3")
    ap.add_argument("--chunk", type=int, default=40)
    ap.add_argument("--band-lag", default="results/calib/band_lag_v3.json", help="band-lag relations (calibrate_band_lag.py)")
    a = ap.parse_args()
    per_class = {k: int(v) for k, v in (x.split("=") for x in a.per_class.split(","))}
    import json
    lag = json.loads(Path(a.band_lag).read_text())
    nt = pd.read_parquet(a.noise_table, columns=["has_M", "s_red_s", "rw_rms_s"])
    noise = nt[nt.has_M][["s_red_s", "rw_rms_s"]].dropna().to_numpy()
    s = sample_table()
    s = s[s.macho_id.fillna("").str.count(r"\.").eq(2)]
    s = s.sample(n=min(a.n_stars, len(s)), random_state=a.seed).sort_values("ogle_id")
    jobs = [(o, mi, float(P), float(T0), ra, de, a.seed * 100000 + k, per_class, noise, lag)
            for k, (o, mi, P, T0, ra, de) in enumerate(zip(s.ogle_id, s.macho_id, s.P, s.T0, s.ra, s.dec))]
    print(f"stars {len(jobs)}; sims per class {per_class}", flush=True)
    chunked.run_chunked(one_star, jobs, a.out, workers=a.workers, chunk=a.chunk, flatten=True)
    d = chunked.merge_parts(a.out)
    print(d.kind.value_counts().to_dict())


if __name__ == "__main__":
    main()
