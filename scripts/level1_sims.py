"""Level 1, simulations: for each star, simulate its light curves in every available band (OGLE I, MACHO B, MACHO R) at the
real epochs and errors, from the star's own fitted templates, with ONE delay/amplitude realization shared by all bands; then
run exactly the same timing fits as for real stars and save the season series.

Classes (inject_recover.draw): null, blazhko, jump, rwalk, jump_big, rwalk_big, ltte, and
  empirical : white per-season jitter s and random-walk phase with rms r, with (s, r) drawn jointly from the real stars'
              H0 noise fits (--noise-table, columns s_red_s, rw_rms_s; candidates excluded with --exclude). The jitter
              variance is split 64% common to all bands / 36% per instrument group (MACHO B+R shared; OGLE separate), as
              measured from same-season residual correlations of real stars; the random walk is common to all bands.

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
from rrlbin.simulate import delay_and_amplitude, delay_gp, simulate_lc  # noqa: E402
from rrlbin.timing import year_labels  # noqa: E402
from inject_recover import draw  # noqa: E402
from level1_real import KEYS, sample_table  # noqa: E402

DAY = 86400.0
F_COMMON = 0.64      # fraction of the empirical season-jitter variance common to all bands (see noise_origin_test.py)


def one_star(args):
    oid, mid, P, T0, ra, dec, seed, per_class, noise, lag = args
    rng = np.random.default_rng(seed)
    out = []
    try:
        lcs = load_star(oid, mid, ra, dec)
        fails = {}
        fits = fit_star(lcs, P, T0, failures=fails)
        if "I" not in fits:
            return [dict(ogle_id=oid, kind="error", msg="no OGLE I fit")]
        bands = list(fits)
        t_all = np.concatenate([lcs[b][0] for b in bands])
        for kind, n in per_class.items():
            for _ in range(n):
                if kind == "empirical":
                    # noise drawn RELATIVE to the per-season timing error (real star's s/sigma, rw/sigma), rescaled by
                    # this star's own median season error: preserves the real distribution of normalized scatter
                    # smooth wander (squared-exponential GP, common to all bands) + white season jitter, with
                    # (s/sigma, A/sigma, l) drawn jointly from the real stars' REML noise fits (gp_null); real excess
                    # noise is red, not Brownian (scripts/noise_timescale_test.py)
                    sig = np.median(np.concatenate([fits[b].tau_err for b in bands])) * DAY
                    s_rel, A_rel, ell = noise[rng.integers(len(noise))]
                    s_s, A_s = s_rel * sig, A_rel * sig
                    p, info = dict(pdot=rng.normal(0, 0.3)), dict(kind=kind, s_inj_s=s_s, A_gp_inj_s=A_s, ell_inj_d=ell)
                    tau, A = delay_and_amplitude(t_all, P, rng, pdot=p["pdot"])
                    tau = tau + delay_gp(t_all, A_s / DAY, ell, rng)
                    # season jitter split as measured on real stars (scripts/noise_origin_test.py: same-season residual
                    # correlation MACHO B-R 0.58, MACHO-OGLE 0.37 -> ~64% of the excess variance common to all bands,
                    # ~36% shared within an instrument group (MACHO B+R) but independent between OGLE and MACHO)
                    yl = year_labels(t_all)
                    grp = np.concatenate([np.full(lcs[b][0].size, "M" if b.startswith("M") else "O") for b in bands])
                    jc = {y: rng.normal(0, np.sqrt(F_COMMON) * s_s / DAY) for y in np.unique(yl)}
                    jg = {(g, y): rng.normal(0, np.sqrt(1 - F_COMMON) * s_s / DAY) for g in ("M", "O") for y in np.unique(yl)}
                    tau = tau + np.array([jc[y] + jg[(g, y)] for g, y in zip(grp, yl)])
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
    ap.add_argument("--sample", choices=["macho", "ogle_only"], default="macho",
                    help="simulate stars with MACHO data, or OGLE-only stars (separate calibration; noise table from the same group)")
    ap.add_argument("--exclude", default=None, help="CSV with ogle_id column: stars excluded from the empirical noise table")
    a = ap.parse_args()
    per_class = {k: int(v) for k, v in (x.split("=") for x in a.per_class.split(","))}
    import json
    lag = json.loads(Path(a.band_lag).read_text())
    nt = pd.read_parquet(a.noise_table, columns=["ogle_id", "has_M", "s_gp_s", "A_gp_s", "ell_gp_d", "err_med_s"])
    if a.exclude:   # e.g. the candidate list: their fitted 'noise' may contain the orbital signal (Codex review)
        ex = set(pd.read_csv(a.exclude).ogle_id)
        nt = nt[~nt.ogle_id.isin(ex)]
        print(f"noise table: excluded {len(ex)} stars", flush=True)
    nt = nt[nt.has_M if a.sample == "macho" else ~nt.has_M].dropna(subset=["s_gp_s", "A_gp_s", "ell_gp_d", "err_med_s"])
    noise = np.column_stack([nt.s_gp_s / nt.err_med_s, nt.A_gp_s / nt.err_med_s, nt.ell_gp_d])   # relative to the season error
    s = sample_table()
    has_mid = s.macho_id.fillna("").str.count(r"\.").eq(2)
    if a.sample == "macho":
        s = s[has_mid]
    else:   # OGLE-only: no MACHO light curve (no ID, or tile/epochs missing -> check the real Level-1 band list if available)
        s = s[~has_mid]
    s = s.sample(n=min(a.n_stars, len(s)), random_state=a.seed).sort_values("ogle_id")
    jobs = [(o, mi, float(P), float(T0), ra, de, a.seed * 100000 + k, per_class, noise, lag)
            for k, (o, mi, P, T0, ra, de) in enumerate(zip(s.ogle_id, s.macho_id, s.P, s.T0, s.ra, s.dec))]
    print(f"stars {len(jobs)}; sims per class {per_class}", flush=True)
    chunked.run_chunked(one_star, jobs, a.out, workers=a.workers, chunk=a.chunk, flatten=True,
                        manifest=dict(script="level1_sims", sample=a.sample, ids=[j[0] for j in jobs], per_class=per_class, seed=a.seed,
                                      band_lag=a.band_lag, noise_table=a.noise_table, exclude=a.exclude))
    d = chunked.merge_parts(a.out)
    print(d.kind.value_counts().to_dict())


if __name__ == "__main__":
    main()
