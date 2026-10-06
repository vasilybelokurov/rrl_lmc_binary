"""Part C3: Keplerian LTTE fits of the candidates (results/partC/candidates_v3_prov.csv).

Season series from Level 1 (results/real/series_v3), common mode and band-lag priors from Part A, robust unwrapping and
whole-cycle band alignment as in Level 2. White jitter fixed at the circular-fit value (jit1_s). Keplerian and circular fits;
Delta BIC = BIC_circ - BIC_kepler (> 0 favours eccentric); f(M), M2,min (M1 = 0.65, edge-on), K1 = 2 pi c A / (P sqrt(1-e^2)).

Usage
-----
    python scripts/partC_kepler.py --workers 2
"""
from __future__ import annotations

import argparse
import json
import sys

sys.path.insert(0, __import__("os").path.dirname(__file__))
import chunked  # noqa: E402,F401
from multiprocessing import Pool  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.optimize import brentq  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from level2 import load, priors_for, to_series  # noqa: E402
from rrlbin.kepler_fit import fit_keplerian  # noqa: E402
from rrlbin.oc import align_bands, apply_common_mode, robust_unwrap  # noqa: E402
from rrlbin.timing import year_labels  # noqa: E402
sys.path.insert(0, __import__('os').path.dirname(__file__))
from version import CM, LAG, PARTB, PARTC, PLOTS_C, PREV, PREV_CANDS, SERIES, SFX, SIM_M, SIM_O, STATS, V  # noqa: E402,F401

DAY = 86400.0
G = {}


def prep(r):
    s = apply_common_mode(to_series(r), G["cm"], year_labels)
    t, y, e, band = s["t"], s["tau"], s["err"], s["band"].astype(int)
    o = np.argsort(t)
    t, y, e, band = t[o], y[o], e[o], band[o]
    pri = priors_for(s["P"], G["lag"])
    y = robust_unwrap(t, y, band, s["P"])
    y = align_bands(t, y, e, band, s["P"], pri)
    return t, y, e, band, pri, s["P"]


def one(args):
    r, P0, jit1 = args
    try:
        t, y, e, band, pri, Pp = prep(r)
        k = fit_keplerian(t, y, e, band, P0, s_jit=jit1 / DAY, priors=pri, n_boot=100)
        c = fit_keplerian(t, y, e, band, P0, s_jit=jit1 / DAY, priors=pri, n_boot=0, circular=True)
        n = k["n"]
        bic = lambda f: f["chi2"] + f["k"] * np.log(n)
        fM = 4 * np.pi ** 2 * (k["A_s"] * 299792458.0) ** 3 / (6.6743e-11 * (k["P"] * DAY) ** 2) / 1.98841e30
        m2 = brentq(lambda m: m ** 3 / (0.65 + m) ** 2 - fM, 1e-4, 100) if fM > 0 else np.nan
        out = dict(ogle_id=r["ogle_id"], P_kep=k["P"], A_kep_s=k["A_s"], e_kep=k["e"], omega_kep=k["omega"], tp_kep=k["t_p"],
                   chi2_kep=k["chi2"], chi2_circ=c["chi2"], n=n, dBIC_ecc=bic(c) - bic(k), P_circ=c["P"], A_circ_s=c["A_s"],
                   fM_kep=fM, M2min_kep=m2, K1_kep=2 * np.pi * k["A_s"] * 299792.458 / (k["P"] * DAY * np.sqrt(1 - k["e"] ** 2)),
                   chi2nu_kep=k["chi2"] / max(n - k["k"], 1))
        out.update({kk: v for kk, v in k.items() if kk.endswith(("_p16", "_p50", "_p84"))})
        return out
    except Exception as ex:
        return dict(ogle_id=r["ogle_id"], err_msg=str(ex)[:200])


def _init(cm, lag):
    G.update(cm=cm, lag=lag)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--cands", default=f"{PARTC}/candidates_{V}_prov.csv")
    a = ap.parse_args()
    c = pd.read_csv(a.cands)
    st = pd.read_parquet(STATS, columns=["ogle_id", "jit1_s"]).set_index("ogle_id")
    ser = load(SERIES)
    ser = ser[ser.ogle_id.isin(set(c.ogle_id))].set_index("ogle_id")
    cm = {tuple(map(int, k.split(","))): v / DAY for k, v in json.loads(Path(CM).read_text()).items()}
    lag = json.loads(Path(LAG).read_text())
    jobs = [(dict(ser.loc[o].to_dict(), ogle_id=o), float(P0), float(st.loc[o, "jit1_s"])) for o, P0 in zip(c.ogle_id, c.P_best)]
    with Pool(a.workers, initializer=_init, initargs=(cm, lag)) as pool:
        rows = pool.map(one, jobs, chunksize=1)
    R = pd.DataFrame(rows)
    R.to_csv(f"{PARTC}/kepler.csv", index=False)
    ok = R[R.err_msg.isna()] if "err_msg" in R else R
    print(f"fitted {len(ok)} / {len(R)}")
    q = lambda x: np.nanpercentile(x, [10, 50, 90]).round(2)
    print("e p10/50/90", q(ok.e_kep), "| dBIC_ecc > 6 (eccentric preferred):", int((ok.dBIC_ecc > 6).sum()),
          "| P_kep/P_circ p10/50/90", q(ok.P_kep / ok.P_circ), "| A_kep/A_circ", q(ok.A_kep_s / ok.A_circ_s))
    print("M2,min p10/50/90", q(ok.M2min_kep), "| K1 p10/50/90", q(ok.K1_kep), "| chi2_nu (with jitter) p10/50/90", q(ok.chi2nu_kep))


if __name__ == "__main__":
    main()
