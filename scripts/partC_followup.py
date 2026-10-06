"""Part C5: follow-up feasibility of the candidates from their Keplerian fits (results/partC/kepler.csv).

Predicted centre-of-mass radial velocity of the RR Lyrae: v_r(t) = c dtau_LTTE/dt (self-consistent with the timing fit),
evaluated over 2027-2030; reported: K1, the predicted systemic-velocity change over that window, its sign, and the
next O-C turning point. V magnitude from OGLE-IV. Pulsation RV amplitude for RRab ~60-70 km/s (Borissova+2004;
Haschke+2012, docs/reviews/partC_literature.md), so the orbit requires per-epoch precision well below K1 after modelling the
pulsation curve.

Usage
-----
    python scripts/partC_followup.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.time import Time

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from rrlbin.ltte import ltte_delay  # noqa: E402
sys.path.insert(0, __import__('os').path.dirname(__file__))
from version import CM, LAG, PARTB, PARTC, PLOTS_C, PREV, PREV_CANDS, SERIES, SFX, SIM_M, SIM_O, STATS, V  # noqa: E402,F401

C_KMS = 299792.458
DAY = 86400.0


def rv_from_delay(t, P, A_d, e, w, tp, h=0.5):
    """Centre-of-mass RV [km/s] of the pulsator: c * d tau / dt (tau in days, t in days)."""
    return C_KMS * (ltte_delay(t + h, P, A_d, e, w, tp) - ltte_delay(t - h, P, A_d, e, w, tp)) / (2 * h)


def hjdp(year):
    return Time(f"{year}-01-01T00:00:00", scale="utc").jd - 2450000.0


def main():
    k = pd.read_csv(f"{PARTC}/kepler.csv")
    k = k[k.err_msg.isna()] if "err_msg" in k else k
    v = pd.read_fwf("data/raw/ogle4_lmc_rrlyr/RRab.dat", colspecs=[(0, 20), (22, 28), (29, 35)], names=["ogle_id", "I", "V"], header=None)
    for col in ("I", "V"):
        v[col] = pd.to_numeric(v[col], errors="coerce")
    k = k.merge(v, on="ogle_id", how="left")
    t0, t1 = hjdp(2027), hjdp(2031)
    tt = np.linspace(t0, t1, 400)
    rows = []
    for r in k.itertuples():
        rv = rv_from_delay(tt, r.P_kep, r.A_kep_s / DAY, r.e_kep, r.omega_kep, r.tp_kep)
        tau = ltte_delay(np.linspace(t0, t0 + r.P_kep, 2000), r.P_kep, r.A_kep_s / DAY, r.e_kep, r.omega_kep, r.tp_kep)
        rows.append(dict(ogle_id=r.ogle_id, V=r.V, I=r.I, P_orb_yr=r.P_kep / 365.25, e=r.e_kep, K1_kms=r.K1_kep,
                         dv_2027_2030=float(rv.max() - rv.min()), rv_2027=float(rv[0]), rv_2030=float(rv[-1]),
                         next_tau_max_yr=2027 + np.argmax(tau) / 2000 * r.P_kep / 365.25,
                         precision_needed_kms=r.K1_kep / 3))
    F = pd.DataFrame(rows).sort_values("dv_2027_2030", ascending=False)
    F.to_csv(f"{PARTC}/followup.csv", index=False)
    print(F.head(12).round(2).to_string(index=False))
    print(f"\ncandidates with predicted systemic-velocity change > 5 km/s over 2027-2030: {(F.dv_2027_2030 > 5).sum()} / {len(F)}; "
          f"median V {F.V.median():.2f}")


if __name__ == "__main__":
    main()
