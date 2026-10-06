"""Part C1: provisional v3 candidate list and comparison with the v2 list.

Cuts (as v2; plot_summary_stats.cut_flags): D > 40, chi2_nu(alpha) < 2 (now pooled over OGLE and MACHO bands), A/sigma > 3,
>= 1.5 cycles, A below the LTTE ceiling. Additional v3 flags (reported, not cut): harmonic coherence chi2_nu, predictive
score, alias Delta D, 2-harmonic (eccentric) Delta D. Thresholds are provisional until Part B is complete.

Outputs: results/partC/candidates_v3_prov.csv, results/partC/v2_vs_v3.csv; prints a summary.

Usage
-----
    python scripts/partC_candidates.py
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import brentq

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_summary_stats import cut_flags  # noqa: E402
from rrlbin.ltte import mass_function  # noqa: E402
sys.path.insert(0, __import__('os').path.dirname(__file__))
from version import CM, LAG, PARTB, PARTC, PLOTS_C, PREV, PREV_CANDS, SERIES, SFX, SIM_M, SIM_O, STATS, V  # noqa: E402,F401

OUT = Path(PARTC)


def m2_min(f):
    return brentq(lambda m: m ** 3 / (0.65 + m) ** 2 - f, 1e-4, 100) if np.isfinite(f) and f > 0 else np.nan


def select(stats: pd.DataFrame) -> pd.DataFrame:
    s = stats.join(cut_flags(stats, "amp_s", stats.baseline))
    return s


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    v3 = select(pd.read_parquet(STATS))
    v3 = v3[v3.err_msg.isna()] if "err_msg" in v3 else v3
    c = v3[v3["all"]].copy()
    c["fM"] = mass_function(c.P_best, c.amp_s)
    c["M2min"] = [m2_min(v) for v in c.fM]
    c["K1_kms"] = 2 * np.pi * c.amp_s * 299792.458 / (c.P_best * 86400)
    c["snr_tot"] = c.amp_s / np.hypot(c.err_med_s, c.jit0_s)
    c = c.sort_values("snr_tot", ascending=False)
    cols = ["ogle_id", "has_M", "bands", "n_season", "D", "D_2h", "dD_harm2", "P_best", "amp_s", "amp_2h_s", "err_med_s", "jit0_s",
            "alpha_chi2nu", "coh_chi2nu", "pred_score", "alias_dD", "fM", "M2min", "K1_kms", "snr_tot",
            "s_gp_s", "A_gp_s", "ell_gp_d"]
    c[[k for k in cols if k in c]].to_csv(OUT / f"candidates_{V}_prov.csv", index=False)

    v2 = pd.read_csv(PREV_CANDS)
    v3all = v3.set_index("ogle_id")
    rows = []
    for oid in sorted(set(v2.ogle_id) | set(c.ogle_id)):
        r = dict(ogle_id=oid, in_v2=oid in set(v2.ogle_id), in_v3=oid in set(c.ogle_id))
        if oid in v2.ogle_id.values:
            x = v2.set_index("ogle_id").loc[oid]
            r.update(P_v2=x.P_best, A_v2=x.amp_s, D_v2=x.D)
        if oid in v3all.index:
            y = v3all.loc[oid]
            r.update(P_v3=y.P_best, A_v3=y.amp_s, D_v3=y.D, coh=y.get("coh_chi2nu"), pred=y.get("pred_score"),
                     failed=",".join(k for k in ["c1_D", "c2_alpha", "c3_snr", "c4_cycles", "c5_ceiling"] if not y[k]))
        rows.append(r)
    cmp_ = pd.DataFrame(rows)
    # internal labels v2/v3 = previous/current pipeline version (v3: v2 -> v3; v4: v3 -> v4)
    relab = lambda x: re.sub(r"v[23]", lambda m: {"v2": PREV, "v3": V}[m.group()], x)
    cmp_.rename(columns=relab).to_csv(OUT / relab("v2_vs_v3.csv"), index=False)
    import builtins
    _print = builtins.print
    print = lambda *a, **k: _print(*(relab(x) if isinstance(x, str) else x for x in a), **k)  # noqa: E731,A001

    print(f"v3 provisional candidates: {len(c)} (with MACHO {c.has_M.sum()}); v2: {len(v2)}")
    kept = cmp_[cmp_.in_v2 & cmp_.in_v3]
    print(f"v2 candidates kept in v3: {len(kept)}; dropped: {(cmp_.in_v2 & ~cmp_.in_v3).sum()}; new in v3: {(~cmp_.in_v2 & cmp_.in_v3).sum()}")
    if len(kept):
        print(f"  kept: median |dP/P| {np.median(np.abs(kept.P_v3 / kept.P_v2 - 1)):.3f}; median A_v3/A_v2 {np.median(kept.A_v3 / kept.A_v2):.2f}; "
              f"median D_v3/D_v2 {np.median(kept.D_v3 / kept.D_v2):.2f}")
    dropped = cmp_[cmp_.in_v2 & ~cmp_.in_v3]
    print("  dropped v2 candidates, failed cuts:", dropped.failed.fillna("not in v3 stats").value_counts().to_dict())
    print("v3 flags among candidates: coherence chi2_nu > 2:", int((c.coh_chi2nu > 2).sum()),
          "| predictive score available:", int(c.pred_score.notna().sum()), ", > 0:", int((c.pred_score > 0).sum()),
          "| alias ambiguous (dD < 5):", int((c.alias_dD < 5).sum()), "| eccentric preferred (dD_harm2 > 10):", int((c.dD_harm2 > 10).sum()))


if __name__ == "__main__":
    main()
