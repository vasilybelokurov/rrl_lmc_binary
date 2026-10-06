"""Stage 1 summary: separability of LTTE orbits from quasi-periodic (QP) intrinsic modulation (results/stage1/ident.parquet).

For each coherence bound c_max of the alternative H_QP(c <= c_max):
  Lambda = lnL_LTTE - max_{c <= c_max} lnL_QP(c)
  threshold Lambda_q = the (1 - q) quantile of Lambda over QP simulations with true c <= c_max (q = 1%, 5%; all P, snr)
  power = P(Lambda > Lambda_q | LTTE) per (P, e, snr)
  leakage = P(Lambda > Lambda_q | QP with true c > c_max) (signals more coherent than the bound pass as orbits)
Output: results/stage1/ident_summary.csv, plots/stage1_identifiability.png; printed tables.
Usage: python scripts/identifiability_summary.py
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

COH = [0.5, 1.0, 2.0, 4.0, np.inf]


def main():
    R = pd.read_parquet("results/stage1/ident.parquet")
    if "err_msg" in R:
        print("failures:", int(R.err_msg.notna().sum()))
        R = R[R.err_msg.isna()]
    rows = []
    for cmax in COH:
        cols = [f"lnl_qp_c{c}" for c in COH if c <= cmax]
        lam = R.lnl_ltte - R[cols].max(axis=1)
        qp_in = (R.kind == "qp") & (R.c <= cmax)
        qp_out = (R.kind == "qp") & (R.c > cmax)
        for q in (0.01, 0.05):
            thr = np.quantile(lam[qp_in], 1 - q)
            for (P, e, snr), g in R[R.kind == "ltte"].groupby(["P", "e", "snr"]):
                rows.append(dict(c_max=cmax, q=q, thr=thr, P=P, e=e, snr=snr, power=float((lam[g.index] > thr).mean()),
                                 leak=float((lam[qp_out] > thr).mean()) if qp_out.any() else np.nan, n_qp_in=int(qp_in.sum())))
    S = pd.DataFrame(rows)
    Path("results/stage1").mkdir(parents=True, exist_ok=True)
    S.to_csv("results/stage1/ident_summary.csv", index=False)
    for q in (0.01, 0.05):
        T = S[S.q == q].pivot_table(index=["P", "e", "snr"], columns="c_max", values="power")
        print(f"\npower at false-alarm {q:.0%} vs coherence bound c_max (rows: P, e, snr):\n{T.round(2).to_string()}")
        L = S[S.q == q].groupby("c_max")[["thr", "leak"]].first()
        print(f"thresholds and leakage of more-coherent QP:\n{L.round(2).to_string()}")
    # figure: Lambda distributions for c_max = 1 (QP by true c; LTTE by e)
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    for k, cmax in enumerate((1.0, np.inf)):
        cols = [f"lnl_qp_c{c}" for c in COH if c <= cmax]
        lam = R.lnl_ltte - R[cols].max(axis=1)
        bins = np.linspace(np.percentile(lam, 1), np.percentile(lam, 99.5), 40)
        for c in (0.5, 1.0, 2.0, np.inf):
            ax[k].hist(lam[(R.kind == "qp") & (R.c == c)], bins, histtype="step", label=f"QP true c = {c}")
        for e in (0.0, 0.5):
            ax[k].hist(lam[(R.kind == "ltte") & (R.e == e)], bins, histtype="stepfilled", alpha=0.3, label=f"LTTE e = {e}")
        ax[k].set(xlabel="Λ = lnL(LTTE) − lnL(QP, c ≤ c_max)", ylabel="simulations", title=f"alternative QP with coherence c ≤ {cmax}")
        ax[k].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig("plots/stage1_identifiability.png", dpi=110)
    print("figure: plots/stage1_identifiability.png")


if __name__ == "__main__":
    main()
