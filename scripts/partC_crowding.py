"""Part C2: crowding / blending indicators of the candidates vs their parent population.

Indicators (companion project rrl_lmc_rotation, Gaia DR3 match of OGLE-IV RRL): n_gaia_within (Gaia sources within 2"),
local Gaia source density sigma_all (30"), RUWE, ipd_frac_multi_peak, G - I offset from the population median (blends that
Gaia resolves inflate the OGLE flux), and the OGLE I-band amplitude (blending dilutes it). Each candidate is compared with
the parent RRab of the SAME group (MACHO / OGLE-only; MACHO fields are the crowded centre): percentile rank, and a KS
test of candidates vs parent. Null check: KS between two random parent subsamples of the candidates' size.

Usage
-----
    python scripts/partC_crowding.py
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import ks_2samp
import sys  # noqa: E402
sys.path.insert(0, __import__('os').path.dirname(__file__))
from version import CM, LAG, PARTB, PARTC, PLOTS_C, PREV, PREV_CANDS, SERIES, SFX, SIM_M, SIM_O, STATS, V  # noqa: E402,F401

COMP = Path("../rrl_lmc_rotation/data")
OUT = Path(PARTC)
COLS = ["n_gaia_within", "sigma_all", "ruwe", "ipd_frac_multi_peak", "dGI", "amp_I"]


def load():
    g = pd.read_parquet(COMP / "ogle4_lmc_rrl_gaiadr3.parquet",
                        columns=["ogle_id", "n_gaia_within", "ruwe", "ipd_frac_multi_peak", "phot_g_mean_mag", "I", "amp_I", "matched"])
    dens = pd.read_parquet(COMP / "ogle4_lmc_rrl_density.parquet", columns=["ogle_id", "sigma_all"])
    g = g.merge(dens, on="ogle_id", how="left")
    gi = g.phot_g_mean_mag - g.I
    g["dGI"] = gi - np.nanmedian(gi)
    st = pd.read_parquet(STATS, columns=["ogle_id", "has_M"])
    return st.merge(g, on="ogle_id", how="left")


def main():
    rng = np.random.default_rng(0)
    p = load()
    c = pd.read_csv(OUT / f"candidates_{V}_prov.csv")[["ogle_id"]]
    p["cand"] = p.ogle_id.isin(set(c.ogle_id))
    rows = []
    for grp, g in p.groupby("has_M"):
        cg, pg = g[g.cand], g[~g.cand]
        for col in COLS:
            a, b = cg[col].dropna(), pg[col].dropna()
            if len(a) < 3:
                continue
            ks = ks_2samp(a, b)
            null = [ks_2samp(*(lambda x: (x[:len(a)], x[len(a):2 * len(a)]))(rng.permutation(b.to_numpy()))).pvalue for _ in range(200)]
            rows.append(dict(group="MACHO" if grp else "OGLE-only", var=col, n_cand=len(a), med_cand=a.median(), med_parent=b.median(),
                             ks_p=ks.pvalue, null_frac_below=np.mean(np.array(null) < ks.pvalue)))
        # per-candidate percentile ranks within the group
        for col in COLS:
            pct = g[col].rank(pct=True)
            p.loc[g.index, f"pct_{col}"] = pct
    R = pd.DataFrame(rows)
    print(R.round(3).to_string(index=False))
    cand = p[p.cand].copy()
    cand["flag_crowded"] = (cand.n_gaia_within >= 3) | (cand.pct_sigma_all > 0.9) | (cand.ruwe > 1.4) | (cand.ipd_frac_multi_peak > 2) \
        | (cand.dGI.abs() > 0.3)
    cand.drop(columns=["cand"]).to_csv(OUT / "crowding.csv", index=False)
    print(f"\ncandidates with any crowding flag (>=3 Gaia within 2\", Sigma in top 10%, RUWE>1.4, ipd multipeak>2, |dGI|>0.3): "
          f"{int(cand.flag_crowded.sum())} / {len(cand)}")
    fig, ax = plt.subplots(2, 3, figsize=(15, 8))
    for a, col in zip(ax.ravel(), COLS):
        for grp, ls in [(True, "-"), (False, "--")]:
            g = p[p.has_M == grp]
            bins = np.histogram_bin_edges(g[col].dropna(), 40)
            a.hist(g[~g.cand][col].dropna(), bins, density=True, histtype="step", color="0.4", ls=ls,
                   label=f"parent ({'MACHO' if grp else 'OGLE-only'})")
            a.hist(g[g.cand][col].dropna(), bins, density=True, histtype="step", color="C3", ls=ls, lw=1.5,
                   label=f"candidates ({'MACHO' if grp else 'OGLE-only'})")
        a.set(xlabel=col, ylabel="density")
    ax[0, 0].legend(fontsize=7)
    fig.tight_layout()
    Path(PLOTS_C).mkdir(parents=True, exist_ok=True)
    fig.savefig(f"{PLOTS_C}/crowding.png", dpi=110)


if __name__ == "__main__":
    main()
