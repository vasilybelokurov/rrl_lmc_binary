"""How much OGLE-IV LMC photometry exists after 2016, judging from the public light curves that extend beyond it?

For the RR Lyrae with public OGLE-IV I-band data after HJD' 7600 (the 2017-2019 catalogue extensions): epochs per observing
year 2010-2020, their sky positions and OGLE-IV fields, and the per-season timing precision they deliver (Level-1 season errors
from results/real/series_v3). Compare with the same years for the candidates' fields (positions of the candidates).

Outputs: results/partC/ogle_post2016_cadence.csv, plots/partC/ogle_post2016.png

Usage: python scripts/ogle_post2016_cadence.py
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from level2 import load  # noqa: E402
from rrlbin.io import lc_path, read_lc, read_ogle4_ident  # noqa: E402
from rrlbin.timing import year_labels  # noqa: E402

RAW = Path("data/raw")


def main():
    inv = pd.read_parquet("data/lc_inventory.parquet")
    ident = read_ogle4_ident(RAW / "ogle4_lmc_rrlyr/ident.dat")[["ogle_id", "ra", "dec", "ogle4_field_id"]]
    inv = inv.merge(ident, on="ogle_id")
    late = inv[inv.t4_last > 7600].copy()
    cand = pd.read_csv("results/partC/partC_tiers.csv")[["ogle_id", "tier"]].merge(ident, on="ogle_id")
    years = np.arange(year_labels(np.array([5300.0]))[0], year_labels(np.array([8924.0]))[0] + 1)
    rows = []
    for oid in late.ogle_id:
        t = read_lc(lc_path(RAW, "ogle4", oid))[0]
        y = year_labels(t)
        rows.append({int(k): int((y == k).sum()) for k in years} | {"ogle_id": oid})
    E = pd.DataFrame(rows).fillna(0).set_index("ogle_id")
    E.to_csv("results/partC/ogle_post2016_cadence.csv")
    from astropy.time import Time
    yr0 = Time(2450245.0, format="jd").decimalyear          # year label 0 starts at HJD' 245 (1996.4)
    lab = {k: f"{yr0 + k:.1f}" for k in years}
    print("OGLE-IV epochs per observing year (label = start of the year), stars with public data after 2016 (N = %d):" % len(E))
    for k in years:
        col = E[k]
        print(f"  {lab[k]}: stars with data {int((col > 0).sum()):5d}; median epochs (if >0) {col[col > 0].median() if (col > 0).any() else 0:.0f}")
    nseas = (E[list(years)] > 0).sum(axis=1)
    print("number of OGLE-IV observing years per late star: p10/50/90", np.percentile(nseas, [10, 50, 90]))
    # field overlap: which OGLE-IV fields host the late stars vs the candidates
    late["field"] = late.ogle4_field_id.str.extract(r"^(LMC\d+)", expand=False)
    cand["field"] = cand.ogle4_field_id.str.extract(r"^(LMC\d+)", expand=False)
    shared = set(late.field) & set(cand.field)
    print(f"late-data stars in {late.field.nunique()} fields; candidates in {cand.field.nunique()} fields; shared fields: {len(shared)} {sorted(shared)[:10]}")
    rr = np.hypot((late.ra - 80.9) * np.cos(np.radians(-69.75)), late.dec + 69.75)
    rc = np.hypot((cand.ra - 80.9) * np.cos(np.radians(-69.75)), cand.dec + 69.75)
    print(f"distance from LMC centre [deg]: late-data stars p10/50/90 {np.percentile(rr, [10, 50, 90]).round(1)}; "
          f"candidates {np.percentile(rc, [10, 50, 90]).round(1)}")
    # timing precision per season for late stars (Level-1 series, OGLE I band)
    ser = load("results/real/series_v3")
    ser = ser[ser.ogle_id.isin(set(late.ogle_id))]
    pre, post = [], []
    for r in ser.itertuples():
        t, e, b = np.asarray(r.t), np.asarray(r.err), np.asarray(r.band)
        m = b == 0
        pre += list(e[m & (t > 5200) & (t < 7550)] * 86400)
        post += list(e[m & (t > 7550)] * 86400)
    print(f"OGLE-IV per-season delay error for these stars: 2010-2016 median {np.median(pre):.0f} s (N={len(pre)}), "
          f"2016-2020 median {np.median(post):.0f} s (N={len(post)})")
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.5))
    q = E[list(years)]
    ax[0].bar([yr0 + k + 0.5 for k in years], (q > 0).sum(), width=0.8, color="C0", alpha=0.6, label="stars with data")
    ax2 = ax[0].twinx()
    ax2.plot([yr0 + k + 0.5 for k in years], [q[k][q[k] > 0].median() if (q[k] > 0).any() else 0 for k in years], "ko-",
             label="median epochs (stars with data)")
    ax2.set_ylabel("median epochs per year")
    ax[0].set(xlabel="observing year", ylabel="stars with OGLE-IV data", title=f"RRL with public data after 2016 (N={len(E)})")
    ax[0].legend()
    ax[1].scatter(inv.ra, inv.dec, s=0.3, c="0.8", rasterized=True, label="all OGLE-IV LMC RRL")
    ax[1].scatter(late.ra, late.dec, s=2, c="C0", label="public data after 2016")
    ax[1].scatter(cand.ra, cand.dec, s=25, c="C3", marker="*", label="candidates")
    ax[1].invert_xaxis()
    ax[1].set(xlabel="RA [deg]", ylabel="Dec [deg]")
    ax[1].legend(fontsize=8, markerscale=3)
    fig.tight_layout()
    fig.savefig("plots/partC/ogle_post2016.png", dpi=110)


if __name__ == "__main__":
    main()


def candidate_fields_cadence():
    """Epochs per observing year after 2016 for public late-data stars located in the candidates' OGLE-IV fields, and
    the per-season delay precision expected for the candidates from their own 2010-2016 OGLE-IV seasons, scaled by
    sqrt(N_2010-16 per season / N_post-2016 per season)."""
    from astropy.time import Time
    inv = pd.read_parquet("data/lc_inventory.parquet")
    ident = read_ogle4_ident(RAW / "ogle4_lmc_rrlyr/ident.dat")[["ogle_id", "ogle4_field_id"]]
    inv = inv.merge(ident, on="ogle_id")
    inv["field"] = inv.ogle4_field_id.str.extract(r"^(LMC\d+)", expand=False)
    E = pd.read_csv("results/partC/ogle_post2016_cadence.csv").set_index("ogle_id")
    E.columns = E.columns.astype(int)
    cand = pd.read_csv("results/partC/partC_tiers.csv")[["ogle_id", "tier"]].merge(inv[["ogle_id", "field"]], on="ogle_id")
    yr0 = Time(2450245.0, format="jd").decimalyear
    rows = []
    for fld, g in cand.groupby("field"):
        late = inv[(inv.field == fld) & (inv.t4_last > 7600)].ogle_id
        late = [o for o in late if o in E.index]
        ep = E.loc[late] if late else None
        r = dict(field=fld, n_cand=len(g), n_tier1=int((g.tier == 1).sum()), n_late_stars=len(late))
        for k in (21, 22, 23):
            r[f"epochs_{yr0 + k:.1f}"] = float(ep[k][ep[k] > 0].median()) if ep is not None and (ep[k] > 0).any() else 0.0
        rows.append(r)
    F = pd.DataFrame(rows).sort_values("n_cand", ascending=False)
    print(F.to_string(index=False))
    # candidates' own OGLE-IV 2010-2016 epochs per season and season errors
    ser = load("results/real/series_v3")
    ser = ser[ser.ogle_id.isin(set(cand.ogle_id))]
    nI, eI = [], []
    for r in ser.itertuples():
        tI = read_lc(lc_path(RAW, "ogle4", r.ogle_id))[0]
        yl = year_labels(tI)
        n_per = np.median(np.bincount(yl - yl.min())[np.bincount(yl - yl.min()) > 0])
        t, e, b = np.asarray(r.t), np.asarray(r.err), np.asarray(r.band)
        m = (b == 0) & (t > 5200)
        nI.append(n_per)
        eI.append(np.median(e[m]) * 86400)
    nI, eI = np.array(nI), np.array(eI)
    print(f"\ncandidates' OGLE-IV 2010-2016: median epochs per season {np.median(nI):.0f}, median season delay error {np.median(eI):.0f} s")
    for n_post in (106, 39, 5):
        print(f"  a post-2016 season with {n_post} epochs -> expected delay error ~ {np.median(eI * np.sqrt(nI / n_post)):.0f} s")


if __name__ == "__main__":
    candidate_fields_cadence()
