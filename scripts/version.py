"""Pipeline version for the analysis scripts (Part B, Part C): RRL_VERSION=v3 (default; 1992-2016 public data) or v4
(1992-2026 with the OGLE-team OGLE-IV files). Inputs are versioned file names; v4 outputs go to *_v4 directories so that the
v3 products (and the frozen predictions, which depend on them) are never overwritten."""
import os

V = os.environ.get("RRL_VERSION", "v3")
SFX = "" if V == "v3" else f"_{V}"
PREV = {"v3": "v2", "v4": "v3"}.get(V, "v3")
STATS = f"results/real/stats_{V}.parquet"
SERIES = f"results/real/series_{V}"
CM = f"results/calib/common_mode_{V}.json"
LAG = f"results/calib/band_lag_{V}.json"
SIM_M = f"results/inject/stats_{V}_macho.parquet"
SIM_O = f"results/inject/stats_{V}_ogle.parquet"
PARTB = f"results/partB{SFX}"
PARTC = f"results/partC{SFX}"
PLOTS_C = f"plots/partC{SFX}"
PREV_CANDS = {"v3": "results/candidates/candidates_v2_corrgrid.csv", "v4": "results/partC/candidates_v3_prov.csv"}.get(V)
