#!/usr/bin/env bash
# Refit v3, part B (simulations), launched after the part-A checkpoint has been reviewed.
# Two campaigns with separate calibration: stars with MACHO (1500) and OGLE-only stars (800). Auto-commit per stage.
set -euo pipefail
cd "$(dirname "$0")/.."
PY=.venv/bin/python
LOG=logs/refit_v3_B.log
EXC=results/candidates/candidates_v2_corrgrid.csv
LAG=results/calib/band_lag_v3.json
NT=results/real/stats_v3.parquet
PC="null=2,blazhko=2,jump=1,rwalk=1,jump_big=2,rwalk_big=2,empirical=4,ltte=4"
stage() { echo "=== $(date '+%F %T') $1" | tee -a "$LOG"; }
done_() { echo "--- $(date '+%F %T') done: $1" | tee -a "$LOG"; scripts/autocommit.sh "Refit v3: $1"; }

stage "1/4 Level 1 sims, MACHO+OGLE stars (1500 x 18)"
$PY scripts/level1_sims.py --sample macho --n-stars 1500 --seed 21 --per-class "$PC" --band-lag "$LAG" \
    --noise-table "$NT" --exclude "$EXC" --out results/inject/series_v3_macho --workers 6 >> "$LOG" 2>&1
done_ "Level 1 sims (MACHO)"

stage "2/4 Level 2 sims (MACHO)"
$PY scripts/level2.py --series results/inject/series_v3_macho --band-lag "$LAG" --out results/inject/stats_v3_macho.parquet \
    --workers 6 >> "$LOG" 2>&1
done_ "Level 2 sims (MACHO)"

stage "3/4 Level 1 sims, OGLE-only stars (800 x 18)"
$PY scripts/level1_sims.py --sample ogle_only --n-stars 800 --seed 22 --per-class "$PC" --band-lag "$LAG" \
    --noise-table "$NT" --exclude "$EXC" --out results/inject/series_v3_ogle --workers 6 >> "$LOG" 2>&1
done_ "Level 1 sims (OGLE-only)"

stage "4/4 Level 2 sims (OGLE-only)"
$PY scripts/level2.py --series results/inject/series_v3_ogle --out results/inject/stats_v3_ogle.parquet --workers 6 >> "$LOG" 2>&1
done_ "Level 2 sims (OGLE-only)"
echo "PART B COMPLETE" | tee -a "$LOG"
