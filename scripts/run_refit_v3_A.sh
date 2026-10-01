#!/usr/bin/env bash
# Refit v3, part A (real data + noise-table acceptance test). Stops at the checkpoint; part B (full simulations) is launched
# only after the checkpoint has been reviewed. Each stage is logged and auto-committed. Resumable (chunked, manifests).
set -euo pipefail
cd "$(dirname "$0")/.."
PY=.venv/bin/python
LOG=logs/refit_v3_A.log
EXC=results/candidates/candidates_v2_corrgrid.csv
stage() { echo "=== $(date '+%F %T') $1" | tee -a "$LOG"; }
done_() { echo "--- $(date '+%F %T') done: $1" | tee -a "$LOG"; scripts/autocommit.sh "Refit v3: $1"; }

stage "1/5 Level 1 real (all RRab)"
$PY scripts/level1_real.py --out results/real/series_v3 --workers 6 >> "$LOG" 2>&1
done_ "Level 1 real"

stage "2/5 common mode (candidates excluded)"
$PY scripts/level2.py --series results/real/series_v3 --common-mode results/calib/common_mode_v3.json \
    --cm-exclude "$EXC" --no-red --out results/real/stats_v3_cmpass.parquet --workers 6 >> "$LOG" 2>&1
done_ "common mode"

stage "3/5 band-lag calibration"
$PY scripts/calibrate_band_lag.py --series results/real/series_v3 --apply-cm results/calib/common_mode_v3.json \
    --out results/calib/band_lag_v3.json >> "$LOG" 2>&1
done_ "band lag"

stage "4/5 Level 2 real (common mode + band-lag priors)"
$PY scripts/level2.py --series results/real/series_v3 --apply-cm results/calib/common_mode_v3.json \
    --band-lag results/calib/band_lag_v3.json --out results/real/stats_v3.parquet --workers 6 >> "$LOG" 2>&1
done_ "Level 2 real"

stage "5/5 CHECKPOINT: empirical-noise acceptance test (60 stars x 4, new noise table)"
rm -rf results/validation_v3/series_emp_v3
$PY scripts/level1_sims.py --n-stars 60 --seed 44 --per-class empirical=4 --out results/validation_v3/series_emp_v3 \
    --chunk 20 --band-lag results/calib/band_lag_v3.json --noise-table results/real/stats_v3.parquet --exclude "$EXC" >> "$LOG" 2>&1
$PY scripts/compare_noise_real_sim.py --sims results/validation_v3/series_emp_v3 --kind empirical >> "$LOG" 2>&1
done_ "checkpoint (noise acceptance test) - awaiting review"
echo "PART A COMPLETE" | tee -a "$LOG"
