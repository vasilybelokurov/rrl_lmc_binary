#!/usr/bin/env bash
# Refit v3, part A, stages 4-5 re-run after the GP (smooth-wander) noise model (Level 2 only; Level 1 unchanged).
set -euo pipefail
cd "$(dirname "$0")/.."
PY=.venv/bin/python; LOG=logs/refit_v3_A3.log; EXC=results/candidates/candidates_v2_corrgrid.csv
echo "=== $(date '+%F %T') 4/5 Level 2 real (with GP noise fits)" | tee -a "$LOG"
$PY scripts/level2.py --series results/real/series_v3 --apply-cm results/calib/common_mode_v3.json \
    --band-lag results/calib/band_lag_v3.json --out results/real/stats_v3.parquet --workers 6 >> "$LOG" 2>&1
scripts/autocommit.sh "Refit v3: Level 2 real re-run with whole-cycle band alignment"
echo "=== $(date '+%F %T') 5/5 CHECKPOINT (MACHO and OGLE-only)" | tee -a "$LOG"
for S in macho ogle_only; do
  rm -rf results/validation_v3/series_emp_gp_$S
  $PY scripts/level1_sims.py --sample $S --n-stars 60 --seed 44 --per-class empirical=4 --out results/validation_v3/series_emp_gp_$S \
      --chunk 20 --band-lag results/calib/band_lag_v3.json --noise-table results/real/stats_v3.parquet --exclude "$EXC" >> "$LOG" 2>&1
  echo "--- sample $S" >> "$LOG"
  $PY scripts/compare_noise_real_sim.py --sims results/validation_v3/series_emp_gp_$S --kind empirical >> "$LOG" 2>&1
done
echo "--- timescale test: real (MACHO stars) vs empirical sims" >> "$LOG"
$PY scripts/noise_timescale_test.py --n 400 >> "$LOG" 2>&1
$PY scripts/noise_timescale_test.py --series results/validation_v3/series_emp_gp_macho --kind empirical --cm none --n 240 >> "$LOG" 2>&1
scripts/autocommit.sh "Refit v3: checkpoint re-run"
echo "PART A3 COMPLETE" | tee -a "$LOG"
