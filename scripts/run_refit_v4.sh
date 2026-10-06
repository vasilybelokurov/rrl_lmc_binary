#!/usr/bin/env bash
# Refit v4: 1992-2026 baseline (MACHO B/R 1992-1999, OGLE-II/III 1997-2009, OGLE-IV 2010-2026 from the OGLE-team files of
# 2026-10-06) with the v4 timing fit (8-sigma first-iteration clip, grid search in the first 3 iterations, season stability
# check; src/rrlbin/timing.py). Same stages as v3 (run_refit_v3_A*.sh, run_refit_v3_B.sh); the noise-model checkpoint is
# logged and the run continues (user: do the full refit; the checkpoint is reviewed afterwards). Auto-commit per stage.
# The frozen prediction test (tag predictions-2026-10-01, scripts/test_frozen_predictions.py) is NOT touched.
set -euo pipefail
cd "$(dirname "$0")/.."
PY=.venv/bin/python
LOG=logs/refit_v4.log
EXC=results/real/candidates_v3_summary.csv          # v3 candidates: excluded from common mode and noise table
CM=results/calib/common_mode_v4.json
LAG=results/calib/band_lag_v4.json
NT=results/real/stats_v4.parquet
PC="null=2,blazhko=2,jump=1,rwalk=1,jump_big=2,rwalk_big=2,empirical=4,ltte=4"
stage() { echo "=== $(date '+%F %T') $1" | tee -a "$LOG"; }
done_() { echo "--- $(date '+%F %T') done: $1" | tee -a "$LOG"; scripts/autocommit.sh "Refit v4: $1"; }

stage "A1 Level 1 real (all RRab, 1992-2026)"
$PY scripts/level1_real.py --ogle4 extended --out results/real/series_v4 --workers 6 >> "$LOG" 2>&1
done_ "Level 1 real"

stage "A2 common mode (v3 candidates excluded)"
$PY scripts/level2.py --series results/real/series_v4 --common-mode "$CM" --cm-exclude "$EXC" --no-red \
    --out results/real/stats_v4_cmpass.parquet --workers 6 >> "$LOG" 2>&1
done_ "common mode"

stage "A3 band lag"
$PY scripts/calibrate_band_lag.py --series results/real/series_v4 --apply-cm "$CM" --out "$LAG" >> "$LOG" 2>&1
done_ "band lag"

stage "A4 Level 2 real"
$PY scripts/level2.py --series results/real/series_v4 --apply-cm "$CM" --band-lag "$LAG" --out "$NT" --workers 6 >> "$LOG" 2>&1
done_ "Level 2 real"

stage "A5 checkpoint: empirical-noise acceptance (60 stars x 4; MACHO and OGLE-only)"
for S in macho ogle_only; do
  rm -rf results/validation_v4/series_emp_$S
  $PY scripts/level1_sims.py --ogle4 extended --sample $S --n-stars 60 --seed 44 --per-class empirical=4 \
      --out results/validation_v4/series_emp_$S --chunk 20 --band-lag "$LAG" --noise-table "$NT" --exclude "$EXC" >> "$LOG" 2>&1
  echo "--- sample $S" >> "$LOG"
  $PY scripts/compare_noise_real_sim.py --ogle4 extended --sims results/validation_v4/series_emp_$S --kind empirical \
      --out results/validation_v4/noise_real_vs_sim_$S.csv >> "$LOG" 2>&1
done
$PY scripts/noise_timescale_test.py --series results/real/series_v4 --cm "$CM" --n 400 >> "$LOG" 2>&1
$PY scripts/noise_timescale_test.py --series results/validation_v4/series_emp_macho --kind empirical --cm none --n 240 >> "$LOG" 2>&1
done_ "checkpoint (noise acceptance test)"

stage "B1 Level 1 sims, MACHO+OGLE stars (1500 x 18)"
$PY scripts/level1_sims.py --ogle4 extended --sample macho --n-stars 1500 --seed 21 --per-class "$PC" --band-lag "$LAG" \
    --noise-table "$NT" --exclude "$EXC" --out results/inject/series_v4_macho --workers 6 >> "$LOG" 2>&1
done_ "Level 1 sims (MACHO)"

stage "B2 Level 2 sims (MACHO)"
$PY scripts/level2.py --series results/inject/series_v4_macho --band-lag "$LAG" --out results/inject/stats_v4_macho.parquet \
    --workers 6 >> "$LOG" 2>&1
done_ "Level 2 sims (MACHO)"

stage "B3 Level 1 sims, OGLE-only stars (800 x 18)"
$PY scripts/level1_sims.py --ogle4 extended --sample ogle_only --n-stars 800 --seed 22 --per-class "$PC" --band-lag "$LAG" \
    --noise-table "$NT" --exclude "$EXC" --out results/inject/series_v4_ogle --workers 6 >> "$LOG" 2>&1
done_ "Level 1 sims (OGLE-only)"

stage "B4 Level 2 sims (OGLE-only)"
$PY scripts/level2.py --series results/inject/series_v4_ogle --out results/inject/stats_v4_ogle.parquet --workers 6 >> "$LOG" 2>&1
done_ "Level 2 sims (OGLE-only)"
echo "REFIT V4 COMPLETE" | tee -a "$LOG"
