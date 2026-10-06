# PLAN — LTTE binaries among LMC RR Lyrae (living document; update the status column as work proceeds)

Last updated: 2026-10-07. Binary-fraction write-up: docs/binary_fraction/binary_fraction.pdf (f < 5.3%, M2 0.4-1.5 Msun, 1-10 kd). Details and numbers: JOURNAL.md. Design rationale: docs/search_design.md + reviews/2026-10-06_codex_search_design.md.

## Goal
1. Population: the binary fraction f_bin of LMC RRab vs companion mass and orbital period (or an upper limit), with a stated detection efficiency.
2. Individuals: a ranked list of stars with a coherent, LTTE-consistent timing signal, for RV / Gaia follow-up.
   Timing alone can claim "strictly periodic, achromatic, shape-invariant delay consistent with an orbit", not "binary".

## Where we are (2026-10-06)
- Data: MACHO B/R 1992–99, OGLE-II/III 1997–2009, OGLE-IV 2010–2026 (OGLE-team files, ~/data/ogle/, symlinked in data/raw/). 17,490 RRab.
- Pipeline v4 (timing fit with gross-outlier clip, repeated grid, season stability check) run on all stars, all data:
  results/real/series_v4, stats_v4.parquet; common mode, band lag v4. 140 stars pass the old (ad hoc) cuts.
- Frozen-prediction test (blind; tag predictions-2026-10-01): 27 testable → 2 pass (13854 strong, 15158 weaker), 20 fail, 5 inconclusive.
- Old population limit (v3): f < 2.7% (M2 0.4–1.5 Msun, P 1–10 kd), to be replaced by stage 4.
- Decision (user, 2026-10-06): replace the ad hoc cuts and time splits by hypothesis testing with explicit alternatives.

## Stages

| # | stage | output | gate / acceptance test | status |
|---|---|---|---|---|
| 0 | v4 simulations (B1–B4: MACHO 1500 × 18, OGLE-only 800 × 18, 1992–2026 epochs) | results/inject/stats_v4_*.parquet | runs complete, 0 failures | done: B1–B3 (MACHO L1+L2, OGLE-only L1); B4 skipped; old-cut check: real 93 vs ~112 expected from red noise (no excess) |
| 0b | Keplerian fits, all data, 146 candidates | results/partC_v4/kepler.csv | fitted / failed counts | done: 143/146 fitted (e p50 0.39, M2,min p50 0.53 Msun) |
| 1 | **Identifiability**: circular / eccentric LTTE vs quasi-periodic intrinsic modulation (coherence time l from 0.5 to ∞ × P), vs period jumps, vs Blazhko phase-only, on real cadences (22–30 seasons, real errors) | table: Λ distributions, power vs coherence time, cycles, e | defines the coherence bound of H_QP and the period range (≥ 2 cycles) where claims are possible | **done** (1 + 1b): unbounded QP → circular orbits unidentifiable; 2.4% of stars have periodic timing components; coherent intrinsic modulation EXISTS (16 coherent stars above the orbit ceiling) and ~90% of it shows amplitude modulation → discriminant = H_BL (joint τ, α, coherence) + physical ceiling + Keplerian shape, not a coherence bound |
| 2 | **Explicit models** (priority: H_BL joint τ + α + harmonic coherence with season covariances; physical amplitude ceiling as part of H_LTTE) (shared code for data and sims): H_RN (SE GP), H_RN2 (Matérn-3/2), H_J (1–2 ΔP breaks, continuous O−C) + RN, H_BL (phase + amplitude + shape at P_B; joint season covariance of τ, α, coherence), H_QP (bounded coherence), H_LTTE (Keplerian + RN); systematics (common-mode uncertainty, MACHO offsets) as nuisance in all; combinations allowed | src/rrlbin/models.py (+ tests: each model recovers its own simulations; nesting; likelihoods vs brute force) | unit tests pass; per-model parameter recovery on sims | to do |
| 3 | **Per-star evidence**: Λ_k = ln L(H_LTTE) − ln L(H_k) for each alternative; p_k by parametric bootstrap with the full fit (period search included); star p = max_k p_k (intersection-union) | per-star table (all stars, cheap Λ; bootstrap for the top ranks) | p-values uniform under each null across the hyperparameter range (calibration plot); used to RANK, not for FDR over 17k stars (p ~ 3e-6 unreachable) | to do |
| 4 | **Population**: mixture likelihood E[pass] = Σ_i [f ε_i + (1 − f) α_i], or the full mixture over Λ; ε_i (efficiency) and α_i (false pass) per star from end-to-end light-curve simulations with the noise population learned by cross-fitting (candidates excluded) | f_bin(M2, P) or upper limits | **f = 0 mock returns 0**; low-f mocks recovered; robust to the nuisance mix | to do |
| 5 | Validation & bias control: procedure frozen on sims + a random half of real stars (development), run once on the other half (confirmation); light-curve-level check of the delay-level bootstrap; frozen-prediction test kept separate | validation report | all gates passed | to do |
| 6 | Individuals: ranked list with orbits, f(M), M2,min, K1, crowding, band/harmonic consistency; follow-up forecast (RV, Gaia DR4) | candidate table + sheets | — | to do |
| 7 | Write-up update (docs/writeup) and data request / follow-up notes | paper draft | — | to do |

## Main false positives (what the alternatives must cover)
1. Red (slow, random) period wander — dominant; common to all bands; timescales 2–16 yr.
2. Abrupt period changes (ΔP/P ~ 1e-5–2e-4).
3. Long-period Blazhko modulation (P_B 1–8 yr); phase-only variants escape the amplitude veto.
4. Strictly/quasi-periodic intrinsic phase modulation — degenerate with circular orbits; timing alone cannot exclude it.
5. Instrumental: survey joins (MACHO–OGLE ±300 s), yearly common mode, re-reduced fields, cycle counting, wrong season minima.
6. Crowding / blends. 7. Sampling: annual aliases, P > baseline/2.

## Rules we keep
- No conclusions without tests; every new function gets a test; JOURNAL.md entry after each step; commit + push.
- The same code for real data and simulations; selection inside the simulations.
- Frozen products (tag predictions-2026-10-01, v3 results) are never overwritten (RRL_VERSION, *_v4 outputs).
- Codex reviews at design points; every finding verified before acting.
