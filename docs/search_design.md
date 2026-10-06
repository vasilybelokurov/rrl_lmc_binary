# Search design v5: hypothesis testing for LTTE binaries among LMC RR Lyrae (DRAFT for review, 2026-10-06)

Status: proposal. Replaces the ad hoc chain (D > 40, alpha veto, S/N, cycles, ceiling, single time split) used in v2-v4.
Context and all numbers: JOURNAL.md (2026-09-30 to 2026-10-06).

## 0. Data per star
Season delays y_j (OGLE I 1997-2026, ~22 seasons; MACHO B and R 1992-1999, ~8 seasons each for 6,600 stars), errors sigma_j,
band label, plus per-season amplitude scale alpha_j and harmonic coherence c_j (delay of fundamental minus delay of higher
harmonics). Derived from the light curves by fit_timing (template per band, season delays). 17,490 RRab.

## 1. Hypotheses (explicit alternatives, not just "orbit vs noise")
Common to all: trend (quadratic: P, dP/dt), band offsets with band-lag priors, white jitter s, measurement errors.

| label | timing model beyond the common part | auxiliary prediction (alpha_j, c_j) |
|---|---|---|
| H_RN  | smooth red noise: GP, squared-exponential kernel (A, l) | alpha, c constant |
| H_RN2 | rough red noise: GP Matern-3/2 (A, l) | alpha, c constant |
| H_J   | abrupt period change(s): 1-2 breaks in dP (piecewise-quadratic O-C) + H_RN | alpha, c constant |
| H_QP  | intrinsic quasi-periodic modulation: GP quasi-periodic kernel (A, P_q, l_decay, gamma) + white | alpha, c constant |
| H_BL  | Blazhko-like modulation: phase AND amplitude/shape modulation at a common period P_B: delay ~ sinusoid(P_B), alpha_j and c_j modulated at the same P_B with free amplitudes | alpha, c vary coherently with the delay |
| H_LTTE | Keplerian orbit (P, a sin i / c, e, omega, t_p) + H_RN (noise refitted) | alpha, c constant; same delay in all bands and harmonics |

H_LTTE nests H_RN. H_QP with l_decay -> infinity becomes strictly periodic but sinusoidal-shaped; the LTTE differs by
(i) Keplerian harmonic structure for e > 0, (ii) strict coherence over all cycles, (iii) physical amplitude-period-mass
relation (M2 below a ceiling), (iv) no amplitude/shape change.

## 2. Per-star decision: intersection-union test
For each alternative k in {RN, RN2, J, QP, BL}: Lambda_k = max ln L(H_LTTE) - max ln L(H_k) (all hyperparameters fitted by
ML/REML under each hypothesis; joint likelihood of delays and, for H_BL, alpha_j and c_j).
p_k = parametric-bootstrap p-value: simulate the star's delays from its fitted H_k (same epochs, errors, bands), refit H_LTTE
and H_k, recompute Lambda_k (includes the period search = look-elsewhere and the hyperparameter fitting).
Star-level p = max_k p_k (reject the composite null only if EVERY alternative is rejected; valid for a union null).
Cost control: stage 1 cheap Lambda for all stars; stage 2 bootstrap (N ~ 10^3-10^4) only for stars above a stage-1 threshold,
with the same two-stage rule applied inside the bootstrap.

## 3. Multiple testing and population
- Candidate list: Benjamini-Hochberg on the star-level p, FDR q = 0.05 (and 0.1 reported).
- Population: p-value histogram of all stars; Storey's pi_0 estimate -> number of stars inconsistent with every alternative.
- Detection efficiency: inject LTTE (P 400-20000 d, M2 0.05-2 Msun, isotropic, e distribution) into simulated stars whose noise
  is drawn from the fitted alternatives of real stars (empirical noise population) -> eps(P, M2). f_bin (or upper limit)
  from the number of non-null stars / eps, with the false-discovery share subtracted.

## 4. Simulation level and realism
Bootstrap at the delay level (fast). Validated on a subset at the light-curve level (simulate light curves at real epochs with
the star's template, inject the same delay series, run fit_timing): checks that timing extraction (clipping, unstable
seasons, cycle counting, outliers) does not change the null distribution of Lambda.

## 5. Bias control
- Freeze the procedure (models, priors, thresholds) on simulations + a random half of the real stars (development set);
  run once on the other half (confirmation set); report both.
- The 27 frozen predictions (tag predictions-2026-10-01) remain a separate, genuinely blind test.
- Diagnostics reported, not used as cuts: crowding (Gaia), band and harmonic consistency, number of orbital cycles.

## 6. Open questions
1. Is the intersection-union (max p) too conservative? Alternative: Bayesian model comparison with priors on each H_k and
   a population-level mixture (hierarchical), which yields f_bin directly.
2. Is H_QP identifiable against H_LTTE at all with 2-4 cycles? If not, what can the timing data claim (and what needs RVs)?
3. Can a per-star bootstrap be trusted when the alternative's hyperparameters are poorly constrained (few seasons)?
   Use the population distribution of hyperparameters (empirical Bayes) instead of per-star point estimates?
4. Should the MACHO-OGLE offset (+-300 s where no OGLE-II overlap) and the per-year common mode be nuisance parameters in
   every model (they are now fixed / prior-constrained)?
5. Cost: 17,490 stars x 6 models; stage-2 bootstrap for perhaps 1-3k stars x 6 alternatives x 10^3 sims.
