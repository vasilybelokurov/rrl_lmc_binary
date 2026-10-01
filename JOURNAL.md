# JOURNAL — Binary fraction of LMC RR Lyrae (light-travel-time effect)

Companion project: `../rrl_lmc_rotation/` (package `lmcrrl`; OGLE-IV × Gaia DR3 catalogue, photometric [Fe/H], local Gaia density).

---

## 2026-09-30 — Session 1: kick-off, infrastructure review, data, literature check

### Goal (from user notes)
Measure the RR Lyrae binary population in the LMC, as a population inference rather than a candidate list:
per-star LTTE time-delay curves from OGLE (+MACHO) photometry → injection/recovery on the real sampling →
hierarchical f_bin(P_orb, M2, e). Side topics: (i) the eclipsing / short-period deficit, (ii) companion light vs local
source density (blends), (iii) the intermediate-age cluster RRL as a test of the Bobrick+24 binary channel, (iv) f_bin vs [Fe/H], radius and environment.

### LTTE scale check (M1 = 0.65 Msun, circular, edge-on; Kepler's third law)
| P_orb [d] | M2 = 0.6 | 0.2 | 0.07 Msun |
|---|---|---|---|
| 300 | 226 s | 98 s | 38 s |
| 1000 | 505 s | 218 s | 85 s |
| 3000 | 1050 s | 453 s | 177 s |
The 3000 d row reproduces the user's table (a1 = 2.10, 0.91, 0.35 AU).

### Reusable from the companion project
- `lmcrrl/io.py`: OGLE-IV parsers and WSDB crossmatch (sqlutilpy + ~/.pgpass; NULL ints → −9999).
- `data/ogle4_lmc_rrl_gaiadr3.parquet`: 41,471 OGLE-IV RRL × Gaia DR3 (RUWE, ipd_frac_multi_peak, n_gaia_within, G−I).
- `data/ogle4_rrab_feh.parquet`: [Fe/H]_I for 27,781 RRab (MV26 relation).
- `data/ogle4_lmc_rrl_density.parquet`: Gaia source count within 30″ (Σ).
- Conventions: JOURNAL.md, tests before conclusions, private GitHub repo, commit after each tested step, local .venv (Dropbox-ignored).

### Data (downloaded to `data/raw/`, Dropbox-ignored; 1.5 GB incl. tarballs)
- OGLE-IV LMC RRL (https://www.astrouw.edu.pl/ogle/ogle4/OCVS/lmc/rrlyr/): phot.tar.gz 188 MB → 41,209 I, 38,602 V light curves.
- OGLE-III LMC RRL (https://www.astrouw.edu.pl/ogle/ogle3/OIII-CVS/lmc/rrlyr/): phot.tar.gz 108 MB (re-issued 2024-04-13) → 24,906 I, 24,439 V.
  **These files already include OGLE-II epochs** (from HJD' = 447, i.e. 1997) for 7,784 stars.
- Not yet downloaded: MACHO (NCI: https://macho.nci.org.au/, TAP + per-tile photometry files), Gaia DR3 epoch photometry (in WSDB: `gaia_dr3.epoch_photometry`).

### Light-curve inventory (`scripts/lc_inventory.py` → `data/lc_inventory.parquet`, 16 s)
- Shared numbering OGLE-III ↔ IV: 24,818 IDs; separation median 0.05″; 8 cases > 1″ (to inspect).
- OGLE-III + IV I band: RRab 17,492, RRc 6,059, RRd 991; median epochs 380 (O-III) + 664 (O-IV); median baseline 5341 d, p90 7052 d (with OGLE-II).
- OGLE-IV only: 16,653 (median 287 epochs; baseline 2231 d).
- OGLE-IV public data end at HJD' ≈ 7507 (2016) for 95% of stars; outer-field additions reach 8924 (2020).
- MACHO IDs (from OGLE-III ident): 8,557; OGLE-II IDs: 7,551.
- HJD' = HJD − 2450000.

### Literature check (subagent; arXiv/Crossref verified)
- Hajdu+2021, ApJ 915, 50 (arXiv:2105.03750): 87 bulge candidates; companion masses in three groups (~0.6, 0.2, 0.067 Msun); none below P_orb ≈ 1000 d; P_orb peak 3000–4000 d; 27,480 RRab.
  **The "~25% completeness" figure in the user's notes is not in the paper.** The closest statement: Blazhko binaries are recovered in about one third of cases.
- Prša+2008, A&A 489, 1209 (arXiv:0808.3560): the HST blend of the OGLE-II "eclipsing RRL" (Soszyński+2003, AcA 53, 93). Soszyński+2016 lists 5 more LMC candidates with P_orb = 1.5–16 d, likely blends.
- Cusano+2021, MNRAS 504, 1 (arXiv:2103.15492): about 7000 overluminous, low-amplitude VMC RRL in the central LMC, attributed to blending.
- Cuevas-Otahola+2025, MNRAS 541, 1434 (arXiv:2411.12741): 23 RRL in 10 intermediate-age clusters (7 LMC); RVs are needed. The Table 3 count needs rechecking (22 vs 23).
- Bobrick+2024, MNRAS 527, 12196 (arXiv:2208.04332): binary-made RRL have P_orb 1000–1800 d and M2 0.65–1.9 Msun; only systems younger than 3 Gyr have M2 > 1 Msun.
- **No systematic LTTE/O−C search of LMC RRL was found.** Only single-star studies exist: Derekas+2004, MNRAS 354, 821; Li+2022, MNRAS 510, 6050 (BE Dor, not LTTE). Still to run: an ADS full-text check.
- Others: Hajdu+2015 MNRAS 449, L113; Prudil+2019 MNRAS 487, L1; Liška+2016 A&A 589, A94 (TU UMa); Kervella+2019 A&A 623, A116 (PM anomaly, ≥7%).

### Proposed staged plan (to be agreed with the user)
0. Infrastructure: package `rrlbin`, local venv, tests, git. Merged O-II/III/IV light curves with per-survey zero points.
1. Per-star template + time-delay model. Build a Fourier template (order K) and fit τ(t) for all harmonics jointly: m(t) = Σ_k A_k cos[2πkν0(t − τ(t)) + φ_k].
   τ(t) is either piecewise per season (the "O−C" data product) or parametric (quadratic period change + Kepler orbit).
   Test: noiseless and noisy injections on the real cadence recover τ to its CRLB.
2. Nuisance models: period change (Ṗ), Blazhko (amplitude/phase modulation, which shows up as non-common harmonic shifts), stochastic phase jitter, and blending (amplitude dilution).
   Test: inject each on the real sampling and measure the false-positive rate of the LTTE detector.
3. Detection statistic plus injection/recovery → p_det(P_orb, M2, e, i | star). Validate against a mock population.
4. Hierarchical inference of f_bin(P_orb, q) (and e), with a null-mock test and a bulge-like (Hajdu+21) mock test.
5. Science splits: [Fe/H], R_LMC, Σ; cluster RRL (Cuevas-Otahola+25); eclipse-deficit prediction; companion light.
6. Extensions: MACHO (baseline to 1992); Gaia DR3/DR4 epoch photometry.

### Open decisions for the user
- Git repo (private GitHub `rrl_lmc_binary`?), local .venv.
- Whether to fetch MACHO now.
- External critique of the plan (ask-codex) before Stage 1.

---

## 2026-09-30 — User decisions; repo; MACHO fetch; Codex plan review; timing estimator + fail-fast test

### Decisions (user)
Public GitHub repo; fetch MACHO; external critique via ask-codex; run the OGLE pipeline while MACHO downloads.
- Repo: https://github.com/vasilybelokurov/rrl_lmc_binary (public). Local `.venv` (Python 3.13, Dropbox-ignored), `pip install -e .`.

### MACHO fetch (`scripts/fetch_macho.py`, running; log in logs/)
- Uses the NCI TAP service, table public.photometry_view, **one query per (field, tile) restricted to `seqn IN (RRL list)`**.
  Only RR Lyrae rows are returned; the per-tile bulk files (about 15 MB each, about 45 GB for 3038 tiles) are NOT downloaded.
  Checked after 417 tiles: 1,251 distinct stars returned, 0 not requested. Projected size about 330 MB.
- VOTable parsing fails (astropy rejects `character(1)` for sideofpier; the ADQL has no CAST) → CSV sync endpoint.
- Sanity check: 3 stars fold on their OGLE periods (phase-binned scatter 0.33–0.45 → 0.22–0.26 mag in MACHO B).
- **MACHO dateobs is MJD (not heliocentric)** → `rrlbin.io.mjd_to_hjd`. Toward the LMC (ecliptic latitude ≈ −85°) the correction is only ±40 s (tested), not the ±8.3 min on the ecliptic I first quoted. Start vs mid-exposure is undocumented → calibrate against OGLE-II.
- TODO (Codex #8, confirmed): validate completeness (requested vs returned IDs per tile; truncation).

### Disk use (2.4 GB): .venv 0.94 GB; OGLE-IV 0.92 GB and OGLE-III 0.57 GB (each incl. its tarball, 0.19 + 0.11 GB); MACHO about 0.33 GB when done.

### Timing estimator (`src/rrlbin/timing.py`; `tests/test_timing.py`, 8 pass)
Joint fit per star/band: shared Fourier template (K = 8), a zero point per survey segment (O2/O3/O4), and per season a delay τ_j plus a mean offset.
The coarse grid over one cycle runs in the first iteration, then Gauss–Newton; errors are Fisher × sqrt(max(χ²_ν, 1)); 4σ clipping; template ↔ τ iterations until converged.
Tests: dm/dφ vs finite difference; amp/phase round trip; season labels; noiseless per-season delays recovered to < 0.01 s;
noiseless 300-s LTTE sinusoid equals the (dm/dt)²-weighted within-season mean to < 3 s; a 0.3P shift found by the grid; 200 noisy realizations give pull std in (0.85, 1.1), mean < 0.05, O3/O4 zero-point offset recovered to 3 mmag.
Bug caught by the stepwise test: the outer loop stopped as soon as the clipping mask stopped changing (0.8 s error); it now also requires max|Δτ| < 1e-7 d.
Note: the within-season delay is a Fisher-weighted mean, so the orbit fit should use the raw epochs (Codex #3), not the season means.

### Fail-fast: timing precision on 1000 random OGLE-III+IV RRab (`scripts/timing_sample.py` → results/timing_sample/rrab_1000.parquet)
| quantity | p10 | p50 | p90 |
|---|---|---|---|
| seasons | 13 | 14 | 19 |
| per-season τ error [s] | 121 | 253 | 510 |
| best-season error [s] | 50 | 108 | 239 |
| χ²_ν of τ_j about a quadratic (constant Ṗ) | 0.81 | 1.94 | 24.4 |
| rms about the quadratic [s] | 143 | 399 | 1473 |
Error vs amplitude: A_I < 0.4: 476 s; 0.4–0.6: 269 s; 0.6–0.8: 159 s; > 0.8: 113 s.
Fraction with χ²_ν(quad) > 3: 0.34; > 10: 0.16 → **one third of RRab have timing structure beyond a constant Ṗ**
(Blazhko, abrupt Ṗ, underestimated errors, or binaries). This population is the main confounder.
Feasibility reading (NOT an injection result): M2 ≈ 0.6 at P ≳ 1000 d (τ ≳ 500 s) is detectable per star in many stars; 0.2 Msun at 3000 d (453 s)
only for higher-amplitude stars (S/N ≈ 453/253 × √7 ≈ 4.7 for the median star, edge-on); 0.07 Msun is statistical only.

### Codex plan review (docs/reviews/2026-09-30_codex_plan_review.md, effort high) — my verification
| point | verdict | action |
|---|---|---|
| Blazhko and Ṗ injections belong in the first end-to-end test | agree | next step |
| a common harmonic delay is not a Blazhko veto (modulation can shift the phase coherently) | plausible | compare with a free per-harmonic phase/amplitude fit |
| season O−C cannot separate long-P orbits from a quadratic | agree | fit the raw epochs jointly; sensitivity vs P_orb |
| survey-wide false-alarm control (max statistic over the search) | agree | full-search null simulations |
| p_det must be conditioned on cadence, errors, shape, Blazhko, Ṗ, crowding, time system | agree | injection on real LCs |
| MACHO fetch completeness not validated | confirmed | validator after the run |
| MJD vs HJD | already known | implement + test |
| precision "tens of s to 2 min" | refuted by the table above (median 253 s) | — |

---

## 2026-09-30 — Injection–recovery on real OGLE cadences (run 1)

### Code
- `timing.py`: each season now also fits an amplitude scale α_j (Blazhko/blend diagnostic). The season epoch is the Fisher-weighted time.
- `ltte.py`: Irwin delay (Kepler solver), a1 sin i / c, mass function; `oc_search` = quadratic + circular sinusoid on a frequency grid
  (300 d … 2 × baseline, oversampling 5), jitter s profiled under H0 and H1; statistic D = 2 max_P ΔlnL.
- `simulate.py`: light curves on real epochs/errors from each star's fitted template: Ṗ, Blazhko (coherent phase + amplitude modulation), abrupt ΔP/P, random-walk phase, LTTE.
- `io.py`: `mjd_to_hjd` (Mount Stromlo; ±40 s toward the LMC), `read_macho`.
- Tests: 14 pass (timing 8, ltte 5, io 1). The LTTE amplitudes reproduce the table; the eccentric peak-to-peak matches 2A·sqrt(1 − e²cos²ω); the O−C search recovers a 600-s, 2500-d signal; the null D has the expected distribution.
- `scripts/inject_recover.py` (300 random RRab, 21 simulations each; about 12 min on 6 cores) → results/inject/run1.parquet;
  `scripts/analyze_inject.py` → results/inject/run1_summary.txt, plots/inject_run1.png.

### Results (per-star FAP from 1200 null sims: D(5%) = 13.5, D(1%) = 17.7, D(0.1%) = 23.9)
Fraction with D > 17.7: null 0.010; jump 0.055; rwalk 0.21; **Blazhko 0.42**; LTTE 0.29 (over the broad injected prior); **real 0.24** (D > 23.9: 0.147).
- Blazhko false positives vs P_B: 12% (20–100 d), 13% (100–300), **65% (300–1000), 84% (1000–3000)**. Long-period modulation mimics LTTE.
- The amplitude diagnostic separates them: median χ²_ν(α_j) = 3.3 / 6.8 for P_B = 300–1000 / 1000–3000 d vs 0.96 for LTTE and nulls.
  With the veto χ²_ν(α) < 2, Blazhko FP falls 0.42 → 0.13, LTTE recovery is essentially unchanged (0.294 → 0.288; 0.68 → 0.68 for M2 > 0.4, P > 1000 d), and **real 0.24 → 0.08**.
- LTTE recovery (FAP 1%, isotropic inclination, M1 = 0.65, half eccentric) vs amp/σ_season: < 1: ≤ 2%; 1–2: 19%; 2–4: 60%; > 4: 80%.
  vs (P_orb, M2): M2 0.4–1.5: 25% (300–1000 d), 70% (1–3 kd), 67% (3–10 kd); M2 0.15–0.4: 4%, 28%, 35%; M2 0.05–0.15: 0–10%.
  Recovered orbits with P < 6000 d: 74% have |P_best/P − 1| < 0.2. The practical threshold is a1 sin i / c ≳ 500–700 s (figure, middle panel).
- Real stars: jitter under H0 p50/p90 = 125/1636 s vs null p90 159 s → **real timing noise has a heavy tail not in the white-noise nulls.**

### Interpretation (tested vs not)
- Sensitivity: roughly Hajdu+21-like companions (≳ 0.4 Msun, 1–10 kd) are recoverable in about 70% of cases; the 0.2-Msun group about 30%; the 0.067-Msun group is not reachable per star.
- **8% of real RRab still pass D and the α veto**, vs 0.3% bulge candidates in Hajdu+21 (87/27,480). The real excess is therefore dominated by
  timing nuisances that the α veto does not catch (phase-only modulation, random-walk phase, abrupt changes), not by binaries.
  The white-noise null is not a valid null for real stars: false-alarm control must use a realistic timing-noise model.
- Next discriminators: (1) Keplerian shape (eccentric fit) and a sinusoid-vs-red-noise model comparison; (2) out-of-sample prediction
  (fit OGLE-III+IV, predict MACHO/OGLE-II seasons: an orbit predicts, a random walk does not); (3) Blazhko side-peaks in the light-curve spectrum; (4) per-harmonic phase coherence.

### Red-noise null and period windows (same run-1 simulations)
- `ltte.oc_search_red`: H0 and H1 with a random-walk phase covariance q·min(t_j, t_k) plus white jitter (grids over s and q); test: reduces to `oc_search` when q = 0 (15 tests pass).
  **It does not separate random walks from orbits**: 40 random walks (800 s rms) give median D 15.8 (white null) vs 15.7 (red null); 40 LTTE (700 s) give 23.7 vs 22.0.
  With about 14 seasonal points, a quadratic plus a sinusoid describes red phase noise as well as a random walk when P ≳ half the baseline.
- Detections (D > 17.7 and χ²_ν(α) < 2) by best-fit period:
  | P_best | null | jump | rwalk | Blazhko | LTTE (all injected) | real |
  |---|---|---|---|---|---|---|
  | 300–1500 d | 0.006 | 0.027 | 0.107 | 0.10 | 0.11 | 0.020 |
  | 1500–2700 d | 0.000 | 0.000 | 0.010 | 0.026 | 0.058 | 0.010 |
  | > 2700 d | 0.002 | 0.027 | 0.083 | 0.004 | 0.12 | 0.050 |
  LTTE with M2 > 0.4 recovered in its own P window: 0.34 (300–1500), 0.53 (1500–2700), 0.58 (> 2700).
- Reading: 1500–2700 d (about 2–3.5 cycles in OGLE-III+IV) is the cleanest window; P ≳ T/2 is degenerate with phase wander in OGLE-III+IV alone.
  The nuisance rates depend on the assumed nuisance priors (arbitrary here), so real false-positive rates are NOT yet known.
- **Decisive test = out-of-sample prediction**: fit an orbit to OGLE-III+IV (2001–2016), then predict the 1992–2000 MACHO/OGLE-II delays.
  An orbit predicts; a random walk / Blazhko does not. This is why MACHO is needed; it also fixes the long-P degeneracy (baseline → 24 yr).

---

## 2026-09-30 — MACHO fetch restart; MACHO–OGLE delay tie (partial: 1,076 stars in fetched tiles)

### MACHO fetch
The first run stalled (4 workers each stuck on a slow tile under a 600-s timeout). Restarted with a 90-s timeout, 2 retries,
per-tile completeness logging (`data/raw/macho/missing_stars.txt`) and `failed_tiles.txt`. At 12:40: 1,370/3,038 tiles, 146 MB, 0 failures,
1 tile with 2 missing stars (14.8975: seqn 585, 597).

### Season labels and delay gauge (`timing.py`)
- MACHO observed the LMC nearly year-round (gaps 30–80 d scattered through the year), so gap-based seasons merged years.
  `year_labels(t)`: season = observing year with the boundary at HJD′ mod 365.25 = 245 (middle of the OGLE seasonal gap, measured from OGLE-III/IV).
- **Gauge fix:** a constant delay is degenerate with the template phase. Delays are now defined so that the template's fundamental harmonic has
  phase 0 (the time of the fundamental's maximum relative to T0). The template phase absorbs the data's mean delay.
  Test: delays are independent of T0 to < 1 s, the fundamental phase is 0, and the delays equal the truth minus a mean (16 tests pass).

### MACHO B vs OGLE I delay offset (`scripts/macho_ogle_offset.py` → results/macho/offset_b_partial.parquet)
RRab with OGLE-II epochs and a MACHO ID, in fetched tiles: 1,021 fitted; overlap years (1997–1999) p10/50/90 = 3/3/4.
- Per-season timing error: **MACHO B 135 s** (median; p10 66, p90 261) vs OGLE I 215 s. MACHO (about 1000 epochs, 1992–1999) is a strong data set, not just a baseline extension.
- Offset Δ = τ_MACHO,B − τ_OGLE,I (781 stars with σ_Δ < 200 s): phase lag of the fundamental −0.0336 cycles (robust SD 0.0081);
  Δ[s] = −5756 P[d] + 1562; robust residual SD 327 s vs median error 127 s → **intrinsic star-to-star scatter in the band lag of about 300 s**; 1% outliers > 1500 s.
- Consequence: for stars without an OGLE-II overlap, the MACHO segment can be tied to OGLE only to about ±300 s (prior on a free per-star offset).
  For stars with overlap (about 2.4k RRab), Δ is measured per star to about 100–150 s. A season-difference χ²_ν median of 1.67 means the overlap differences are mostly consistent.

### Next
When the fetch completes: rerun the calibration on all tiles; add a per-star MACHO offset column (constrained by the overlap years) to the O−C search;
rerun injection–recovery with MACHO+OGLE cadences (24-yr baseline) and test out-of-sample prediction (orbit fitted on 1997–2016 → predict 1992–1996).

### Is the MACHO–OGLE time offset worrying? Control: OGLE V − I (same time system) (`scripts/band_lag_ogle_vi.py` → results/macho/lag_ogle_vi.parquet)
User question: why a free MACHO–OGLE delay offset, and is it worrying? Answer: delays are template-based, and the fiducial (fundamental phase)
falls at a slightly different time in each band (the band phase lag), plus any timestamp convention offset. The worry was that the intercept
of Δ(P) = −5756 P + 1562 s is non-zero, which could mean a 26-min MACHO timestamp error.
| relation | median lag [cycles] | slope [s/d] | intercept [s] | robust scatter |
|---|---|---|---|---|
| OGLE V − OGLE I (N = 336) | −0.0353 (sd 0.0115) | −6497 ± 277 | +2039 ± 151 | 405 s (err 107) |
| MACHO B − OGLE I (N = 781) | −0.0336 (sd 0.0081) | −5663 ± 214 | +1576 ± 115 | 327 s (err 127) |
→ The intercept and the 300–400 s star-to-star scatter are **physical** (a shape/period-dependent band lag), seen equally within OGLE's own time system.
At P = 0.57 d: V−I = −1664 s vs MACHO B − OGLE I = −1652 s → **MACHO timestamps agree with OGLE to about ±50 s** (assuming MACHO blue ≈ V; not exact).
The free offset is measured per star to 100–150 s where MACHO and OGLE-II overlap (1997–1999); elsewhere it is set by the lag relation only to ±300–400 s.

---

## 2026-09-30 — Injection–recovery with MACHO B + OGLE I (`scripts/inject_recover_macho.py` → results/inject/macho_run1.parquet; summary in results/inject/macho_run1_summary.txt)
292 random RRab with OGLE-III+IV and a fetched MACHO tile (central LMC; 3,579 eligible at run time), 21 sims each; a single delay realization on the union of epochs; year seasons;
free MACHO offset column; OGLE-only statistic D_O recorded on the same simulations. Veto: χ²_ν(α) < 2.
| class (FAP 1% per star) | OGLE only | MACHO + OGLE |
|---|---|---|
| null | 0.010 | 0.010 |
| Blazhko (with veto) | 0.158 | 0.145 |
| jump | 0.074 | **0.192** |
| random walk | 0.236 | 0.303 |
| LTTE (broad prior) | 0.370 | **0.476** |
| real (with veto) | 0.164 | **0.247** |
LTTE recovery, MACHO + OGLE (OGLE only): M2 0.4–1.5 → 0.39 (0.32) at 0.3–1 kd, **0.86 (0.77)** at 1–3 kd, **0.91 (0.73)** at 3–10 kd;
M2 0.15–0.4 → 0.15 (0.07), **0.53 (0.44)**, **0.72 (0.46)**; M2 0.05–0.15 → 0.01, 0.19 (0.10), 0.28 (0.20). Period accuracy (|ΔP/P| < 0.2): 0.90 vs 0.71.
Recovery vs amp/σ: 1–2: 0.41 (0.25); 2–4: 0.82 (0.63); > 4: 0.95 (0.81).
Reading:
- MACHO clearly raises sensitivity, especially for long periods and the 0.2-Msun class, and pins P_orb much better.
- The longer baseline also raises false positives from abrupt period changes (7 → 19%) and random walks (24 → 30%), and the fraction of real stars passing (16 → 25%).
  The real passing fraction ≫ any plausible LTTE fraction (Hajdu+21 bulge candidates ≈ 0.3%) → the **real O−C noise is dominated by intrinsic, non-Keplerian period changes**.
- The detection problem is therefore model selection between Keplerian and stochastic/abrupt period changes, and a Δχ² against a quadratic is not enough.
  Candidate discriminators: (1) Keplerian coherence over ≥ 2 cycles, i.e. predictive tests (fit one part of the 24 yr, predict the rest);
  (2) an explicit red-noise/jump alternative with its population fitted to the real stars (hierarchical), not arbitrary sim priors.

---

## 2026-09-30 — Predictive (out-of-sample) orbit test (`ltte.predictive_score`; results/inject/macho_pred_run1.parquet)
Stars with OGLE-II + MACHO (1,536 eligible at run time; 292 used + 8 errors), 21 sims each. Held out: MACHO seasons before 1997 (about 5).
Training: the rest (the MACHO offset is fixed by the 1997–1999 overlap). Score = lnL_pred(test | best training orbit) − lnL_pred(test | quadratic),
with the parameter uncertainty in the predictive covariance. Unit test: toy orbits score higher than random walks (median > 2 and larger by > 2).
- Detections (D > D1 = 16.8 and χ²_ν(α) < 2) with score > 2: LTTE 53%, Blazhko 37%, jump 27%, random walk 18%, null 2/12.
- Rates over all sims (D ∧ veto ∧ score > 2): LTTE 0.274; random walk 0.074; jump 0.074; Blazhko 0.057; null 0.002; **real 0.041**.
- LTTE retained with score > 2 (vs D ∧ veto alone): M2 0.4–1.5 → 0.19/0.45, 0.56/0.90, 0.63/0.97 for P 0.3–1, 1–3, 3–10 kd; M2 0.15–0.4 → 0.05/0.17, 0.32/0.64, 0.44/0.82.
Reading (tested): the predictive score improves the orbit/nuisance ratio only by about 2–3×. Five held-out years carry little of the shape of multi-kd orbits.
**No per-star cut produces a clean binary list at the LMC**: 4% of real RRab pass all cuts vs about 0.3% bulge candidates.
→ Switch to a population mixture fit: per-star summary statistics (D, P_best, α statistic, predictive score, s, rw_rms) modelled as a mixture of
simulated classes (null / Blazhko / jump / random walk / LTTE), with class fractions fitted to the real distribution. This requires realistic
nuisance simulations: first match their parameter distributions to the real H0 noise fits (s, q) from the science run.

---

## 2026-09-30 — MACHO fetch complete; population mixture fit: first mock tests

### MACHO fetch complete
3,038/3,038 tiles (1 retried), 303 MB. Validation: requested 8,557 stars, returned **8,406 (98.2%)**, missing 151, extra (not requested) 0.

### Mixture fit (`src/rrlbin/mixture.py`; `tests/test_mixture.py`; `scripts/mixture_mock.py`)
Features per star: ln(1 + D), log P_best, ln χ²_ν(α), asinh(jitter/σ_season). Class densities: Gaussian KDE per class (null, Blazhko, jump, rwalk, LTTE)
from simulations; class fractions by EM; errors by bootstrap over stars. Unit test: exact recovery with known Gaussian classes.
Mock test on MACHO+OGLE sims, split by STAR into a density half and a mock half (N = 3000 per mock, nuisance mix 0.55/0.15/0.10/0.20):
| KDE bandwidth / features | f_LTTE = 0 | 0.03 | 0.10 |
|---|---|---|---|
| Scott, 4 features | 0.013–0.015 | 0.042–0.045 | 0.112 |
| 0.15 | 0.094 | 0.114 | — |
| 0.25 | 0.031 | 0.057 | — |
| 0.5 | 0.008 | 0.033 | 0.099 |
| 0.8 | 0.011 | 0.042 | 0.112 |
| no jitter feature, 0.5 | 0.001 | 0.027 | 0.085 |
| ln D + α only, 0.5 | 0.000 | 0.005 | 0.046 |
Bootstrap 68% coverage was poor (0–0.4) with Scott because of the bias.
Diagnosis: the bias depends on the density estimate, because only about 290 stars (21 sims each) build the class densities and mock points from other stars
fall where the nuisance KDEs are thin; the broad LTTE class absorbs them (too narrow a bandwidth → strong positive bias).
The spike at zero jitter is not the cause (the same jitter on the data changes nothing).
→ **Systematic of about ±0.01 in f_LTTE from density modelling alone** at the current simulation size. Fixes: cross-validated bandwidth by held-out stars
(`cv_bandwidth`: selects 0.35–0.5 per class) and many more simulated stars (running: 1,500 MACHO+OGLE stars × 10 sims → results/inject/macho_big.parquet).
The LTTE fraction is conditional on the injected orbital prior (log-uniform P 300–10⁴ d, M2 0.05–1.5 Msun, isotropic, half eccentric).

### Run management (2026-09-30 afternoon)
The first full science run was stopped by the session's 2-h background limit with no output: the laptop was down for about 3 h, so wall time elapsed
while the job was paused (not a slowdown; my thread-oversubscription guess was wrong). The 1,500-star simulation was restarted for the same reason.
Fix: `scripts/chunked.py`, a resumable chunked parallel map (part files per chunk, progress line per chunk, single-threaded BLAS in the workers);
jobs are launched with `nohup caffeinate -i` (detached, idle sleep prevented). Profile: 1.2 s per MACHO+OGLE star single-core
(fit_timing 53%, MJD→HJD 28%, O−C search 17%). Full run: about 24 s per 500 stars on 6 workers.

---

## 2026-09-30 — First look at real candidates (8,000 of 17,492 stars processed; `scripts/plot_oc_candidates.py` → plots/oc_candidates_partial.png)
Cuts: D > 40, χ²_ν(α) < 2, amp/σ_season > 3, baseline/P > 1.5 → 43 stars (0.54%); ranked by amp/sqrt(σ² + jitter²).
- Repeated P_best values (4920, 4430 d) are **period-grid quantization** (nearly identical 8,350-d baselines give the same grid; grid spacing about 500 d near 4900 d); sims show the same.
- Two artefacts (07273, 00870: P ≈ 370–400 d, A ≈ 8,600–12,700 s): **cycle-unwrapping failures** aliasing with annual sampling → add a flag (A > P/8).
- About 8 stars show smooth, near-sinusoidal O−C over about 2 cycles, with MACHO (1992–99) and OGLE (2001–16) joining consistently; the phases differ between stars (no common systematic).
  As LTTE (M1 = 0.65, edge-on): 01106 P = 13.6 yr, a1 sin i = 3.1 AU, f(M) = 0.16, M2,min = 0.65, K1 = 6.8 km/s; 11055: 13.5 yr, 2.6 AU, 0.094, 0.50;
  10449: 12.1 yr, 2.4 AU, 0.094, 0.50; 11047: 12.0 yr, 2.2 AU, 0.078, 0.46; 09732: 12.1 yr, 5.0 AU, 0.87, 1.68 (K1 12.4 km/s); 02150: 13.7 yr, 1.9 AU, 0.035, 0.32; 08205: 12.0 yr, 1.0 AU, 0.006, 0.16.
- NOT claimed as detections: P ≈ T/2, where red phase noise mimics about 2 cycles (rwalk FP 30% in sims); the MACHO offset is free (no OGLE-II → no predictive test);
  the 0.5% rate needs population calibration. Confirmation routes: RVs (K1 4–12 km/s over about 12 yr, after removing the pulsation RV), post-2016 OGLE-IV seasons (predicted turn-over), mixture fit with realistic red noise.

---

## 2026-09-30 — Full real-data run; candidates; common mode; the amplitude ceiling examined

### Full run (`scripts/real_oc.py` → results/real/oc_all.parquet, 16 MB; 26 min on 6 workers)
17,492/17,492 RRab OK; with MACHO 6,614; with the predictive test (OGLE-II overlap) 2,378. Median σ_season 226 s; median extra jitter 136 s.

### Candidates (`scripts/plot_oc_candidates.py` → plots/oc_candidates.png, plots/oc_candidates.csv)
Cuts D > 40, χ²_ν(α) < 2, amp/σ > 3, ≥ 1.5 cycles, amp below the LTTE ceiling → **64 stars (0.37%)**.
Notable: 13854 (OGLE-II overlap, predictive score +7.7, P = 4918 d, A = 1281 s); 13392 and 08275 (P ≈ 2885 d, about 3 cycles, A ≈ 350–370 s); the P ≈ 4400–5000 d group (09732, 10449, 01106, 11055, 11047, …).
Several P ≈ 4400 d stars share a phase (09732, 20052, 10900: minima near HJD′ 1000 and 5300) → prompted the common-mode test.

### Common-mode timing test (`scripts/common_mode.py` → plots/common_mode.png)
O−C residuals (after each star's quadratic + MACHO offset) stacked by year over well-behaved stars (jitter < 500 s):
MACHO years 1992–1999: −48, +9, +19, +30, +40, +39, −67, −113 s (robust s.e. about 3 s) → a **significant common MACHO timing pattern of about 100 s**.
OGLE: −21 … +24 s (s.e. 2–5 s). It cannot produce 1,000–2,500-s candidate signals, but it must be subtracted (per-year common-mode correction) before the final fits.
Caveat: part of the pattern may come from the per-star quadratic + offset removal (to be tested on simulations, where no common mode exists).

### The amplitude ceiling (user concern: "what causes it — we need to establish that") (`scripts/ceiling_diagnostics.py` → results/real/ceiling_diag.parquet)
Ceiling = a1/c for M1 = 0.65, M2 = 2 Msun, edge-on, at the fitted period. Real stars above it: 5,309 (30.4%). Dominant cause:
| cause | share | median D | P_best | amp vs ceiling | jitter |
|---|---|---|---|---|---|
| short_P (P_best < 1000 d; noise fit) | 61.8% | 13.2 (below D1) | 368 d | 1095 vs 524 s | 483 s |
| large_pc (O−C range > P/4) | 12.6% | 30.0 | 4022 d | 4994 vs 2579 s | 2500 s |
| other (P_best at the grid edge, about 10⁴ d; trend beyond a quadratic) | 12.3% | 20.1 | 10257 d | 8387 vs 4815 s | 834 s |
| boundary jump at MACHO→OGLE | 7.6% | 25.2 | 2262 d | 3318 vs 1758 s | 1596 s |
| cycle slip (|Δτ| > 0.35 P between seasons) | 5.8% | 20.5 | 2330 d | 8707 vs 1792 s | 5270 s |
Fraction above vs P_best: 0.39 (300–600 d), 0.16, 0.13, 0.20, 0.41 (> 4000 d).
Simulations (same pipeline): null 3.5%, LTTE 4.5% (all with true amplitude below the ceiling, i.e. estimator overshoot), jump 5.3%, rwalk 15%, **Blazhko 30%**; 'real' rows 21%.
**Conclusion:** the ceiling excludes 4.5% of true LTTE (a measured completeness loss) and is otherwise dominated by insignificant short-P noise fits (62%) and genuine large intrinsic
period changes (25%). It is legitimate only as a cut applied identically to sims and data (in the mixture), not as a pre-filter. Cycle slips (6%) are a pipeline issue → robust unwrapping to do.

### Common mode is real; correction applied (`scripts/common_mode.py --null`; `scripts/reanalyse_oc.py` → results/real/oc_all_cm.parquet)
- Null test (white-noise delays at the real epochs, errors and flags; same stacking): all years consistent with 0 (|median| ≲ 4 s, s.e. 2–3 s) →
  the MACHO −113 … +40 s pattern and the OGLE +10 … +24 s (2006–2008) drift are **real common timing systematics**, not fit artefacts.
- Correction: subtract the per-(survey, year) median residual (≥ 100 stars per bin) from every star, then recompute D, P, amplitude, jitter,
  red-noise H0 and the predictive score (4.3 min, 6 workers). Remaining MACHO common mode after one pass: −15 … +4 s. Median per-star correction rms 11.7 s.
- Candidates (same cuts): 64 → 67; 62 in common, 2 lost, 5 new; for the common ones median |ΔP/P| = 0.000, amplitude ratio 1.001, D ratio 0.999 → **candidate list robust to the common mode**.
- The sims contain no common mode by construction → the corrected data and the sims are consistent. Robust unwrapping (cycle slips) is deferred to the next simulation round so that both pipelines stay identical.

### Real vs simulated timing noise; extended nuisance classes; a propagation bug
- Jitter under H0 [s], p50/p75/p90/p97: real MACHO stars 144/564/**1594/3811**; sims: rwalk 169/384/641/840, Blazhko 388/657/978/1250, jump 31/129/217/311, null 0/27/94/193.
  About 15% of real stars are noisier than any simulated nuisance → a mixture fit would push them into the broad LTTE class (the catch-all bias seen in the mock tests).
- New classes in `inject_recover.draw`: `rwalk_big` (rms logU 1000–15000 s) and `jump_big` (|ΔP/P| logU 2e-5–2e-4), matching the large period changes found in the ceiling diagnostics.
  Tail run: 400 MACHO+OGLE stars × (3 + 3) → results/inject/macho_tail.parquet.
- **Bug:** `--per-class` was ignored because macOS multiprocessing spawns workers, which re-import the module defaults. Fixed by passing the dict in each job.
  Consequence: macho_big runs the default 20 sims per star (not 10). The results are valid, just larger and slower.

---

## 2026-09-30 evening — Mixture fit with all classes; completeness-corrected upper limit
Sims: macho_big (1,500 stars × 20; per-class bug → defaults) + macho_tail (400 × rwalk_big 3 + jump_big 3) + earlier runs (≈ 32k sims).
CV bandwidths: 0.25–0.45. Mock recovery (N = 3000; nuisance mix 0.45/0.15/0.08/0.17/0.05/0.10):
f_true 0 → 0.001 ± 0.003; 0.01 → 0.008 ± 0.006; 0.03 → 0.023 ± 0.008; 0.10 → 0.079 ± 0.013 (bias flips sign vs the small-sim run).
Misspecified mock (rwalk > 400 s only, different mix): 0.03 → **0.002**. → **f_LTTE is only weakly identified in these 4 features; it depends on the nuisance shapes at the factor-≥2 level.
Not a measurement yet.**

### Upper limit (no nuisance model needed; `inline analysis, to be scripted`)
Candidate cuts (D > 40, χ²_ν(α) < 2, amp/σ > 3, ≥ 1.5 cycles, amp < ceiling), pass rates in sims (baseline set to 8350 d for the cycle cut):
LTTE 0.132 (broad prior); jump_big 0.072; Blazhko 0.023; rwalk 0.005; jump 0.003; rwalk_big 0.003; null 0.
Completeness (isotropic, half eccentric): M2 0.4–1.5, P 1–10 kd: 0.40; M2 0.15–0.4: 0.11.
Real MACHO+OGLE stars (common-mode corrected): N = 6614, k = 60 (0.91%); 95% Poisson upper limit 74.4.
→ **f(M2 0.4–1.5 Msun, P 1–10 kd) < 2.8% (95%)**; f(M2 0.15–0.4) < 10.7%, assuming all candidates were of that class (conservative).
Main contaminant: large abrupt period changes (jump_big). Next: an explicit break-vs-sinusoid test per candidate; light-curve discriminators; RVs.

---

## 2026-09-30 evening — Summary statistics and candidate sheets (user request); plots/ convention
- **Convention (user):** all plots as PNG in `plots/` only (former figures/ moved; CSVs → results/candidates/; no PDFs).
- `scripts/plot_summary_stats.py` → plots/summary_stats.png, results/real/candidates.csv. Funnel (all / with MACHO): D > 40: 705/504; + α veto: 306/233; + S/N > 3: 299/226;
  + ≥ 1.5 cycles: 92/71; + ceiling: **67/60**. Sim pass fractions: LTTE 0.132, jump_big 0.072, Blazhko 0.023, rwalk 0.005, jump 0.003, rwalk_big 0.003, null 0.
- `scripts/plot_candidate_sheets.py` → plots/candidates/NN_<id>.png (67 sheets): A/B delay-corrected folds (OGLE I, MACHO B) with templates;
  C rising branch of the orbit-extreme OGLE seasons with the period change removed (raw points; the shift equals the fitted orbit signal, e.g. 57 min for 20052);
  D O−C with circular orbit + residuals (jitter-inflated errors); E D(P); F season amplitudes α_j.
  Panel C is a visualisation (it confirms the delay is in the light curve), not a test against red phase noise.

---

## 2026-10-01 — Corrected-grid numbers; write-up §5 "Selection of candidates" (user request)
- Sims v2 (`results/inject/macho_v2.parquet`: 1,472 MACHO+OGLE stars; 2 null, 2 Blazhko, 1 jump, 1 rwalk, 2 jump_big, 2 rwalk_big, 4 LTTE each; Nyquist-safe grid;
  season series saved). `scripts/selection_numbers.py` → results/real/selection_numbers.json (replaces the earlier inline upper-limit calculation).
- Null D thresholds (2,944 nulls): 12.4 (5%), 15.6 (1%), 20.2 (0.1%); max 23.0.
- Funnel real (all / MACHO): 691/496 → 301/232 → 294/225 → 88/71 → **69/62**. Sim fractions after all cuts: LTTE 0.136, jump_big 0.056, Blazhko 0.025, rwalk_big 0.005, rwalk 0.003, jump 0.001, null 0.
- Completeness: M2 0.4–1.5: 0.46 (0.8–3 kd), 0.34 (3–10 kd), 0.42 (1–10 kd); M2 0.15–0.4: 0.07, 0.12, 0.10; M2 < 0.15: ≤ 0.01.
- Upper limit (N = 6614, k = 62, k95 = 76.6): **f < 2.8% (M2 0.4–1.5, 1–10 kd), f < 11% (M2 0.15–0.4)**, unchanged by the grid fix.
- Candidates: P 2.3–13.8 yr (median 10.9; 81% at 2.5–5 kd); 52% above the edge-on 0.5-Msun amplitude; M2,min 0.11–2.0; K1 1.8–16.5 km/s. Sheets regenerated (plots/candidates/, 69).
- Write-up: new §5 (panel-by-panel description of plots/summary_stats.png on a landscape page; Table 2 of cuts with per-class pass fractions; the motivation for each cut); §6 status updated with the corrected numbers. 13 pp.

---

## 2026-10-01 — Pipeline v3: all fixes implemented and validated at small scale (full refit awaits user approval)
User: implement all fixes; test before the refit; **the full refit needs the user's approval**.
Aliasing question (user: "the sampling does not repeat exactly from year to year"): `scripts/alias_test.py` on 200 real MACHO+OGLE stars:
season epochs scatter by 40 d rms within the year; spectral window at 1/yr = 0.85; injected orbits on a grid from 300 d are recovered at the TRUE period
86–100% of the time (alias ≤ 10% at A = 2σ, ≤ 1% at 4σ) → the 800-d cut was over-conservative (my mock had near-degenerate sampling).
Real limit = season-mean smearing, |sinc(π·240/P)| = 0.23 (300 d), 0.50 (400 d), 0.66 (500 d) → **grid from 400 d**.

### New code (tests: 29 pass)
- Level 1 `src/rrlbin/pipeline.py` (load → fit per band → series): **MACHO R** added (band 2); `timing._season_coherence`: per-season delay of the fundamental
  minus that of the higher harmonics (**harmonic coherence**; 0 for a pure time shift). Test: pure delay χ²_ν < 2.5; a shape change gives > 4.
- Level 2 `src/rrlbin/oc.py` (shared by data and sims): **robust unwrapping** (local linear prediction from ±2 neighbours; repairs slips, keeps a 1.5-cycle drift);
  separate B/R offsets with **band-lag priors** (pseudo-rows without jitter); circular and **2-harmonic (eccentric)** orbit fits (the eccentric amplitude is recovered better);
  **alias flag** (D_best − max D at ±1/yr aliases); red-noise null; **predictive test** (the period is searched on the TRAINING seasons only; leakage found and fixed);
  requires each MACHO band tied to OGLE by OGLE-II overlap or a prior (bug found: 08101 had a score without such a tie); pooled per-band α and coherence χ²_ν;
  **iterative common mode**, with the degenerate directions (common quadratic, per-band constants) projected out and a clipped mean
  (the median version did not converge: 19-s per-iteration drift; now a synthetic test converges; real data settle into a 7-s cycle, negligible).
- Scripts: level1_real.py, level1_sims.py (all bands from one delay realization; new class **empirical** = white season jitter + random walk with (s, rw) drawn
  from the real H0 fits), level2.py (identical stats for both), calibrate_band_lag.py.

### Validation (results/validation_v3/: 200 real OGLE-II+MACHO stars; 30 sim stars × 18)
- Real: 198/200 have I + MACHO B + R. Common mode (200 stars): MACHO 1998–99 about −90 … −290 s; OGLE ±10–70 s (noisy with 200 stars).
- Band lag: B −6730 P + 2001 s, intrinsic sd 325 s (lag −0.036 cycles); **R −2920 P + 851 s, sd 161 s (lag −0.016 cycles)**.
- MACHO R per-season errors are comparable to B (e.g. 08101: I 101, B 59, R 81 s; D 56 → 79 with R).
- **Bug found by the null check:** with priors, null D reached 39.5 (60 sims), because simulated MACHO bands had zero offset (gauge-fixed templates),
  in 5σ tension with the lag prior. Fix: sims draw each star's B/R offsets from the lag relation + intrinsic scatter. After the fix: null D p50/p90/max = 7.7/12.5/15.4.
- Coherence: all sim classes about 1.0 (the sim Blazhko is coherent by construction); real p90 1.8 vs sims about 1.4 → real shape changes exist → a useful veto.

### Proposed full refit (awaiting approval)
Order: real L1 (17,492 stars, ~75 min) → common mode → band lag → real L2 (~30 min) → sims L1 with the real noise table + lag (1,500 stars × 18, ~2.5 h) → sims L2 (~45 min).

### Pre-refit review (user: "discuss with codex … don't trust codex blindly") — docs/reviews/2026-10-01_codex_pipeline_v3.md
| Codex finding | Codex severity | my verdict | action |
|---|---|---|---|
| predictive test: full-series unwrapping before the train/test split leaks held-out info | high | confirmed in code (tiny in practice) | fixed: train-only unwrap; held-out block joined by one cycle shift; test that a held-out cycle shift leaves P_train, D_train unchanged |
| band-lag priors calibrated on the same stars | med/high | correct in principle; each star ≈ 1/1500 of the relation; overlap stars dominated by their own data | none (noted) |
| empirical noise drawn from fits that include binaries | medium | by design; conservative | `--exclude candidates.csv` for the noise table |
| common mode may absorb coherent signals | medium | needs phase-coherent orbits across many stars; implausible | `--cm-exclude candidates.csv` anyway |
| unwrapping of clustered slips | medium | mostly addressed by the conservative rule | test: clustered slip never increases jumps |
| stale chunk resumption | medium | confirmed | run manifest (ids + options, digest); mismatch refused |
| silent band-fit failures; real ok without I | medium | confirmed | ok requires I (as in sims); failure reasons recorded and summarized |
Codex verified (I agree): prior pseudo-rows consistent between H0 and H1; no Level-2 row misalignment; the 2-harmonic amplitude definition is self-consistent.

### My own pre-refit checks
- Robust unwrap on 300 real "slip" stars: the original version changed 17% and made the worst jumps worse (p90 0.49 → 0.58 P): these are near-half-cycle ambiguities, not slips.
  → conservative rule (accept only if the adjacent-season jumps decrease): 0% changed, none worse; the synthetic slip test still passes.
- Orbit recovery v3 vs v2 (150 stars × 4 LTTE + 2 null): D > 40 for M2 0.4–1.5, P 1–10 kd: 0.72 (v2 0.65); M2 0.15–0.4: 0.30 (0.24); P 400–1000 d at D > 25: 0.38 (0.19).
  Null D p99/max 18.2/21.2 (v2 15.6/23.0).
- **Bug found by this comparison:** fitted amplitudes 19% low (0.81) → name collision in `oc.design`: the quadratic 'c1' vs the orbit 'c1' made `orbit_amplitude` read the
  linear O−C coefficient as the orbit cosine. Fixed (names q0–q2, sin_h/cos_h) + phase-independence regression test. Now: 0.99 (detections), 0.98 (no D selection).
  Eccentric (e > 0.4): half-ptp/injected 0.87, as expected (half-ptp = a sin i/c · sqrt(1 − e² cos² ω)). D, P and the earlier detection results were unaffected.
- Open (minor): the real-data common mode settles into a 7-s limit cycle (likely clipping-set flips).
Tests: 32 pass.

### Light-curve and timing-noise realism of the simulations (user question)
- `scripts/compare_noise_real_sim.py` (150 stars; real vs null sims of the same stars): median per-season delay error sim/real = 0.98 (I), 1.03 (B), 0.98 (R);
  p10–p90 ≈ 0.8–1.3 → **the photometric noise model reproduces the real timing precision.** The scatter about H0 (χ²_ν, no jitter): real p50/p90 2.0/41 (I),
  2.3/33 (B), 2.1/25 (R) vs nulls 1.1/1.7, 1.0/1.9, 1.3/2.5 → real stars carry extra timing noise (the reason for the empirical class).
- `scripts/noise_origin_test.py`: same-season normalized H0 residuals, real (200 OGLE-II+MACHO stars): ρ(B, R) = 0.58, ρ(MACHO, OGLE) = 0.37 (Spearman);
  nulls: −0.05, +0.03. → the excess is largely intrinsic timing noise common to all bands, plus a MACHO-shared part.
- Empirical class updated: season-jitter variance 64% common to all bands, 36% per instrument group (MACHO B+R shared; OGLE separate); the random walk is common;
  noise drawn RELATIVE to the season error (s/σ, rw/σ from the real table; rescaled by the target star's σ; candidates excluded).
  Check (60 stars × 4): ρ(B, R) = 0.63, ρ(M, O) = 0.37 (real 0.58, 0.37); χ²_ν p50 2.6/2.5/2.3 vs real 1.8/2.4/2.9; **p90 OGLE 65 vs real 18 (tail too heavy)**,
  probably because the noise table still comes from the old pipeline (MACHO B only; large rw absorbing MACHO–OGLE discontinuities).
  → **Acceptance test in the refit:** regenerate the noise table from the new real Level 2, re-run this comparison, and only then use the empirical class.

---

## 2026-10-01 — Refit v3 part A (real data) and the noise-model checkpoint
- User: automated regular commits/pushes → `scripts/autocommit.sh` + `.claude/settings.json` Stop hook (async); orchestrators commit after each stage; chunk parts gitignored.
- User: OGLE-only detections → part B now has two calibration campaigns: MACHO stars (1500) and OGLE-only stars (800), with separate noise tables (`--sample`).
- Part A (scripts/run_refit_v3_A.sh): Level 1 for all 17,492 RRab (21 min): 17,490 ok (2 OGLE I fit failures); bands I only 10,878, I+MB+MR 6,599, I+MB 13.
  Common mode (full sample, candidates excluded): converged (update 95 → 2.7 s). MACHO B and R per-year patterns agree within a few s
  (1992–99 B: −37, +26, +34, +39, +56, +53, −60, −110 s; R: −38, +24, +32, +41, +55, +49, −63, −99) → a timing systematic of the MACHO data. OGLE −66 … +46 s.
  Band lag (full): B −6648 P + 1966 s, intrinsic sd 293 s, −0.036 cycles (N = 1582); R −2656 P + 718 s, sd 218 s, −0.016 cycles (N = 1552).
- **Checkpoint 1 failed** (empirical sims χ²_ν p50 8 vs real 1.8). Cause 1 (**bug**): ~13% of MACHO stars had a band series a whole cycle off the band-lag prior →
  jitter at the grid ceiling and inflated D (median 40 vs 14). Fix: `oc.align_bands` (whole-cycle alignment to the prior; training-only in the predictive test) + regression test.
- **Checkpoint 2 failed** (MACHO sims still ×4). Cause 2: the real excess noise is RED — `scripts/noise_timescale_test.py`: χ²_ν per band p50/p75/p90 = 2.0/4.1/13.7
  vs joint 2.6/10.8/67 (slow wander absorbed within each survey's quadratic), whereas Brownian sims gave per band 2.9/29.5/132.
  Fix: `oc.gp_null` (squared-exponential GP wander, l ∈ {700, 1500, 3000, 6000} d, plus white jitter, **REML**; ML biased A low: injected 900 s → 95 s);
  `simulate.delay_gp`; the empirical class draws (s/σ, A/σ, l) from the real REML fits. Tests (34 pass): recovery A 950 vs 900 s at l = 1500 d; at l = 3000 d poorly constrained (566 vs 900, wide).
  Real fits: A/σ p50/p75/p90 = 2.05/11.0/21.5 (MACHO), 1.46/7.8/21.5 (OGLE-only); l = 700 d most common.
- **Checkpoint 3 (accepted, with documented residual mismatch)**: χ²_ν about H0 per band, real vs empirical sims (60 stars each; p90 rests on ~6 stars):
  MACHO stars I 1.8/18 vs 3.4/33, B 2.4/37 vs 3.6/19, R 2.9/17 vs 2.5/11; timescale test per band 2.0/4.1/13.7 vs 2.2/9.1/36, joint 2.6/10.8/67 vs 3.1/18/107
  (was ×7–10 with Brownian; now within ×1.5–2.5, sims slightly noisier → conservative for false positives); OGLE-only I 2.0/89 vs 2.3/18 (median OK, real tail heavier;
  the jump_big / rwalk_big classes cover the tail). → part B launched.

---

## 2026-10-01 — Part C: investigation of the candidates (plan: docs/partC_plan.md; outputs results/partC/, plots/partC/)
Provisional: Part-A (v3) real data; thresholds and contamination rates to be finalized with Part B. ≤ 2 workers (Part B running).
- **C1** `scripts/partC_candidates.py` (tests: v2 cuts on v2 stats reproduce the frozen 69): **75 provisional v3 candidates** (69 with MACHO).
  v2 → v3: 51 kept (|ΔP/P| median 0.004, A ratio 1.00, D ×1.16 from MACHO R), 18 dropped (12 by the α veto now pooled over MACHO bands → MACHO-era amplitude
  changes, i.e. Blazhko-like; others cycles/ceiling/D), 24 new. Flags: coherence > 2: 0; predictive score available 69, > 0: 31; alias-ambiguous 2; eccentric ΔD > 10: 35.
- **C2** `scripts/partC_crowding.py` (companion Gaia DR3 match): MACHO candidates indistinguishable from MACHO-field parents (KS p: n_Gaia(2″) 0.99,
  Σ 0.87, RUWE 0.53, ΔG−I 0.20; ipd multipeak 0.09); amplitude higher (0.68 vs 0.57 mag, p < 1e-3; selection, opposite to blending dilution).
  The 6 OGLE-only candidates are more crowded (Σ p = 0.04; ipd multipeak median 5.5 vs 0) → caution. 18/75 carry an individual crowding flag.
- **C3** `src/rrlbin/kepler_fit.py`, `scripts/partC_kepler.py` (variable projection + multi-start LS + residual bootstrap; test: e = 0, 0.4, 0.7 injected orbits recovered):
  73/75 fitted (2 with absurd f(M)); e p10/50/90 = 0.17/0.37/0.86; 14 prefer eccentric (ΔBIC > 6); M2,min 0.17/0.52/1.97 Msun; K1 3.0/7.1/19 km/s; χ²_ν 0.63/1.15/1.7.
  6 fits at the e bound (0.95) → flagged unreliable.
- **C4** literature (subagent; docs/reviews/partC_literature.md, 31 refs Crossref-checked): no systematic Magellanic RRL binary search found (ADS full-text query
  could not be run: HTTP 405; query given in the report). Irregular period changes are common in all populations (LMC RRc ~10% strong changes in 6.5 yr, Alcock+2000;
  M3, M5, ω Cen). MACHO epochs = exposure START (Michalska & Pigulski 2005: +150 s) → absorbed in the band offset. Published LMC RRL RVs: 10–35 km/s per epoch
  (FORS/GMOS); confirming K1 ~ 5–10 km/s needs ~1–2 km/s.
- **C5** `scripts/partC_followup.py` (RV = c dτ/dt; test vs the analytic Keplerian RV): 41/73 predicted Δv_sys > 5 km/s over 2027–2030; median V 19.37.
- **C6** `scripts/partC_sheets.py` → plots/partC/sheets/ (75 PNG): O−C with Keplerian + circular curves, residuals, D(P), α_j and coherence by band, flags;
  results/partC/partC_table.csv. Example 13854: B, R and I follow one orbit (P = 12.7 yr, e = 0.2, pred +13.9). One bad MACHO R season (α = 0.45, coherence −6 ks)
  → TODO: season-level outlier rejection in Level 2.
- **Tiers** `scripts/partC_tiers.py` (≥ 1.5 cycles also at the Keplerian period — 11058 had moved to P = 13.2 kd, e = 0.84): **Tier 1: 11, Tier 2: 17, Tier 3: 47.**
  Tier 1 (all MACHO+OGLE, pred > 2, clean vetoes, no crowding flag): 11538, 10449, 09642, 05821, 15158, 16187, 03269, 16750, 16755, 17610, 22630;
  P 2.0–5.3 kd, e 0.15–0.57, M2,min 0.14–0.96 Msun, K1 2.5–11.5 km/s. Note: OGLE-only stars cannot reach Tier 1 (no predictive test).

---

## 2026-10-01 — Our candidates vs the Hajdu+2021 bulge candidates (user request)
- Follow-up status (subagent, docs/reviews/hajdu_followup.md; key quote verified by me on arXiv:2603.28684, Salinas+2026, co-authors Hajdu & Prudil:
  "only a single RRL is confirmed as belonging to a binary system" — TU UMa, a field star): **none of the 87 bulge candidates has been confirmed OR refuted**;
  no published RV follow-up of them. Field-star RV tests: Barnes+2021 (arXiv:2106.05208) 15/19 no binarity; Poretti+2025 rules out KIC 2831097.
  Hajdu+2026 (arXiv:2512.15636) adds 1 candidate, no re-test. A non-refereed note (Zenodo 10.5281/zenodo.23048122) finds ΔBIC favours an orbit on pure noise in 72%.
  → The user's recollection "1–2 confirmed" is not supported: the honest status is "untested".
- Hajdu+2021 Table 1 parsed from the arXiv source (`scripts/parse_hajdu2021.py` → data/external/hajdu2021_binprop.csv: 87 rows; Q1 25, Q2 32, Q3 30).
- `scripts/compare_hajdu.py` → results/partC/hajdu_comparison.csv, plots/partC/hajdu_comparison.png (medians):
  | | N | P_orb [d] | A [s] | e | M2,min | K1 | σA/A | σP/P | σe | cycles |
  |---|---|---|---|---|---|---|---|---|---|---|
  | Hajdu Q1 | 25 | 3437 | 779 | 0.27 | 0.44 | 6.5 | 0.017 | 0.005 | 0.024 | 1.9* |
  | Hajdu Q2 | 32 | 3477 | 355 | 0.31 | 0.14 | 2.6 | 0.037 | 0.015 | 0.054 | 1.9* |
  | Hajdu Q3 | 30 | 4910 | 534 | 0.35 | 0.20 | 3.4 | 0.092 | 0.053 | 0.065 | 1.4* |
  | LMC Tier 1 | 11 | 2914 | 804 | 0.34 | 0.49 | 6.2 | 0.067 | 0.011 | 0.102 | 2.9 |
  | LMC Tier 2 | 17 | 3131 | 736 | 0.32 | 0.44 | 6.2 | 0.070 | 0.016 | 0.127 | 2.4 |
  (*Hajdu baseline assumed 6700 d.) Our candidates occupy the high-amplitude end (A ≳ 350 s; sensitivity), similar to Hajdu Q1 in A, M2,min and K1.

### Gaia photometry as a test (user question) — `scripts/gaia_epochs_check.py`, `gaia_timing.py`, `gaia_timing_control.py`, `gaia_dr4_forecast.py`
- DR3 epoch G for 62/75 candidates: median 34 usable transits (p10–p90 27–41), HJD′ 6864–7875 (2014.6–2017.4), per-transit flux S/N ≈ 59, median G 19.33.
- Test: G template (K = 4) + one delay per half of the DR3 window (split at 7470); the G–I lag cancels in Δτ; compared with each Keplerian prediction.
- Control (228 quiet non-candidates, G 18.8–19.8): Gaia Δτ vs OGLE's own Δτ at the same epochs. With Gaia errors only: robust sd(z) 2.7 (→ the OGLE reference error
  was missing); with OGLE errors included: **robust sd 1.64**, |z| > 5 2.2% (Gaia-own template); a fixed OGLE-I template is worse (1.82). → Gaia errors inflated by 1.64.
- DR3 result: 4/53 informative (|pred| > 2σ): 16755 (T1) z = 0.8 and 13469 z = 0.3 agree; **13392 (T2) z = −6.7 and 01106 (T3) z = +3.5 disagree** (~0.5 expected by chance).
- DR4 forecast (assumed ~2× DR3 epochs; bins 2014.6–2017.4 / 2017.4–2020.1, the latter entirely after public OGLE-IV): testable > 2σ for 6/9 Tier 1 (5 > 3σ),
  8/13 Tier 2, 13/31 Tier 3.
- Write-up: new §7 "Part C: the candidates compared with the Galactic-bulge sample" (tiers, Hajdu follow-up status, comparison table + figure, strengths and weaknesses,
  expectation, Gaia test + DR4 forecast); 16 pp. Note: one DOI (Alcock+2000) was first written from memory; now verified via Crossref.
