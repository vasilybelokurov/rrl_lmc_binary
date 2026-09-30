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
- **MACHO dateobs is MJD (not heliocentric)**: an HJD conversion is required (up to ±8.3 min), to be implemented and tested.
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
