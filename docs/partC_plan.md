# Part C — investigation of the candidates

Goal: for each LTTE candidate decide how plausible a binary it is, using information that is independent of the O−C
detection statistic, and prepare the best ones for follow-up. Runs on the Part-A (v3) real-data results; thresholds become
final once Part B is complete (provisional until then). Light on CPU (≤ 2 workers) while Part B runs.

| step | what | method | test | output |
|---|---|---|---|---|
| C1 | provisional v3 candidate list + comparison with v2 | v2 cuts (D > 40, α veto, A/σ > 3, ≥ 1.5 cycles, ceiling) on stats_v3; new flags: harmonic coherence χ²_ν, pooled α (incl. MACHO), predictive score, alias ΔD, eccentric ΔD | unit test of the cut function on a toy table; consistency: v2 cuts on v2 stats reproduce the 69 v2 candidates | results/partC/candidates_v3_prov.csv, v2_vs_v3.csv |
| C2 | crowding / blending | companion-project Gaia DR3 match: neighbours within 2″, RUWE, ipd_frac_multi_peak, G−I offset, local source density Σ (30″); compare candidates with the parent RRab population (percentiles, KS test, rank of each candidate) | KS on two random parent subsamples gives p ≈ U(0,1) (no false signal) | results/partC/crowding.csv, plots/partC/crowding.png |
| C3 | Keplerian orbit fits | full Irwin LTTE (P, A = a1 sin i/c, e, ω, T_p) + quadratic + MACHO band offsets (band-lag priors) + white jitter; multi-start nonlinear least squares; uncertainties by residual bootstrap; ΔBIC vs circular | recover injected eccentric orbits (e = 0.0, 0.4, 0.7) on a real cadence within the bootstrap errors | results/partC/kepler.csv (P, A, e, ω, f(M), M2,min, K1, errors) |
| C4 | literature | ADS / arXiv full-text search for earlier LMC RR Lyrae O−C / LTTE / binary work and for any of the candidate IDs | – (sources cited with links; unverified items labelled) | docs/reviews/partC_literature.md |
| C5 | follow-up feasibility | V magnitude, K1 (with √(1−e²)), predicted systemic-velocity change over 2027–2030 from the Keplerian ephemeris; required RV precision after removing the pulsation RV (RRab ~50–70 km/s) | ephemeris prediction reproduces the fitted O−C at the observed epochs | results/partC/followup.csv |
| C6 | candidate sheets v3 | as the v2 sheets + Keplerian curve, crowding and veto annotations | visual check | plots/partC/sheets/NN_<id>.png (PNG only) |

Acceptance for "strong candidate" (to be finalized with Part B rates): passes the C1 cuts, coherence χ²_ν < 2, no crowding
flag (C2), Keplerian fit consistent across MACHO and OGLE (C3), and, where available, predictive score > 0.
