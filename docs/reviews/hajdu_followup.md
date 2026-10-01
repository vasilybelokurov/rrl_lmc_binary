# Follow-up status of the Hajdu et al. 2021 bulge RR Lyrae LTTE binary candidates

Compiled 2026-10-01. Sources opened: arXiv PDFs and abstracts (arXiv API), OpenAlex citation list (27 citing works, checked one by one). ADS was not reachable (bot wall).

## Executive summary

- **The claim that "1 or 2 were confirmed" is REFUTED.** I found no primary source reporting confirmation of any of the 87 candidates by RV, Gaia astrometry, or longer-baseline O-C. I also found none that refutes any of them.
- Papers from 2026, two co-authored by Hajdu, still say only **one** RR Lyrae is a confirmed binary, and that one is TU UMa, a field star. The Gaia studies say no genuine RRL binary has been robustly identified.
- **Radial velocities:** no RV follow-up of the bulge LTTE candidates has been published. The modern RV searches looked at field stars and found nothing.
- **Longer baselines:** the only update is from the same group. Hajdu et al. 2026 used OGLE, KMTNet, MOA, EROS-2 and MACHO data and added one candidate (OGLE-BLG-RRLYR-32226), making 88. They did not re-test the 87.
- **Gaia DR3:** no genuine RR Lyrae has an astrometric orbital solution. This was not a targeted test of the 87.

## (a) Radial-velocity follow-up

- None of the bulge candidates has published RVs. A 2018 Carnegie project summary called RV follow-up "ongoing" ([Carnegie](https://carnegiescience.edu/continued-search-rr-lyrae-binary-systems-towards-galactic-bulge)), but I found no result from it.
- **Barnes, Guggenberger & Kolenberg 2021** ([arXiv:2106.05208](https://arxiv.org/abs/2106.05208)) measured 19 *field* RRL.
  - 15 show no evidence of binarity.
  - 3 are suspicious (SS Leo, ST Leo, AO Peg).
  - For TU UMa, the RVs alone leave the binary detection uncertain.
- **Poretti et al. 2025** ([A&A, DOI 10.1051/0004-6361/202556279](https://doi.org/10.1051/0004-6361/202556279)) tested the Kepler LTTE candidate KIC 2831097. RVs taken at the predicted quadratures were identical (−203 km/s), so the LTTE is ruled out.
- **Salinas et al. 2026** ([arXiv:2603.28684](https://arxiv.org/abs/2603.28684)), with Hajdu and Prudil as co-authors:
  - "only a single RRL is confirmed as belonging to a binary system"
  - "modern dedicated [RV] searches have come back empty handed (Barnes et al. 2021; Poretti et al. 2025)"
- **Prudil et al. 2025, Bulge Exploration V** ([arXiv:2506.19074](https://arxiv.org/abs/2506.19074)) derived systemic velocities for 8456 bulge RRL. It cites Hajdu 2021 only as a cause of phase shifts and does not test binarity. Whether any candidate has multi-epoch velocities in BRAVA-RR or APOGEE is **UNVERIFIED** (I found no paper that checks this).

## (b) Longer-baseline O-C

- **Hajdu et al. 2026**, "RR Lyrae stars with variable mean magnitudes" ([arXiv:2512.15636](https://arxiv.org/abs/2512.15636); [A&A](https://doi.org/10.1051/0004-6361/202557724)):
  - Adds OGLE-BLG-RRLYR-32226 (Appendix E): P_orb = 1869±19 d, a₁sin i = 0.202 AU, K₁ = 1.21 km/s, f(m) = 3.1×10⁻⁴ M☉. It had been discarded in 2021 because its mean-magnitude variation was mistaken for the Blazhko effect.
  - Finds 5 of 86 candidates among 71 RRab with variable mean magnitudes. The quoted chance probability is about 0.02% (about 3.5σ), which they call suggestive, not proof.
  - Does not re-fit the 87 against the new data.
- Hajdu 2021 is itself a longer-baseline test of earlier candidates. Of 38 previous bulge candidates (Hajdu+2015 and +2018, Prudil+2019 [arXiv:1905.11878](https://arxiv.org/abs/1905.11878)), **12 were discarded** for being LTTE-incompatible or Blazhko, and 26 were kept.
- Secondary, not refereed: Nguyen 2026 on Zenodo ([10.5281/zenodo.23048122](https://doi.org/10.5281/zenodo.23048122)) rebuilt the O-C curves from public OGLE data.
  - On pure noise, a negative ΔBIC favours an orbit in 72% of cases.
  - 68 of the 87 curves are not described by the stochastic-noise model.
  - Its own conclusion: "No candidate is confirmed and none is refuted."

## (c) Gaia DR3

- **Nagarajan et al. 2026** ([arXiv:2602.21289](https://arxiv.org/abs/2602.21289)):
  - Models predict about 202 RRL with DR3 orbital solutions; none are observed.
  - SB1 and EclipsingBinary solutions in the NSS catalogue for RRL are spurious.
  - Conclusion: "no genuine RR Lyrae received astrometric binary solutions in DR3."
- **Iorio et al. 2026** ([arXiv:2603.20429](https://arxiv.org/abs/2603.20429)): "no genuine RRL binaries have been robustly identified, including in the Gaia DR3 astrometric binary catalogues."
- Caveat: most bulge candidates have P_orb ≳ 3000 d, longer than the DR3 baseline of about 1000 d. DR3 therefore has little power on them, so this null result does not refute them.

## (d) Confirmation rate

- Confirmed: TU UMa only, a field star (Salinas+2026; Hajdu+2021 §1).
- Bulge confirmations: **0 of 87, or 0 of 88** with 32226.
- Refuted field LTTE claims: Z CVn by RV (Skarka et al. 2018, cited in Hajdu+2021 and Barnes+2021) and KIC 2831097 (Poretti+2025).
- IY Lyr, a field RRc, has an LTTE signal plus RV residuals and Gaia proper motions ([arXiv:2605.05708](https://arxiv.org/abs/2605.05708)). The authors call the companion "possible", not confirmed.

## Hajdu et al. 2021 itself ([arXiv:2105.03750](https://arxiv.org/abs/2105.03750); [ApJ 915, 50](https://doi.org/10.3847/1538-4357/abff4b))

- **Data:** OGLE-III and OGLE-IV I-band photometry, RRab only. The search used data to 2017; the fits added 2018–2019. 27,480 O-C diagrams were built, with mean (median) baseline 10.5 (7.4) yr. One O-C point per season, or per segment of a season.
- **Selection:**
  1. Visual inspection kept about 400 "even vaguely" LTTE-like curves; Blazhko stars were not excluded.
  2. Iterative LTTE+parabola fits, with bootstrap (500 resamplings) O-C errors.
  3. Stars were removed for amplitude modulation at P_orb, for departing from the fit when 2018–19 data were added, or when the curve was indistinguishable from a nonlinear period change. Fewer than 100 survived.
  4. Final emcee MCMC fits left 87.
- **No numerical detection threshold** (no false-alarm probability or ΔBIC) is given, and **no typical O-C residual or precision** is stated. Per-point bootstrap errors are shown only in the Fig. 1 panels.
- **Quality classes:** Q1 = 25 (σ(P_orb) < 1%), Q2 = 32, Q3 = 30 (problematic fits or data coverage).
- **Fitted parameters (Table 1; 0 / 16 / 50 / 84 / 100th percentiles):**
  - P_orb: 1076 / 2744 / 3718 / 5228 / 8819 d
  - a₁sin i: 0.14 / 0.39 / 0.99 / 2.34 / 3.97 AU, which is 0.8 / 2.3 / 5.7 / 13.5 / 22.9 × 10⁻³ d of light time
  - e: 0.06 / 0.18 / 0.30 / 0.51 / 0.78
  - K₁: 0.74 / 1.1 / 3.6 / 7.4 / 20.8 km/s
  - f(m): 10⁻⁴ / 4×10⁻⁴ / 0.012 / 0.12 / 0.67 M☉
  - These are computed from the parsed table `data/external/hajdu2021_binprop.csv`; light time = a₁sin i × 499.0 s/AU.
- **Orbital cycles covered:** not stated by the authors; my estimate follows.
  - The longest possible baseline is 2001–2019, about 19 yr ≈ 6900 d.
  - That gives about 1.9 cycles at the median P_orb (6900/3718).
  - Stars with OGLE-IV data only (2010–2019, about 3500 d) cover less than one cycle at the median P_orb.
  - The 3 stars with P_orb > 7000 d have less than one cycle even on the full baseline.
- **Machine-readable data:**
  - No VizieR entry for J/ApJ/915/50, and my guessed IOP `_mrt` URL returned 404.
  - Table 1 is in the arXiv LaTeX source (https://arxiv.org/src/2105.03750); a local copy is at `data/external/hajdu2021_arxiv2105.03750.tex`.
  - Code is at https://github.com/gerhajdu/rrl_binaries_1.

## Caveats

- ADS full-text search was not possible. A conference abstract or a non-arXiv paper could have been missed; the gap is the AAS240 abstract [BAAS 54(6) 415.04](https://baas.aas.org/pub/2022n6i415p04), which was blocked to me.
