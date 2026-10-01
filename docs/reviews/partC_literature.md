# Part C: Literature check for the LMC RR Lyrae LTTE binary project

Date: 2026-10-01. How I searched: WebSearch, the arXiv API (abstract-field Boolean queries), arXiv PDFs and Crossref for the bibliographic data. ADS full-text search was **not accessible**: the ADS UI returned HTTP 405 and no API token was available. One residual check is still open: run the ADS full-text query `full:"RR Lyrae" full:("light travel" OR "light-time") full:(LMC OR Magellanic) year:2000-2026`.

## Executive summary
- **Q1:** No published systematic LTTE/O-C binary search among LMC or SMC RR Lyrae was found, up to Oct 2026. Only single-star LMC studies exist (BE Dor / MACHO J050918.7-695015). The closest analogues are Prudil+2019 and Hajdu+2021 for the Galactic bulge (OGLE), and Rathour+2024 for **Magellanic Cepheids** (197 LTTE candidates).
- **Q2:** Irregular, non-evolutionary period changes are common. They affect about 10% of LMC RRc stars in 6.5 yr of MACHO data, more than one third of RRL in M5, 57/134 (+23 slightly irregular) in M3, and most RRc in omega Cen. LMC Cepheids show the same: 16.5% (F) and 68% (1O) are irregular.
- **Q3:** No dedicated MACHO-timing systematics paper was found. Two practical notes are documented: MACHO epochs mark the **start of exposure**, so add 150 s and apply the heliocentric correction; and the data contain a spurious 1-yr instrumental term. Poleski 2008 found that MACHO- and OGLE-derived Pdot for LMC Cepheids are uncorrelated, with opposite signs more common than matching ones.
- **Q4:** Published LMC RRL RVs come from low or intermediate resolution, with per-epoch errors of **about 10-35 km/s** (mean 25 km/s). The RRab pulsation RV amplitude is about 60-70 km/s, and about 30-40 km/s for RRc. The expected LTTE orbital K is a few km/s, so existing LMC RRL RVs cannot test binarity.
- **Q5:** No Gaia astrometric or PMa companion work on LMC RRL was found. All PMa, speckle and Gaia-detectability work is Galactic.

## Q1. Binary searches among Magellanic RRL
- **Derekas+2004:** MACHO J050918.7-695015 (BE Dor), an RRc star in front of the LMC. MACHO+OGLE-III data give an O-C cycle of about 8 yr.
- **Li, Qian & Zhu 2022:** used DASCH+MACHO+OGLE+ASAS-SN+TESS. Concluded the changes are "quasi-periodic and abrupt" and rejected LTTE in favour of convection-magnetic interaction.
- **Derekas+2021 (TESS Sci. Conf. 2 abstract, Zenodo):** BE Dor follow-up; adds about 10 yr of data plus spectroscopy.
- **Soszynski+2016 (OCVS, more than 45,000 Magellanic RRL):** searched for eclipsing RRL. The candidates (for example OGLE-LMC-RRLYR-03541, Porb = 16.23 d) have Porb of 1.5-16 d. That is too short for a classical RRL, so they are probably binary evolution pulsators or blends. These details come from a search snippet of the full paper; I opened only the abstract.
- **No match on Magellanic RRL:** arXiv abstract queries combining "RR Lyrae" × (LMC|SMC|Magellanic) × (O-C | light-travel | companion | binarity | period change) returned nothing beyond the above. The only other hit was Alcock+2000, below.
- **Analogues:**
  - Prudil+2019: about 9000 bulge O-C diagrams (OGLE+KMTNet) gave 20 candidates with Porb 3-15 yr.
  - Hajdu+2021: 87 bulge candidates (61 new), with a trimodal mass function at about 0.6, 0.2 and 0.067 Msun.
  - Rathour+2024: more than 7200 OGLE Cepheids; 52 LMC and 145 SMC LTTE candidates; Porb mostly 2000-4000 d, e of 0.2-0.5.
- **Related 2026 OGLE-bulge work:** Hajdu+2026 report mean-magnitude variations in 72 RRL (about 0.9% of RRL). Abdollahi+2025 found no LTTE in TESS data for metal-rich RRL.

## Q2. Period-change statistics
- **LMC, Alcock+2000 (MACHO, more than 1300 RRc):** "10% of the stars show strong period changes" over 6.5 yr.
- **LMC, Alcock+2003 (6391 RRab):** 11.9% are Blazhko stars. Pdot statistics are not in the abstract.
- **LMC, Nagy & Kovacs 2006:** 7.5% (99/1332) of RRc are Blazhko stars. The paper separates period-change stars from Blazhko stars.
- **M3, Jurcsik+2012 (134 RRL, about 120 yr):** 54 regular, 23 slightly irregular, 57 irregular. About 50% of RRab are Blazhko stars.
- **M5, Szeidl+2011 (86 RRL):** more than one third show irregular Pdot (RRc 50%, RRab 34%). Irregularity correlates with the Blazhko effect. The mean Pdot is about -0.006 d/Myr.
- **omega Cen, Jurcsik+2001 (126 stars):** 44 RRab have monotonic Pdot with a mean of +0.15 d/Myr. The RRc stars are mostly irregular.
- **GEOS field RRab, Le Borgne+2007 (123 stars):** median beta is +0.14 d/Myr for the 27 increasing stars and -0.20 d/Myr for the 21 decreasing stars. Erratic O-C patterns are included.
- **Cautionary LTTE work:**
  - Liska+2016: 11 GEOS RRab stars modelled with LiTE.
  - Skarka+2018 (arXiv:1710.06709): Z CVn period instability mimics a companion. UNVERIFIED journal details: I saw only the arXiv listing and DOI 10.1093/mnras/stx2737.
- **LMC Cepheids:**
  - Poleski 2008: Pdot detected in 18% of F and 41% of 1O stars; random O-C fluctuations dominate on timescales of a few thousand days.
  - Rathour+2025: 33.5% irregular overall (F 16.5%, 1O 68.1%).

## Q3. MACHO timing
- **Michalska & Pigulski 2005, Sect. 2.2:** MACHO epochs "corresponded to the beginning of exposures, half the exposure time (150 s) was added". The heliocentric correction was applied separately, which implies the archive epochs are not heliocentric; check this in your files. They also removed a spurious 1-yr instrumental variation and many outliers.
- **Allsman & Axelrod 2001:** describes the archive but not the time standard (no HJD/UTC statement found).
- **Poleski 2008:** compared MACHO RM/VM with OGLE-II/III for 1515 Cepheids. Same versus opposite Pdot sign: OGLE vs RM 23/36 (F) and 53/87 (1O); OGLE vs VM 24/48 and 65/100. The paper gives no timing-offset diagnosis.
- **Gap:** no published MACHO-vs-OGLE zero-point timing test was found. Recommendation: test it yourself with stable RRab or EBs.

## Q4. RV precision for LMC RRL (V about 19-19.5)
- **Minniti+2003 / Borissova+2004 (VLT/FORS1, R about 1000, S/N about 15):** individual RVs good to 10-30 km/s. The internal line-to-line spread is 1-33 km/s. They quote RRab RV amplitudes of 60-70 km/s and RRc of 30-40 km/s.
- **Borissova+2006 (FORS1/2 with R = 780-1688, and GMOS):** mean per-star error 25±1 km/s. The repeat-epoch dispersion is 18±3 km/s (FORS) and 25±4 km/s (GMOS). The phase term is about 9 km/s.
- **Haschke+2012 (Magellan/MagE, R about 2000-6000, 6 LMC + 3 SMC):** RRab RV varies by up to 60-70 km/s over a cycle.
- **Galactic benchmark, Barnes+2021 (bright field RRL):** the gamma-velocity precision is 0.16-2.5 km/s (mean 0.92). This is not achievable at V about 19.5 without 8-m high-resolution spectroscopy and many epochs.
- **Requirement estimate (my calculation, not from the literature):**
  - K = 2π c A / (Porb √(1−e²)).
  - For a light-time semi-amplitude A = 0.01 d and Porb = 10 yr: a₁ sin i = 864 s × c ≈ 1.7 AU, giving K ≈ 5 km/s.
  - Confirming such an orbit needs ≲1-2 km/s per gamma epoch after removing the pulsation curve.
  - No published LMC RRL RVs at that level were found. A VLT/FLAMES-GIRAFFE or MUSE RRL RV paper on the LMC was not found: UNVERIFIED that none exists.

## Q5. Gaia/PMa
- **Kervella+2019 (Hipparcos-Gaia DR2 PMa):** 13/198 Galactic RRL significant plus 61 candidates; binary fraction at least 7%.
- **Salinas+2026 (speckle):** 10 companions in 81 RRL; binary fraction >12%.
- **Iorio+2026:** no robust Gaia RRL binaries among 100 candidates.
- **LMC RRL:** none found. At 50 kpc the expected photocentre wobble is ≲0.03 mas, so this is unsurprising.

## References (all URLs opened unless noted)
1. Derekas A. et al. 2004, MNRAS 354, 821. https://arxiv.org/abs/astro-ph/0407420 ; doi:10.1111/j.1365-2966.2004.08242.x
2. Li L.-J., Qian S.-B., Zhu L.-Y. 2022, MNRAS 510, 6050. https://arxiv.org/abs/2110.11615 ; doi:10.1093/mnras/stab3808
3. Derekas A. et al. 2021, TESS Sci. Conf. 2 abstract. https://zenodo.org/records/5123547
4. Soszynski I. et al. 2016, AcA 66, 131. https://arxiv.org/abs/1606.02727
5. Prudil Z. et al. 2019, MNRAS 487, L1. https://arxiv.org/abs/1905.11878 ; doi:10.1093/mnrasl/slz069
6. Hajdu G. et al. 2021, ApJ 915, 50. https://arxiv.org/abs/2105.03750 ; doi:10.3847/1538-4357/abff4b
7. Rathour R. S. et al. 2024, A&A 686, A268. https://arxiv.org/abs/2403.14039 ; doi:10.1051/0004-6361/202349117
8. Rathour R. S. et al. 2025, A&A 695, A114. https://arxiv.org/abs/2503.00661 ; doi:10.1051/0004-6361/202453392
9. Hajdu G. et al. 2026, A&A 707, A58. https://arxiv.org/abs/2512.15636 ; doi:10.1051/0004-6361/202557724
10. Abdollahi H., Molnar L., Varga V. 2025, A&A (accepted). https://arxiv.org/abs/2503.01018
11. Alcock C. et al. 2000, ApJ 542, 257. https://arxiv.org/abs/astro-ph/0005361
12. Alcock C. et al. 2003, ApJ 598, 597. https://arxiv.org/abs/astro-ph/0308019 ; doi:10.1086/378689
13. Nagy A., Kovacs G. 2006, A&A 454, 257. https://arxiv.org/abs/astro-ph/0602485 ; doi:10.1051/0004-6361:20054538
14. Jurcsik J. et al. 2012, MNRAS 419, 2173. https://arxiv.org/abs/1109.4525 ; doi:10.1111/j.1365-2966.2011.19868.x
15. Szeidl B. et al. 2011, MNRAS 411, 1744. https://arxiv.org/abs/1010.1115 ; doi:10.1111/j.1365-2966.2010.17815.x
16. Jurcsik J. et al. 2001, AJ 121, 951. doi:10.1086/318746 (Crossref-verified; I did not open the abstract page)
17. Le Borgne J. F. et al. 2007, A&A 476, 307. https://arxiv.org/abs/0710.1846 ; doi:10.1051/0004-6361:20077957
18. Liska J. et al. 2016, MNRAS 459, 4360. https://arxiv.org/abs/1504.05246 ; doi:10.1093/mnras/stw851
19. Skarka M. et al. 2018, MNRAS (arXiv:1710.06709). doi:10.1093/mnras/stx2737. UNVERIFIED volume and pages.
20. Poleski R. 2008, AcA 58, 313. https://arxiv.org/abs/0901.0884
21. Michalska G., Pigulski A. 2005, A&A 434, 89. https://arxiv.org/abs/astro-ph/0501380 ; doi:10.1051/0004-6361:20042343
22. Allsman R. A., Axelrod T. S. 2001, arXiv:astro-ph/0108444. https://arxiv.org/abs/astro-ph/0108444
23. Minniti D. et al. 2003, Science 301, 1508. https://arxiv.org/abs/astro-ph/0309351 ; doi:10.1126/science.1088529
24. Borissova J. et al. 2004, A&A 423, 97. https://arxiv.org/abs/astro-ph/0405377 ; doi:10.1051/0004-6361:20034494
25. Borissova J. et al. 2006, A&A 460, 459. https://arxiv.org/abs/astro-ph/0609209 ; doi:10.1051/0004-6361:20054132
26. Gratton R. G. et al. 2004, A&A 421, 937. https://arxiv.org/abs/astro-ph/0405412 ; doi:10.1051/0004-6361:20035840 (metallicities; RV errors not extracted)
27. Haschke R. et al. 2012, AJ 144, 88. https://arxiv.org/abs/1206.4999 ; doi:10.1088/0004-6256/144/3/88
28. Barnes T. G. III et al. 2021, AJ 162, 117. https://arxiv.org/abs/2106.05208 ; doi:10.3847/1538-3881/ac09f2
29. Kervella P. et al. 2019, A&A 623, A116. https://arxiv.org/abs/1903.03632 ; doi:10.1051/0004-6361/201834210
30. Salinas R. et al. 2026, A&A 711, A102. https://arxiv.org/abs/2603.28684 ; doi:10.1051/0004-6361/202659481
31. Iorio G. et al. 2026, A&A 712, A223. https://arxiv.org/abs/2603.20429 ; doi:10.1051/0004-6361/202659978
