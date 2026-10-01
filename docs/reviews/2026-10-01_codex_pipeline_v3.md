**I would not start the full refit on the evidence in the packet.** I found a real train/test leak in the predictive statistic and several ways the calibration can reuse information from the stars it evaluates. I did **not** find a row-order bug in Level 2, and the prior pseudo-rows are internally consistent under the stated Gaussian-prior assumption.

### Findings

1. **Predictive scores use held-out data during unwrapping — high severity for the predictive test.** `oc_stats` unwraps the entire series before defining the training and test sets; `robust_unwrap` uses neighboring seasons on both sides of each point. Thus an early held-out MACHO delay can affect the integer-cycle assignments of training delays, even though the period search itself uses training rows only. A test-season slip could change the selected training period and its reported predictive score. ([oc.py:229–232](src/rrlbin/oc.py:229), [oc.py:250–265](src/rrlbin/oc.py:250), [oc.py:26–55](src/rrlbin/oc.py:26))

2. **Band-lag priors are calibrated on stars that can receive those same priors — medium/high severity for candidate scores and calibration.** The calibration script fits offsets for stars with MACHO–OGLE overlap and uses those measurements to construct the lag relation. Level 2 then applies that relation to any star with the relevant band. There is no exclusion or leave-one-star-out step, so calibration stars are partly evaluated using priors informed by their own offsets. This can make their offsets appear better constrained than an independent prior would. ([calibrate_band_lag.py:38–59](scripts/calibrate_band_lag.py:38), [level2.py:51–55](scripts/level2.py:51), [level2.py:99–102](scripts/level2.py:99))

3. **The “empirical” simulation class reuses real H0 noise fits that may contain binary signal — medium severity for false-positive/completeness calibration.** The simulation script draws white jitter and random-walk amplitudes from a real-data noise table; its own docstring says those fits include any real binaries’ contribution. If a candidate orbit raises a star’s fitted noise, it can make the simulated nuisance population broader, potentially changing both false-positive and recovery rates. ([level1_sims.py:5–8](scripts/level1_sims.py:5), [level1_sims.py:50–57](scripts/level1_sims.py:50), [level1_sims.py:99–102](scripts/level1_sims.py:99))

4. **Common-mode estimation can absorb population-coherent candidate signals — medium severity, conditional on phase coherence.** The estimator uses selected real stars’ H0 residuals, and the shown invocation selects by season count rather than excluding candidate-like signals. Clipping limits outliers, but a signal shared by enough stars in a given band/year can enter the common-mode table and then be subtracted from those stars. That would attenuate real signals. ([level2.py:89–96](scripts/level2.py:89), [oc.py:281–316](src/rrlbin/oc.py:281), [oc.py:336–340](src/rrlbin/oc.py:336))

5. **Unwrapping can choose the wrong cycle for clusters of ambiguous seasons — medium severity for O−C fits.** The repair is an unweighted local linear fit to nearby delays, with no use of delay uncertainties or outlier-resistant fit. If multiple neighboring seasons have incorrect or noisy cycle assignments, they can support one another’s wrong branch; that can create or erase a smooth trend that Level 2 then fits. The visible unit test covers one isolated slip and a smooth quadratic drift, not clustered slips, large jumps, or low-S/N real-like series. ([oc.py:26–55](src/rrlbin/oc.py:26), [test_oc.py:17–29](tests/test_oc.py:17))

6. **Chunk resumption can silently reuse stale output — medium severity for a refit.** Existing `part_NNNN.parquet` files are skipped solely because they exist; the code does not check whether the input rows, options, or code version match the current run. Reusing an output directory after changing the sample, chunk size, or configuration can merge old and new parts without warning. ([chunked.py:21–35](scripts/chunked.py:21), [chunked.py:38–39](scripts/chunked.py:38))

7. **Fit failures are silently converted into missing bands/stars — medium severity for completeness accounting.** `fit_star` catches every exception and drops that band; `level1_real.one` can then mark a row `ok=True` even if the resulting series is empty or lacks OGLE I. Simulations explicitly reject a missing OGLE I fit instead. Level 2 later drops series with fewer than six points, but the real-data success count can therefore overstate usable output, and the real/sim failure paths differ. ([pipeline.py:38–48](src/rrlbin/pipeline.py:38), [pipeline.py:70–71](src/rrlbin/pipeline.py:70), [level1_real.py:38–47](scripts/level1_real.py:38), [level1_sims.py:42–46](scripts/level1_sims.py:42), [level2.py:83–88](scripts/level2.py:83))

### Areas I checked

- **Prior pseudo-rows:** The implementation adds the same offset-prior rows to H0 and H1, includes their Gaussian normalization in the likelihood, and gives them a separate covariance block in the red-noise fit. That is internally consistent if the priors are fixed, independent Gaussian information. The reuse of calibration stars above violates the “independent information” interpretation for those stars, but I did not find an H0/H1 mismatch in the code. ([oc.py:76–106](src/rrlbin/oc.py:76), [oc.py:160–176](src/rrlbin/oc.py:160))

- **Amplitude and two harmonics:** The reported two-harmonic amplitude is half the peak-to-peak range of the fitted two-harmonic curve; that is a coherent amplitude definition for that fitted curve. It is not automatically the fundamental’s amplitude or the exact Keplerian LTTE amplitude, so comparisons must use the same definition. ([oc.py:142–149](src/rrlbin/oc.py:142), [oc.py:233–242](src/rrlbin/oc.py:233))

- **Level 2 row ordering:** I found no row-misalignment bug in the shown path: records and metadata come from the same filtered dataframe, and `Pool.map` returns results in input order before concatenation. ([level2.py:83–88](scripts/level2.py:83), [level2.py:100–106](scripts/level2.py:100))

### Verified versus assumed

**Verified by reading the code:** the full-series unwrapping precedes the predictive split; lag calibration and application use the same population without a holdout; empirical noise is drawn from real H0 fits; common-mode estimation uses selected stars’ H0 residuals; chunk resumption skips existing parts without validating them; Level 2 preserves `Pool.map` order.

**Assumptions:** the refit may reuse an output directory; some candidate or binary signals may contribute to the empirical noise table; and enough real stars could share a timing pattern for common-mode estimation to remove part of it. The common-mode concern depends on that last condition; it is not proof that the current candidates are being attenuated.

### Strongest argument against this assessment

The packet reports shared Level 1 fitting for data and simulations and small-scale validation, including a null simulation check after adding band-offset draws. The Level 2 multiprocessing path also preserves order. Those are meaningful safeguards. They do not test predictive unwrapping leakage or calibration-star reuse, and the packet’s validation set is small relative to the proposed run. ([pipeline.py:38–48](src/rrlbin/pipeline.py:38), [level1_sims.py:68–80](scripts/level1_sims.py:68), [level2.py:100–106](scripts/level2.py:100))

### What I would check before the run

- Recompute predictive scores with cycle unwrapping performed using training seasons only; compare selected periods and scores.
- Calibrate lag priors on a disjoint star set, or use leave-one-star-out priors for calibration stars.
- Refit empirical noise distributions after excluding candidate-like stars, and compare resulting class pass rates.
- Estimate the common mode with candidate-like stars excluded and compare per-star corrections and candidate statistics.
- Exercise unwrapping on clustered slips, abrupt period changes, and low-S/N series drawn from the real cadence and errors.
- Ensure the refit uses fresh output directories, or add a manifest that rejects stale parts.
- Report usable Level 1 rows by required band and fit-failure reason, rather than relying on `ok` alone.