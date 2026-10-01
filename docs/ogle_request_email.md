**To:** Igor Soszyński <soszynsk@astrouw.edu.pl>
**Cc:** Andrzej Udalski (OGLE PI) — *address to be confirmed*; optionally Przemek Mróz (lead of the 2022–24 high-cadence Magellanic analysis)
**Subject:** Request: post-2016 OGLE-IV photometry of LMC RR Lyrae for a light-travel-time binary search

Dear Igor,

I am writing to ask whether the OGLE team would be willing to share OGLE-IV photometry of the LMC RR Lyrae taken after
the end of the public OCVS light curves (April 2016).

**The project.** We are searching for binary companions to LMC RR Lyrae through the light-travel-time effect. We combine the
public OCVS light curves (OGLE-II/III/IV, 1997–2016; Soszyński et al. 2009, 2016) with MACHO B and R photometry from 1992–1999.
For 17,490 RRab we measure one pulsation delay per observing season. The pipeline uses a template per band, a common-mode
timing correction, band-lag priors to tie MACHO to OGLE, and a red-noise model of the intrinsic period wander. Every step is
calibrated with simulations on the real cadences.

**Results so far.** Over the 24-yr baseline we find 75 orbit-like O−C signals. Eleven are particularly clean: P_orb ≈ 5–15 yr,
coherent in MACHO and OGLE, about three cycles covered, consistent amplitudes and harmonic phases, and an out-of-sample
prediction of the MACHO seasons passed. Our noise simulations show, however, that the *number* of candidates is consistent with
the tail of the RR Lyrae's own irregular period changes. With the public data alone we can therefore only set an upper limit:
fewer than ~2.7% of LMC RRab have 0.4–1.5 Msun companions on 1–10 kd orbits. Telling real orbits from period wander needs more
cycles.

**Why your data would be decisive.** Before seeing any new data, we have frozen and dated the predicted 2016–2026 O−C curves
of the 28 strongest candidates (orbit vs. no-orbit, with red-noise uncertainties). The predictions are archived in our public
repository with a checksum and a git tag. With the 2017–2020 and 2022–2026 seasons, every one of these candidates would be
either confirmed or rejected at high significance. Seven of the 11 best lie in fields LMC502–516, so the 2022–24 high-cadence
data would give especially precise timings (~20 s per season, against ~130 s now).

**The request.** If possible, OGLE-IV I-band (and V, if convenient) light curves after 2016 for all RR Lyrae in the LMC OCVS
collection, in the OCVS format, with the HJD convention stated (and BJD-TDB if that is now standard). Having the full sample,
rather than just the candidates, lets us re-run the whole calibrated analysis: the noise model, the false-alarm rates and the
population limits. It also avoids biasing the test towards stars we already selected. If the full set is too much, the 75
candidates plus a random control sample of ~2000 RRab from the same fields would already go a long way.

We would of course treat the data as you prefer: kept within the project, with no redistribution and with any embargo you set.
We would also acknowledge OGLE appropriately, and we would be very happy to involve you and any interested members of the team
as co-authors. The code is open (https://github.com/vasilybelokurov/rrl_lmc_binary), and we can share the candidate list, the
frozen predictions and a short write-up of the method whenever useful.

Many thanks for considering this, and for the wonderful OCVS.

Best wishes,
Vasily Belokurov
Institute of Astronomy, University of Cambridge

---
*Notes for the sender (delete before sending):*
- The address soszynsk@astrouw.edu.pl is from the OCVS LMC RR Lyrae README. I have not verified addresses for Udalski or Mróz;
  please check them on the OGLE web pages.
- Facts used, with sources: the public RRL light curves end 2016-04-17; OGLE-IV observed the LMC until 2020-03 and again from
  2022-08; there was a high-cadence run on the central LMC in 2022-10 to 2024-05 (Mróz et al. 2024, arXiv:2410.06251). The list of
  high-cadence fields comes from a subagent's reading of that paper, and I have not checked it myself.
- Numbers: 75 candidates, 11 Tier 1, the upper limit and the decisiveness forecast are all from this session's v3 analysis and
  are provisional (JOURNAL.md, results/partB/, results/predictions/).
