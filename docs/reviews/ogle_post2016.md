# OGLE-IV LMC monitoring after 2016 (checked 2026-10-01)

## Summary
- OGLE-IV observed the LMC continuously from March 2010 to 15 March 2020. It stopped for COVID-19 and **restarted on 12 Aug 2022**. It is **still running in 2026**: XROM LMC light curves run to 2026-05-26, and EWS bulge alerts reach 2026-BLG-0914.
- Regular LMC fields get about **30–50 I-band epochs per season** from 2014/15 onward, including after 2022. In 2010–2013 they got about 110–270 per season.
- From Oct 2022 a **high-cadence campaign** observed selected central fields at about 20-min cadence, roughly 1500–2400 epochs per season.
- **Public post-2016 LMC photometry is sparse.** Most OCVS RR Lyrae and Cepheid light curves end in **April 2016**. The OCVS eclipsing-binary release (2025) runs to **May 2024**. Data after 2022 for other stars has been "non-public", provided on request (Rathour+2024).
- I found no evidence of an "OGLE-V" or a camera upgrade. This is UNVERIFIED either way.

## Timeline
| Date | Event | Source |
|---|---|---|
| 2010 Mar 4–5 | First OGLE-IV LMC epochs (HJD 2455260.6) | XROM/OCVS files, tested below; Glowacki+2025 says "between March 2010 and May 2024" |
| 2016 Apr 17 | Last epoch in OCVS LMC RR Lyrae light curves (Soszyński+2016) | tested below |
| 2020 Mar 15 | Last pre-COVID LMC epoch (JD 2458924) | Mróz+2024 ApJS |
| 2022 Aug 12 | "OGLE-IV project resumed regular observations after 2.4 year long break" | OGLE homepage HTML notice |
| 2022 Oct 4 – 2024 Apr 27 | High-cadence run on 5 LMC fields (LMC502/503/509/510/516); 4868–4880 epochs each; median 20 min, ≤23 visits/night | Mróz+2024 ApJL |
| 2024/25–2025/26 | Continued high cadence on other fields (LMC531, LMC552: ~2060–2200 epochs/season) | XROM, tested |
| 2026 | Homepage shows "OGLE-IV IN OPERATION"; LMC data to 2026-05-26 | OGLE homepage, XROM |

## Data sources
Column definitions: Span = epochs covered; Cadence = LMC I-band epochs; Public = whether light curves can be downloaded.

| Source | Span | Cadence (LMC, I) | Public? |
|---|---|---|---|
| OCVS LMC RR Lyrae (Soszyński+2016, AcA 66,131) | 2010-03 → 2016-04-17 | ~400 epochs per star (3 sampled) | Yes (ftp) |
| OCVS LMC classical Cepheids (Soszyński+2015) | 2010-03 → 2016-04 for 4645/4656 stars; 11 stars run to 2020-03 | — | Yes |
| OCVS LMC T2Cep / δ Sct (2023) | to 2017-12 / 2020-03-11 | δ Sct: ~130 epochs after 2017 | Yes |
| OCVS MC eclipsing binaries (Glowacki+2025) | 2010-03 → 2024-05 | LMC median 879 epochs, maximum 5808 | Yes; ~75–220 epochs after Aug 2022 in samples |
| Mróz+2024 Nature/ApJS (20-yr LMC microlensing) | OGLE-III 2001–2009 plus OGLE-IV 2010-06-29 → 2020-03-15 | 121–911 epochs per field (median 368); 3–10 d cadence | Only photometry of the event candidates |
| Mróz+2024 ApJL (high-cadence PBH search) | 2022-08-09 → 2024-05 | see timeline | Only the event candidates (ftp LMC_FFP_PBH) |
| XROM (19 LMC X-ray counterparts) | 2010-03 → 2026-05-26 | 30–50 per season; ~1500–2400 in high-cadence fields | Downloadable, but real-time with zero-point ±0.15 mag; "contact us" before publishing |
| Rathour+2024 A&A 686, A268 (MC Cepheid binaries) | 1997–2020 public data plus ~1.3 yr non-public data after 2022-08-12 | — | Data after 2022 came from the OGLE team |
| OGLE-IV Transient Detection System page | last updated 2020-03-16 | — | Not resumed publicly |

**Tests run.** I downloaded OCVS light curves (3 per class, plus the full LMC Cepheid tarball) and all 19 XROM LMC `phot.dat` files. Then I computed min/max HJD and epochs per season. Example from XROM, field LMC530.18: 2016/17: 40, 2018/19: 58, 2022/23: 50, 2024/25: 48, 2025/26: 46. Two fields, LMC510.17 and LMC511.17, have no XROM data after 2019-10. That may just be how XROM is maintained, so it is not evidence those fields stopped being observed (UNVERIFIED).

## Public availability
- **OCVS LMC RR Lyrae:** April 2016 is the end of the public time series. The `phot.tar.gz` was re-dated 2023-09-25, but the light curves I sampled still end in 2016.
- **Public data after 2016:** only in newer OCVS classes (δ Sct to 2020; ecl to 2024-05), the microlensing-event files and XROM.
- **Anything else after 2016:** request it from the OGLE team, as Rathour+2024 did.
- **VizieR mirrors:** not checked (UNVERIFIED).

## Not verified
- Soszyński+2019 "AcA 69, 87" as an extension of the RR Lyrae catalogue. arXiv:2001.00025 is AcA 69, 321 (bulge/disk), not the Magellanic Clouds.
- Exact OGLE-IV start date quoted from Udalski+2015. The abstract I opened gives no date.
- The OGLE Magellanic photometric maps.
- Any post-2016 RR Lyrae light-curve update.

## References (all opened)
- OGLE homepage HTML: https://ogle.astrouw.edu.pl/main/
- Udalski+2015, AcA 65, 1: https://arxiv.org/abs/1504.05966
- Mróz+2024 Nature 632, 749: https://arxiv.org/abs/2403.02386
- Mróz+2024 ApJS 273, 4: https://arxiv.org/abs/2403.02398 ; data https://www.astrouw.edu.pl/ogle/ogle4/LMC_OPTICAL_DEPTH/
- Mróz+2024 ApJL (PBH high-cadence survey): https://arxiv.org/abs/2410.06251 ; data https://ftp.astrouw.edu.pl/ogle/ogle4/LMC_FFP_PBH/
- Glowacki+2025: https://arxiv.org/abs/2503.15596
- Rathour+2024 A&A 686, A268: https://arxiv.org/abs/2403.14039v2
- Marković+2026 (WR stars; OGLE-IV "2010–2025", COVID gap HJD 2458900–2459800): https://arxiv.org/html/2602.05820v1
- OCVS LMC ftp: https://ftp.astrouw.edu.pl/ogle/ogle4/OCVS/lmc/
- XROM: https://ogle.astrouw.edu.pl/ogle4/xrom/xrom.html
- EWS 2026: https://ogle.astrouw.edu.pl/ogle4/ews/2026/ews.html
