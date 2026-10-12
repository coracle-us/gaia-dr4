# Euclid DR1 brown-dwarf companions: literature and overlap check

Written 2026-10-11 for Phase 2 scoping (see CLAUDE.md §2b).

**How this was checked.** A research agent opened every arXiv abstract marked [V]; [S] means seen only in search results or on a non-arXiv page. I checked the closest precedent myself (Mohandasan et al. 2025, §4.4, full text). Reference numbers are at the end.

## Bottom line

- **No published systematic search exists** for brown-dwarf companions to Gaia stars in Euclid data.
- **Two caveats:**
  - **The Euclid Consortium (EC) has a head start.** It has held DR1 internally since October 2025, with one year of proprietary access (Lorentz Center workshop report of the Milky Way & Resolved Stellar Populations working group). Its ultracool-dwarf (UCD) team has already run a small Gaia companion search on Q1 [2].
  - **Backyard Worlds is independently mining public Euclid data** [4], and turned Q1 into a paper within about 9 months.
- **The gap:** a specific, rigorous measurement that only DR1 + Gaia DR4 allow. That means faint late-T companions to Gaia DR4 stars at 50–200 pc, at small separations (1–10″), with calibrated chance-alignment statistics and companion-fraction limits. A generic catalogue of "new T dwarfs" or "new benchmarks" is not the gap.

## 1. Euclid brown-dwarf work so far

**Žerjal et al. 2026** [1], Q1 Deep Fields, 63 deg²
- Selection: a single colour cut, I_E − Y_E > 2.5, plus S/N.
- 5,306 UCD candidates (~1,200 L/T); 546 spectroscopically confirmed (329 L, 26 T).
- Only 13 photometric T candidates, so the selection is biased against T dwarfs.
- No proper motions ("repeated EDF visits" needed). No companions.
- Forecast for the full Wide Survey: ≥ 1.4 M UCDs, ≥ 300k L, > 2,600 T (lower limits).

**Mohandasan, Smart, Reylé et al. 2025** [2], EDF-N, 22.9 deg², J_E = 24.5
- 33 new UCDs (M7–T1) typed from NISP spectra.
- Proper motions from Euclid against AllWISE.
- **§4.4 "Possible companions"** (verified in the full text): Gaia DR3 search.
  - Criteria: spectroscopic parallax > 8 mas, projected separation < 200,000 au, parallaxes agreeing within 3σ, proper motions agreeing within 20%.
  - Candidate companions for 2 UCDs; one is called a probable false positive.
  - **This is the direct precedent, by the EC.**

**Dominguez-Tagle et al. 2025** [3]: NISP spectral classification of hundreds of UCDs, including ≥ 11 new late-L/T. No companions.

**Kiwy, Kirkpatrick, Schneider, Meisner, Faherty, Kuchner et al. 2025** [4] (Backyard Worlds, not EC)
- Y_E − J_E vs J_E − H_E colour selection plus NISP spectra.
- 15 validated mid/late-T dwarfs, 8 of them new.

**Early Release Observations, young regions**
- σ Ori, Martín et al. 2025 [6]: no mass-function cutoff down to ~4 M_J; used Gaia.
- Taurus LDN 1495, Bouy et al. 2025 [7]: proper motions from Euclid against more than 20 years of archival ground-based imaging. That is our technique.
- Barnard 30, Barrado et al. 2026 [8].

**Other items**
- **First field brown-dwarf binary resolved by Euclid:** L/T transition, ~129 pc, ~66 au (talk at BD30, September 2025). UNVERIFIED; no paper found.
- **Forecasts:**
  - Solano et al. 2021 [9]: ~1–2 M L, 0.5–1 M T and a handful of Y dwarfs over 15,000 deg².
  - Marocco et al. 2017 [10]: ~24,000 resolved UCD companions to Gaia primaries at |b| > 10°.
  - Martín (2020 talk): "tens of benchmark young wide substellar companions" [S].
- **Bedin 2025** [11] proposes repeating the Wide Survey for proper motions in an extended mission (not approved). **DR1 is effectively single-epoch.**

## 2. Consortium plans

- **ESA DR1 timeline** (15 June 2026) [29]:
  - DR1-Foundation, November 2026: ~1,900 deg² of raw data, calibrated images, catalogues and spectra.
  - Full DR1, mid-2027: higher-level clustering and weak-lensing products.
  - Details are in `notes/euclid_scoping.md`.
- **EC DR1 Key Project list** (73 projects, euclid-ec.org, updated 16 September 2025) [29]:
  - The words "dwarf", "ultracool", "brown" and "substellar" do not appear.
  - The Milky Way working group has 3 Key Projects: globular clusters, resolved stellar peripheries, and the Magellanic Clouds.
  - A UCD Key Project is absent, but EC UCD papers are likely anyway. The same group (Martín, Žerjal, Zapatero Osorio, Smart, Reylé, Dominguez-Tagle, Muñoz Torres, Zhang J.-Y.) wrote the Q1 papers.
  - No public list of their DR1 papers exists.

## 3. Existing wide brown-dwarf-companion searches

| Work | Data / depth | Area | Result |
|---|---|---|---|
| Faherty+2010 [13] | SDSS / 2MASS | — | 9 systems (common proper motion) |
| Burningham+2013 [14] | UKIDSS LAS (K ≈ 18) | ~7,000 deg² | 2 benchmark T companions |
| Deacon+2014 [15] | Pan-STARRS1 | 3π sky | 61 new companions (K7–L9) |
| Marocco+2017 [10] | UKIDSS + SDSS × TGAS | — | 13 benchmarks (false-alarm probability < 1e-4) |
| dal Ponte+2020 [17] | DES DR1 + Gaia DR2 | ~5,000 deg² | 255 candidates; wide fraction 2–4% |
| Marocco+2021, CatWISE2020 [18] | W1 17.7 / W2 17.5 (90% complete) | all sky | motion precision ~20 mas/yr at W1 = 15, ~100 at W1 = 17 |
| Kirkpatrick+2021 [19] | 20 pc census | — | 525 L/T/Y; complete except ≥ T8 |
| El-Badry, Rix & Heintz 2021 [20] | Gaia eDR3 | ~1 kpc | ~1.3 M pairs (both components seen by Gaia) |
| Rothermich+2024 [22] | Backyard Worlds; Gaia × CatWISE | all sky | 89 new companions (M7–T9) |
| Baig+2024 [23] | Gaia UCD sample, 100 pc | — | 278 systems |
| Zhang Z.H.+2025 [24] | CatWISE + Gaia | — | benchmark brown dwarfs (series) |
| Schneider+2026 [25] | WISE / Backyard Worlds | all sky | 3,006 L/T; 28 new comoving companions |
| Cheng+2026 [28] [S] | VHS + DES | — | 10 T dwarfs with Gaia hosts (5 new) |

Notes on the table:
- **dal Ponte+2020** matched photometric distance to the host and estimated chance-alignment probabilities; 73% also pass common proper motion. This is the closest methodological template for our faint regime.
- **Cheng+2026** used 3–13 yr proper-motion baselines.
- **Rubin:** Honaker et al. 2026 [26] (89 UCDs in LSST DP1) includes Žerjal and Martín, so the EC UCD team is also positioned on LSST.

## 4. What DR1 newly enables

**Depth** (5σ point source):

| Survey | Depth |
|---|---|
| CatWISE | W2 ≈ 17.5 (Vega, 90% completeness) |
| VHS | J ≈ 20.2, Ks ≈ 18.1 (Vega, design) [S] |
| UKIDSS LAS | K ≈ 18 |
| Euclid Wide | Y_E / J_E / H_E ≈ 24.5 AB (J_Vega ≈ 23.6); I_E ≈ 26.2 AB [12] |

**Distance reach** (rough, using M_J ≈ 15 for T5, 16.3 for T8, ~19.5 for Y0; detailed table in `notes/euclid_scoping.md`):

| Type | VHS | CatWISE | Euclid |
|---|---|---|---|
| T5 | ~110 pc | ~90 pc | ~500 pc |
| T8 | ~60 pc | ~50 pc | ~280 pc |
| Y0 | ~15 pc | — | ~60–80 pc |

**Resolution:** VIS at 0.1″ per pixel opens 1–10″ separations (100–2,000 au at 100–200 pc), where WISE (~6″ PSF) and ground-based seeing fail.

**Catch:** no earlier survey sees T dwarfs at J > 21, so common proper motion is possible only for J ≲ 20–21, exactly where competitors are strongest. Fainter companions need photometric or spectroscopic distance plus statistics.

## 5. Where we can contribute (and what to avoid)

**Go for:**
1. **Faint, close companion statistics.** T5–T9 companions at 1–10″ around Gaia DR4 stars at 50–200 pc. Use Euclid photometric distances against DR4 parallaxes and an empirical chance-alignment calibration (El-Badry-style shifted-position control samples). Deliver a wide-companion fraction with upper limits. CatWISE and VHS cannot do this.
2. **Host-class niches enabled by Gaia DR4:**
   - white-dwarf hosts (ages);
   - metal-poor / thick-disc hosts;
   - **exoplanet hosts.** This ties directly to Phase 1: the same transiting-planet host list, giving brown-dwarf companions plus our DR4 astrometric limits as a combined architecture constraint.
3. **Proper-motion confirmation** of the brighter subset: Euclid against Legacy Surveys / DES / VHS / UKIDSS over 7–12 yr baselines (Bouy et al. 2025 technique).

**Avoid:**
- generic DR1 T-dwarf catalogues (EC UCD team, Backyard Worlds);
- broad "benchmark companion" discovery lists (Hertfordshire group);
- anything that needs the EC's proprietary head start.

**Expect fast competition.** Speed matters, and so does the niche.

## References

1. Žerjal M. et al. 2026, A&A 711, A37, arXiv:2503.22497
2. Mohandasan A., Smart R.L., Reylé C. et al. 2025, arXiv:2503.22559 (§4.4 verified)
3. Dominguez-Tagle C. et al. 2025, ApJ, arXiv:2503.22442
4. Kiwy F. et al. 2025, arXiv:2512.17494
5. Zhang J.-Y., Lodieu N., Martín E.L. 2024, A&A 686, A171, arXiv:2403.15288 [S]
6. Martín E.L. et al. 2025, A&A 697, A7, arXiv:2405.13497
7. Bouy H., Martín E.L. et al. 2025, A&A 696, A80, arXiv:2502.16349
8. Barrado D. et al. 2026, A&A (accepted), arXiv:2608.07048
9. Solano E. et al. 2021, MNRAS, arXiv:2010.16392
10. Marocco F. et al. 2017, MNRAS, arXiv:1706.06038
11. Bedin L.R. 2025, A&A 704, A193, arXiv:2510.23694
12. Scaramella R. et al. 2022, A&A 662, A112, arXiv:2108.01201
13. Faherty J.K. et al. 2010, AJ 139, 176, arXiv:0911.1363
14. Burningham B. et al. 2013, MNRAS, arXiv:1304.7246
15. Deacon N.R. et al. 2014, ApJ 792, 119, arXiv:1407.2938
17. dal Ponte M. et al. 2020, MNRAS 499, 5302, arXiv:2001.11015
18. Marocco F. et al. 2021, ApJS, arXiv:2012.13084
19. Kirkpatrick J.D. et al. 2021, ApJS, arXiv:2011.11616
20. El-Badry K., Rix H.-W., Heintz T.M. 2021, MNRAS, arXiv:2101.05282
22. Rothermich A. et al. 2024, AJ, arXiv:2403.04592
23. Baig S. et al. 2024, MNRAS 533, 3784, arXiv:2408.07024
24. Zhang Z.H. et al. 2025, MNRAS 542, 656, arXiv:2505.24560
25. Schneider A.C. et al. 2026, AJ (accepted), arXiv:2604.01323
26. Honaker E.J. et al. 2026, arXiv:2604.08387
28. Cheng G. et al. 2026, Cool Stars 23 poster, Zenodo 20181679 [S]
29. ESA Euclid DR1 timeline page (15 June 2026); EC DR1 Key Project list, euclid-ec.org/?p=4955 (16 September 2025); Lorentz Center workshop report (2026)
