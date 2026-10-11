# Literature review: outer companions to transiting giant-planet hosts (Phase 0, step 4)

Written 2026-10-11.

**How this was checked.** Research agents opened each paper's arXiv abstract, journal page or full text. Several key numbers were spot-checked here against the extracted PDF text:
- Lammers & Winn's thresholds
- Ranshaw+26's DR4 forecast
- Espinoza-Retamal+23's yields
- Greenfield & Owen's turnover

Items not verified are marked **UNVERIFIED**. Numbered references are in §8.

## 0. What this means for the project

1. **The rate we are testing.**
   - RV surveys give 50–75% of hot-Jupiter (HJ) hosts an outer companion of 1–20 M_J at 1–20 AU [1,2].
   - The field rate is about 15–20% [19,21].
   - DR4 can turn M sin i into true mass for the nearest and brightest hosts. Over our whole sample it can only constrain this occurrence through completeness-corrected upper limits.
2. **Mutual inclinations are now an active question.**
   - The pre-DR4 population result is moderately coplanar: Rayleigh σ ≈ 16° [30], and σ ≈ 12° for super-Earths [28].
   - The one HJ with a direct measurement is HD 118203, at about 10° [34].
   - Coplanar high-e migration predicts ΔI ≲ 20°; planet–planet secular migration predicts about 40°; Kozai by stellar companions needs ≳ 40° [14,37–39].
3. **Our sensitivity is revised down.**
   - Lammers & Winn's DR4 criterion [L1] is Δχ² > 50 with P ≲ 4 yr.
   - Companions aligned with a transiting planet are **edge-on**, the worst case: they give half the face-on Δχ².
   - On that basis, 73 hosts are sensitive to ≤ 5 M_J at P = 4 yr (a ≈ 2–2.6 AU), 370 to ≤ 13 M_J, and only 6 to ≤ 2 M_J.
   - These counts come from the `mmin_mj_edge_*` columns of `data/targets.csv`. This supersedes the face-on figure of ~125 quoted earlier.
4. **The known companions mostly sit at P > 4 yr.**
   - In `data/validation_known_companions.csv`, K2-232 c, WASP-132 d, TOI-2537 c, WASP-53 c, HAT-P-17 c, HD 118203 c, HAT-P-2 c and WASP-8 c all give large Δχ² (10²–10⁴).
   - Their periods are 5–23 yr, so DR4 alone sees curvature or a partial orbit.
   - Only WASP-81 c (a brown dwarf, P = 3.6 yr) is a clean astrometry-only orbit.
   - **Phase 1 design consequence:** the pipeline needs, besides the P ≤ 4 yr periodogram,
     - (a) acceleration (7- and 9-parameter) models;
     - (b) a joint astrometry + literature-RV mode for systems with a known RV companion or trend.
   - Ranshaw+26's ~75% forecast for K2-232 c [16] assumes exactly that joint fit.
5. **False positives are the main risk** (§5). Most DR3 astrometric "planets" around faint stars were near-twin binaries. Four DR3 orbits were retracted for a software bug. Long-period binaries alias into short-period orbits. Every claim needs independent refitting plus epoch-RV and colour–magnitude checks.
6. **Competition** (see also `dr4_paper_overlap.md`):
   - Xiao, Feng et al. [30] plan DR4 population ΔI.
   - Ranshaw, S. Wang et al. [16] have a published DR4 forecast for K2-232.
   - Espinoza-Retamal, Zhu & Petrovich [E1] already forecast our exact experiment for TESS hosts.
   - The Do Ó+26 candidate list [27] contains **no transiting giant-planet hosts**: its transiting hosts all have small planets.

## 1. Companion occurrence (RV, imaging)

**Hot-Jupiter hosts**
- **Knutson+14 (Friends of Hot Jupiters I)** [1]. Sample: 51 transiting HJs.
  - 51 ± 10% have companions of 1–13 M_J at 1–20 AU.
  - No difference between misaligned/eccentric and aligned systems.
- **Bryan+16** [2]. Rate 52 ± 5% for 1–20 M_J at 5–20 AU. Split by the inner planet's orbit:
  - hot (< 0.1 AU): 75 (+4/−6)%;
  - warm (0.1–1 AU): 49 ± 10% (2.4σ difference, tentative);
  - cold (1–5 AU): 54 (+7/−8)%.
- **Zink & Howard 2023** [7], California Legacy Survey.
  - Completeness-corrected outer giants to HJs are "ubiquitous" and at least 3× the HJ mass.
  - They are more eccentric: mean e = 0.34 ± 0.05, versus 0.19 ± 0.02 without an HJ.
  - They argue for coplanar high-e migration.

**Stellar companions to HJ hosts**
- **Ngo+15/16** [3,5]: 47 ± 7% have stellar companions at 50–2000 AU, 2.9× field. Only about 4% have one at 1–50 AU, versus 16% in the field. Fewer than 20% of HJs have a stellar companion capable of Kozai–Lidov migration.
- **Shariat+26** [8]: 147 HJs with obliquities; intrinsic stellar-companion fraction 62 ± 5% at 50–50,000 AU. Hosts with companions are misaligned 46% of the time, versus 24% without.

**Warm Jupiters (WJs)**
- About half have close small neighbours [11]. For comparison, 12 ± 6% of HJs and ≥ 70% of WJs have nearby companions from TTVs [9].
- Kozai migration of a WJ needs a Jupiter-mass perturber within about 3 AU [12].
- Planet-driven secular migration could make about 20% of WJs [14].
- Ranshaw+26 [16]: outer companions are more eccentric in HJ systems than in WJ systems, 0.37 ± 0.02 versus 0.22 ± 0.01. The headline significance is 5.75σ, but an Anderson–Darling test gives only p = 0.04.

**Field baseline**
- Fulton+21 [19]: giant occurrence peaks at 1–10 AU, with 14 ± 2 giants per 100 stars at 2–8 AU.
- Rosenthal+22 [21]: 17.6 (+2.4/−1.9)% of all stars have an outer giant.

## 2. Mutual inclinations and theory

**Population results**
- Masuda+20 [28]: super-Earths with cold Jupiters, Rayleigh σ ≈ 12°.
- Becker+17 [29]: outer companions to HJs around cool stars have σ ≲ 20°.
- Xiao+26 [30]: 19 transiting systems with RV plus Hipparcos–Gaia (HG) data. σ = 15.8 (+2.8/−2.6)°, with isotropic rejected; giant-inner-planet systems are possibly more coplanar.

**Individual systems**
- HD 118203 c: 11.8 M_J at 6.3 AU, ΔI = 10 (+16/−9)° [34].
- Kepler-419 c, a photodynamical measurement: 9 (+8/−6)° [35].
- π Men (inner super-Earth): misaligned, ΔI ≈ 35–140° [31–33].
- HAT-P-11 (inner Neptune): misaligned [31].

**What each migration channel predicts**
- Disk migration: ΔI ≈ 0.
- Coplanar high-e migration [39]: companion at ≳ 5 AU, 1–3× the HJ mass, ΔI ≲ 20°, e ≈ 0.2–0.5.
- Planet–planet secular migration [14]: ΔI clusters near 40° for HJs and ≳ 60° for WJs.
- Stellar Kozai [37,38]: needs ≳ 40°.
- Review: Dawson & Johnson 2018 [41].

**Measurement note.** For transiting systems, astrometry gives the companion's i_c but not its node, so the line-of-sight value |i_c − 90°| is a *minimum* ΔI. A population analysis (as in [30]) is needed to infer ΔI itself.

## 3. Pre-DR4 astrometry of planet hosts

**Catalogues**
- **HGCA EDR3** [A1]: VizieR `J/ApJS/254/42/catalog`, 115,346 stars. The column `chi2` is for a constant proper motion with 2 d.o.f., so 11.8 means 3σ and 28.7 means 5σ. Gaia errors are inflated by 1.37.
- **Kervella+22 proper-motion anomaly (PMa)** [A3]: VizieR `J/A+A/657/A7/tablea1`.
  - `snrPMaH2EG3b` is the PMa signal-to-noise; `BinH2EG3b = 1` when it exceeds 3.
  - `M23au`, `M25au` and `M210au` give the companion mass needed at 3, 5 and 10 AU.

**Gaia DR3 non-single-star (NSS) solutions**
- Pipeline papers: Halbwachs+23 [A5]; Holl+23 [A6], the exoplanet pipeline with 1,162 solutions, 198 "validated".
- Arenou+23 [A7]:
  - 72 candidates of ≤ 20 M_J, 9 known planets validated.
  - HD 114762 b and HD 162020 b turn out to be stars.
  - WASP-18 b was recovered as an SB1 from Gaia RVs.
- **ESA known issue, 2024-05-27** [A8]: 4 false-positive orbits from a software bug, found while validating the DR4 epoch astrometry. The affected systems are WD 0141-675, HIP 64690, 54 Cas and HIP 66074.

**Follow-up of DR3 candidates**
- Winn 2022 [A9]: of 7 RV planets, 3 consistent and 4 inconsistent.
- Marcussen & Albrecht 2023 [A10]: Gaia and RV disagree in 6 of 12 systems.
- Marcussen+26 [A11]: 13 of 20 are near-twin binaries.
- Stefansson+25 [A12]: 21 of 28 M/K candidates are stellar binaries. This paper confirmed Gaia-4b (11.8 M_J, P = 571 d), which is in the DR4 pre-release sample, and Gaia-5b (a brown dwarf).
- Sahlmann & Gómez 2025 [A16]: long-period binaries mimic short-period orbits.
- An+25 [A18]: DR3 accelerations are offset from predictions by a median 1.85σ (tail to 10σ); the DR3 long-period orbits are spurious.

**HG companion studies**
- Feng+22 [A25]: 167 cold giants.
- Thompson+26, the G23H catalogue [A26]: supports 94 of 120 Jovian systems.
- Do Ó+26 [27]: 29 outer-companion candidates. **None is a transiting giant-planet host.** The transiting hosts in the list (ν² Lup, LHS 1903, HD 63433, HD 96735) all have small planets.
- de Beurs+23: HAT-P-2 c is 10.7 M_J [A29].
- Yang+26: KELT-20 has a companion of about 34 M_J [32].

**Hosts that are binaries**
- Belokurov+20 [A33]: hints of elevated RUWE among massive HJ hosts.
- Ziegler+20/21/26 [A34–36]: close binaries (< 100 AU) are suppressed about 3–7× among TESS hosts.
- Michel & Mugrauer 2024 [A41]: 16–19% of hosts are multiples.

## 4. DR4 forecasts and how to calibrate our sensitivity

**Lammers & Winn 2025** [L1], arXiv:2511.04673 (AJ). Read in full; the key numbers were checked here.
- **Noise:** σ_fov = 54 µas · max(1, 10^{0.2(G−14)}), fit to Holl+23 [A6]. Their mocks use `gaiamock` (σ_ccd/√8, the same curve as ours), with 10% of epochs rejected and U(0, 40) µas extra noise for G < 13.
- **Detection:** Δχ² > 50 for the 12- versus 5-parameter model. Calibrated median relations: Δχ² ≈ 22·SNR1² (SNR1 = α/σ_fov, about N/3 for N ≈ 66 epochs, i.e. averaged over inclination); 20% period accuracy fails beyond P ≈ 4.0 yr in DR4. Pure-noise fits reach Δχ² > 30 in 5% of cases.
- **Yields:**
  - 7,500 ± 2,100 DR4 planets, 1,900 ± 540 of them well characterized.
  - Typical detection: ≳ 3 M_J at 1–3 AU, 50–250 pc.
  - Impostor binaries: about 13% of candidates overall; about 50% for M < 0.6 M☉ hosts at a < 1 AU; ≲ 10% contamination for G < 12.
- **Mock catalogues:** `CalebLammers/GaiaForecasts`, MIT licence. Files:
  - `DR4_mock_exoplanet_catalog.csv`: 7,545 rows, real DR3 source_ids;
  - `DR4_mock_planet_impostor_catalog.csv`: 1,151 rows;
  - `DR4_master_orbital_fits.csv`: 29,129 rows.

**Other forecasts**
- **Espinoza-Retamal, Zhu & Petrovich 2023** [E1], AJ 166, 231, arXiv:2309.08665. **This is the closest prior work: our experiment, forecast for TESS hosts.**
  - Assumed cold-Jupiter rates: 0.75 for HJ hosts, 0.49 for WJ hosts.
  - About 206 detections, with about 118 having cos i measured well enough (σ_cos i < 0.2); 46 of those are around large planets.
  - Optimistic: 34 µas noise, a Fisher-matrix method, and no period cap.
- **Greenfield & Owen 2026** [E2], arXiv:2610.07149, submitted. Detectability turns over at about 3 AU (about 5 yr) in DR4. Real BH3 epochs give 46 µas per transit, matching the Holl curve at G ≈ 11.
- **Perryman+14** [E3]: 34 µas per FoV for G ≤ 12, about 69 FoV transits per 5 yr, Δχ² > 30 "marginal" and > 50 "reliable". They noted that Gaia might clarify outer companions to HJs.
- **El-Badry+24 `gaiamock`** [A19]: forecasts about 1 million DR4 orbits, and that spurious orbits come mainly from marginally resolved pairs.

**Our model versus theirs** (`scripts/build_targets.py`)
- Same σ_ccd curve.
- Our N_DR4 = 1.75 × the DR3 matched-transit count, consistent with 1.96 × 0.9 rejection.
- Our old `snr_*` columns are the face-on best case (Δχ² = snr²). The new `mmin_mj_edge_*` columns use edge-on Δχ² = N α²/(4σ²) = 50 at P = 4 yr, matching their criterion.
- Our σ is about 20% smaller than their semi-analytic 54 µas at G ≈ 12, and the pre-release data suggest the real errors are smaller still (`dr4_data_model.md`). Injection–recovery in Phase 1 decides.

## 5. False-positive lessons, turned into vetting rules

1. **Near-twin binaries** cause most faint-star "planets" [A10–A12]. Check:
   - DR4 epoch RVs and SB2 flags;
   - colour–magnitude position (over-luminous for its colour);
   - `nss_masses.exoplanet_candidate_hrd_proba`.
2. **Pipeline bugs** [A8]: refit independently from the epochs; never adopt a catalogue orbit unchecked.
3. **Long-period aliasing** [A16,A18]: always compare against 7- and 9-parameter acceleration models, and use literature RV trends.
4. **Partially resolved pairs** [A19]: check neighbours, `blended`/`multipeak`, `nss_resolved_pair` and `optical_pair`.
5. **Biased e and i in real planets** (HIP 66074, LHS 1610) [A13,A17]: be cautious with ΔI from a single Gaia fit, and prefer joint RV fits.
6. **Underestimated errors** [A1,A18]: calibrate significance on null injections and quiet stars, not formal χ².
7. **Stellar activity and variability** inflate RUWE [A33]. HAT-P-11's RV trend tracks activity [1].

## 6. Our cross-match: pre-DR4 companion evidence for our hosts

`scripts/crossmatch_prior_astrometry.py` → `data/target_prior_astrometry.csv` (781 hosts with DR3 IDs).

**Summary counts**

| Evidence | Hosts |
|---|---|
| RUWE > 1.4 | 52 |
| In HGCA / Kervella+22 (Hipparcos stars only) | 29 |
| HGCA χ² 3σ | 5 |
| Kervella+22 PMa S/N > 3 | 8 |
| DR3 NSS solution | 5 |
| DR3 acceleration solution | 0 |
| DR3 RV-trend solution | 0 |

**HGCA 3σ and PMa signals**
- KELT-20, χ² = 62. This matches the companion of Yang+26 [32].
- TOI-4603, χ² = 62.
- HIP 94235, χ² = 59.
- TOI-201, χ² = 40.
- HD 118203, χ² = 12. This is the known companion c [34].
- PMa S/N > 3 only: HD 189733, DS Tuc A, TOI-1898 and HD 17156.
  - HD 189733 and DS Tuc A have known wide stellar companions, which could explain the anomaly.
  - **To check before treating any of these as an outer-planet hint.**

**DR3 NSS solutions**
- Kepler-16: AstroSpectroSB1 plus OrbitalTargetedSearchValidated, P = 41 d. This is the circumbinary host binary, as expected.
- WASP-2: OrbitalTargetedSearch, P = 38 d. A known close stellar companion; to check.
- HATS-26: OrbitalTargetedSearch, P = 194 d. To check.
- WASP-18: SB1, P = 0.94 d. This is the hot Jupiter's own RV signal [A7].
- TIC 172900988: SB2, P = 19.7 d. A circumbinary-planet host binary.
- "OrbitalTargetedSearch" solutions are from the DR3 exoplanet pipeline's targeted search of known systems [A6].

**What this tells us**
- Prior astrometric evidence is sparse: most of our hosts are too faint for Hipparcos.
- **The 52 hosts with RUWE > 1.4 are already suspected binaries.** Flag them so that a DR4 "planet" there is not over-interpreted.

## 7. Validation set

`data/known_outer_companions.csv` is curated, with a reference on every row. `scripts/validation_systems.py` → `data/validation_known_companions.csv`. Masses are M sin i unless marked as true masses, so the expected signals are lower limits there.

**Strong signals at P > 4 yr** (curvature or a partial orbit; joint RV fit needed):

| Companion | Edge-on Δχ² | Period |
|---|---|---|
| WASP-8 c | 8908 | 11.8 yr |
| HAT-P-2 c | 2423 | 23 yr |
| HD 118203 c | 1274 | 14 yr |
| WASP-53 c | 1195 | 7.8 yr |
| HAT-P-17 c | 442 | 13 yr |
| WASP-132 d | 145 | 5.0 yr |
| TOI-2537 c | 124 | 5.2 yr |
| K2-232 c | 123 | 5.3 yr |

**Astrometry-only orbit expected**
- WASP-81 c: 57 M_J M sin i (a brown dwarf), P = 3.6 yr, Δχ²_edge = 616. **This is our best end-to-end pipeline test.**

**Below threshold** (useful for false-alarm checks: we must *not* detect these):
- HAT-P-13 c: Δχ²_edge 17, Δχ²_face 33.
- HATS-59 c: 14.
- KELT-6 c: 8.
- HAT-P-44 c: 4.
- WASP-41 c: 3.
- TOI-2158 c: 0.9.
- WASP-47 c: 0.2.
- Kepler-424 c: 0.2.
- Kepler-419 c: 0.5.

**Still to add:** linear-trend-only hosts (HAT-P-4, -7, -22, -29, -32; WASP-10, -22, -34) [1]. They are acceleration-model tests. HAT-P-10 and XO-2 are not in our target list.

## 8. References

**RV, occurrence and theory**

| # | Reference |
|---|---|
| 1 | Knutson+2014, ApJ 785, 126, arXiv:1312.2954 |
| 2 | Bryan+2016, ApJ 821, 89, arXiv:1601.07595 |
| 3 | Ngo+2015, ApJ 800, 138, arXiv:1501.00013 |
| 5 | Ngo+2016, ApJ 827, 8, arXiv:1606.07102 |
| 7 | Zink & Howard 2023, ApJL, arXiv:2310.01567 |
| 8 | Shariat+2026, arXiv:2609.19249 |
| 9 | Wu, Rice & Wang 2023, AJ 165, 171, arXiv:2302.12778 |
| 11 | Huang, Wu & Triaud 2016, ApJ 825, 98, arXiv:1601.05095 |
| 12 | Dong, Katz & Socrates 2014, ApJL 781, L5, arXiv:1309.0011 |
| 14 | Petrovich & Tremaine 2016, ApJ 829, 132, arXiv:1604.00010 |
| 16 | Ranshaw+2026, AJ (accepted), arXiv:2609.09077 |
| 19 | Fulton+2021, ApJS, arXiv:2105.11584 |
| 21 | Rosenthal+2022, ApJS, arXiv:2112.03399 |
| 27 | Do Ó+2026, arXiv:2609.37160 |
| 28 | Masuda, Winn & Kawahara 2020, AJ, arXiv:1912.00387 |
| 29 | Becker+2017, AJ, arXiv:1710.01737 |
| 30 | Xiao+2026, ApJL (accepted), arXiv:2609.01377 |
| 31 | Xuan & Wyatt 2020, MNRAS 497, 2096, arXiv:2007.01871 |
| 32 | Yang, Lu, Yahalomi & Winn 2026, arXiv:2606.02705 |
| 33 | Damasso+2020, A&A 642, A31, arXiv:2007.06410 |
| 34 | Zhang+2024, AJ 168, 295, arXiv:2407.21377 |
| 35 | Dawson+2014, ApJ 791, 89, arXiv:1405.5229 |
| 37 | Petrovich 2015a, arXiv:1405.0280 |
| 38 | Anderson, Storch & Lai 2016, MNRAS, arXiv:1510.08918 |
| 39 | Petrovich 2015b, ApJ 805, 75, arXiv:1409.8296 |
| 41 | Dawson & Johnson 2018, ARA&A, arXiv:1801.06117 |

Also cited: De Rosa+2020, A&A 640, A73, arXiv:2007.08549.

**Astrometry**

| # | Reference |
|---|---|
| A1 | Brandt 2021, ApJS 254, 42, arXiv:2105.11662 |
| A3 | Kervella+2022, A&A 657, A7, arXiv:2109.10912 |
| A5 | Halbwachs+2023, A&A 674, A9, arXiv:2206.05726 |
| A6 | Holl+2023, A&A 674, A10, arXiv:2206.05439 |
| A7 | Gaia Collab., Arenou+2023, A&A 674, A34, arXiv:2206.05595 |
| A8 | ESA DR3 known issue, NSS false positives, 2024-05-27 (cosmos.esa.int) |
| A9 | Winn 2022, AJ, arXiv:2209.05516 |
| A10 | Marcussen & Albrecht 2023, AJ, arXiv:2305.08623 |
| A11 | Marcussen+2026, arXiv:2609.08590 |
| A12 | Stefansson+2025, AJ 169, 107, arXiv:2410.05654 |
| A13 | Sozzetti+2023, A&A 677, L15 (arXiv ID **UNVERIFIED**) |
| A16 | Sahlmann & Gómez 2025, MNRAS, arXiv:2404.09350 |
| A17 | Fitzmaurice+2023, arXiv:2310.07827 |
| A18 | An+2025, arXiv:2508.08374 |
| A19 | El-Badry+2024, OJAp, arXiv:2411.00088 |
| A25 | Feng+2022, ApJS 262, 21, arXiv:2208.12720 |
| A26 | Thompson+2026, arXiv:2602.00235 |
| A29 | de Beurs+2023, AJ 166, 136, arXiv:2309.03256 |
| A33 | Belokurov+2020, MNRAS, arXiv:2003.05467 |
| A34 | Ziegler+2020, AJ, arXiv:1908.10871 |
| A35 | Ziegler+2021, AJ, arXiv:2103.12076 |
| A36 | Ziegler+2026, arXiv:2608.26250 |
| A41 | Michel & Mugrauer 2024, MNRAS 527, 3183 |

**Forecasts**

| # | Reference |
|---|---|
| L1 | Lammers & Winn 2025, AJ (accepted), arXiv:2511.04673 |
| E1 | Espinoza-Retamal, Zhu & Petrovich 2023, AJ 166, 231, arXiv:2309.08665 |
| E2 | Greenfield & Owen 2026, arXiv:2610.07149 |
| E3 | Perryman+2014, ApJ 797, 14, arXiv:1411.1173 |

**Sources for the validation table**
- Neveu-VanMalle+2016, A&A 586, A93, arXiv:1509.07750
- Endl+2014, ApJ 795, 151, arXiv:1409.8156
- Triaud+2017, MNRAS 467, 1714, arXiv:1612.04166
- Knudstrup+2026, arXiv:2609.30915

## Open questions for Phase 1

- Add a joint astrometry + RV mode for known RV companions and trends, or scope it out explicitly. Most validation systems need it.
- Sample inclination explicitly in injection–recovery. Aligned (edge-on) companions are the hardest, and they are the science case.
- Decide the threshold: Δχ² > 50 (Lammers) versus our own null-calibrated false-alarm probability (FAP). With ~790 hosts, a 5% per-star false-alarm rate would give about 40 false detections, so we need FAP ≪ 1/790 per star.
