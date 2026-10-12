# Euclid × Gaia brown-dwarf companion pilot: design

Written 2026-10-11. Design only; nothing has been run.

**Goal.** Build the cross-match pipeline on public Euclid Q1 data and Gaia DR3 now, so that on release it is a data swap to Euclid DR1 and Gaia DR4. **Priority rule** (CLAUDE.md §2b): build only if Gaia Phase 1 is on track in late November; run after the Gaia release-day work.

**Target regime** (`notes/euclid_overlap.md` §5):
- late-L to late-T companions, at 1″ to ~10,000 au projected separation;
- around Gaia stars within ~200 pc;
- with a calibrated chance-alignment probability for every candidate;
- delivering both candidates and a companion fraction with limits.

## 1. Inputs

| | Now (pilot) | At release |
|---|---|---|
| Euclid | Q1, `astroquery.esa.euclid`, environment `PDR`, `catalogue.mer_catalogue` (+ `catalogue.mer_morphology`, spectra) | DR1-Foundation MER catalogue (environment and table names to verify) |
| Gaia hosts | DR3, `gaiadr3.gaia_source` | DR4 `gaia_source` (better parallaxes and proper motions; epoch astrometry to vet binary hosts) |
| Earlier epoch | CatWISE2020, Legacy Surveys DR10, VHS DR6/DR7, UKIDSS DR11+, Pan-STARRS1, 2MASS (VizieR / NOIRLab TAP / ESO / WSA) | same |
| Calibrators | Q1 UCDs with spectral types (Žerjal+2026, Dominguez-Tagle+2025, Mohandasan+2025 catalogues); UCDs with Gaia parallaxes in Q1 | DR1 equivalents |

**Host selection**
- Parallax > 5 mas, parallax_over_error > 10, inside the Euclid footprint (MER tile cores or the DR1 MOC).
- Keep all RUWE values but flag RUWE > 1.4. Unresolved binary hosts are interesting but complicate ages.
- Q1 host counts are in `results/euclid_q1_hosts.csv` (`scripts/euclid_q1_hosts.py`).

## 2. Candidate selection, per host

1. **Search annulus.**
   - Outer radius: projected 10,000 au, i.e. 10,000 × ϖ[″] (10 pc → 1,000″; 200 pc → 50″). Cap it at 5′ to keep the chance alignment tractable.
   - Inner radius: set by the host's PSF wings and saturation, an empirical function of host G and I_E measured on Q1 (the radius where spurious_prob and point_like_prob degrade). Expect ≥ 1–3″ for G = 10–14 hosts and more for brighter ones.
2. **Clean point sources:** point_like_prob high, spurious_flag = 0, det_quality_flag clean, and NISP detection (Y_E, J_E, H_E at S/N ≥ 5).
3. **Brown-dwarf colours:**
   - Locus in I_E − Y_E, Y_E − J_E and J_E − H_E.
   - L/early T: I_E − Y_E > 2.5 (Žerjal+2026).
   - Mid/late T: methane-blue J_E − H_E with red I_E − J_E (Kiwy+2025 colour space).
   - Galaxies and QSOs are rejected via point-likeness and the QSO locus (Mohandasan+2025 §3).
   - Calibrate the locus and its width on the spectroscopically typed Q1 UCDs.
4. **Photometric distance.**
   - Assign a spectral type from colours by template matching (or from NISP spectra when available).
   - Get M_JE from an absolute-magnitude–spectral-type relation built from Q1 UCDs with Gaia or literature parallaxes. Fallback: MKO relations (Dupuy & Liu 2012; Kirkpatrick+2021) transformed to the Euclid bands.
   - Then d_phot with a scatter of σ(M) ≈ 0.3–0.5 mag.
   - Require |μ_phot − μ_host| < 2.5 σ_tot in distance modulus (σ_tot includes the host parallax error and the relation scatter).
5. **Proper motion** where an earlier epoch detects the source (roughly J_Vega ≲ 20–21):
   - Position offset between epochs over the baseline gives μ.
   - Require consistency with the Gaia host proper motion within 3σ. The host's motion over a 5–15 yr baseline is usually ≫ the measurement error for hosts within 200 pc.
   - Use the earlier-survey astrometric precision and Euclid MER position errors (`notes/euclid_scoping.md`).
   - A non-detection at the predicted earlier-epoch position is not informative for faint candidates. Mark them "photometric only".
6. **NISP spectra** for H_E ≲ 21–22.5: spectral type and classification from the combined 1D spectra. Spectral and photometric types must agree.
7. **Chance-alignment probability, for every candidate.**
   - Count unrelated sources passing steps 2–4 (and 5 where applicable) around shifted host positions: El-Badry-style mock hosts, displaced by 0.5–2° at the same Galactic latitude, inside the footprint.
   - That gives the expected number N_exp of interlopers within the candidate's separation.
   - P_chance = 1 − exp(−N_exp). Report N_exp per candidate.
   - Use the same calibration for the population: expected background counts are subtracted when computing companion fractions.
   - Candidates with P_chance < 1% are "probable"; < 0.1% with common proper motion are "confirmed".

## 3. Outputs

**Per-candidate table**
- host Gaia id, separation (″ and au), Euclid id, I/Y/J/H, colour type, spectral type (spectrum);
- d_phot vs d_host, proper motion and its χ² against the host (or "none");
- N_exp and P_chance, flags.

**Population**
- Wide-companion fraction per spectral-type bin versus separation, with completeness from injected synthetic brown dwarfs in Euclid images or catalogues.
- Completeness needs care: either image-level injection (MER pipeline unavailable to us) or catalogue-level, using the depth maps (`q1.vmpz_healpix_depthmap`).

## 4. Validation

- **Known benchmark companions in Q1 fields.** Query SIMBAD and the literature for L/T companions to stars within the Q1 tiles (EDF-N, EDF-S, EDF-F, LDN 1641). Expect very few or none: Q1 is only 63 deg² at high latitude. Each one recovered tests steps 2–7.
- **The 2 Mohandasan+2025 candidate companions** (EDF-N; §4.4, Table 4). Our pipeline should recover them and assign P_chance. They flag one as a probable false positive, which is a direct check of our chance-alignment estimate.
- **Injection tests:** synthetic L/T dwarfs at host distances at catalogue level (colour, magnitude, noise from the MER error model) → completeness against spectral type and separation.
- **Null tests:** run the full selection around shifted (fake) host positions. The candidate rate must match N_exp.

## 5. Code plan (when green-lit)

| Module | Contents |
|---|---|
| `companions/euclid/archive.py` | Euclid and earlier-survey queries; cached raw downloads in `data/raw/euclid/` |
| `companions/euclid/select.py` | annulus, cleaning, colour locus, photometric distance, proper motion |
| `companions/euclid/chance.py` | shifted-host chance-alignment calibration |
| `scripts/euclid_pilot.py` | Q1 + DR3 end to end. A config switch (environment and table names, Gaia release) makes DR1 + DR4 a data swap |

Tests use the TDD style of Phase 1: synthetic catalogues for selection and chance alignment, plus the Mohandasan candidates as a real-data check.

## 6. Risks

- **The EC proprietary head start** and Backyard Worlds speed. Mitigate with the niche (faint, close, statistical; exoplanet and white-dwarf hosts) and by preparing now.
- **Proper motions exist only for the bright subset.** For faint companions, association is statistical (photometric distance and P_chance). Say so honestly.
- **The bright-host halo** limits small separations around G < 10 hosts.
- **DR1 environment and table names, and the DR4 source-id mapping,** are unknown until release.
