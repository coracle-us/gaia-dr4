# Gaia DR4 data model: epoch astrometry and epoch RVs

Written 2026-10-11 (Phase 0, step 1). Main source: ESA's **draft DR4 data model PDF**
(pre-released June 2026). Checked against the pre-release epoch-astrometry VOTable and
ESA's `gaiasupdate` code. The model is a draft, so re-check field names on release day.

## Sources

| What | Where |
|---|---|
| DR4 overview + content table (rows/sizes per table, updated 14 Sep 2026) | https://www.cosmos.esa.int/web/gaia/dr4 |
| Previews index | https://www.cosmos.esa.int/web/gaia/dr4-previews |
| Pre-release page (sample, data model, tools) | https://www.cosmos.esa.int/web/gaia/dr4-prerelease |
| Draft data model (zip → PDF) | https://anonftp.cosmos.esa.int/pub/GAIA_PUBLIC_DATA/Gaia_DR4/dr4-prerelease/gaia-dr4-prerelease-draft-data-model_2026-06-26.zip |
| Pre-release epoch astrometry, 12 sources | https://anonftp.cosmos.esa.int/pub/GAIA_PUBLIC_DATA/Gaia_DR4/dr4-prerelease/gaia-dr4-prerelease-epoch-astrometry_2026-06-26.zip |
| DR4 software tools | https://www.cosmos.esa.int/web/gaia/dr4-software-tools |
| `gaiasupdate` (ESA single-star fitter) | https://github.com/esa/gaia-supdate · `pip install gaiasupdate` (v0.1.2, 2026-06-26) · docs https://esa.github.io/gaia-supdate/ |
| ESA DR4 tutorial notebooks | https://github.com/esa/gaia-jupyter-notebooks/tree/main/data-release-4-tutorials |
| DR3 DataLink reference (formats, limits) | https://www.cosmos.esa.int/web/gaia-users/archive/datalink-products |

Local copies are in `data/raw/dr4_prerelease/` (gitignored). SHA-256 checksums:
- `d807eae9…f42c66`: data model zip
- `07f0e8d9…5fa0b`: epoch astrometry zip

Full hashes are in `data/raw/dr4_prerelease/SHA256SUMS.txt`.

## Headline facts

- **Coverage:** 25 Jul 2014 to 20 Jan 2020 (66 months). 2.8 billion sources were processed, about 2 billion of them in `gaia_source`.
- **Astrometric solution:** reference epoch **J2017.5 (TCB)**, frame Gaia-CRF4. The core model is a 6-parameter solution: 5 astrometric parameters plus a "constrained colour" term (pseudo-colour).
- **`epoch_astrometry`:** 129.6 billion rows, 62 TB, **DataLink only, not TAP** (data model §3.1).
- **`rvs_epoch_parameters_single`:** 6.9 billion rows, 323 GB, **DataLink only**. Epoch RVs are computed only for **G_RVS ≤ 14** (§17.2). For fainter stars the table holds only metadata.
- **Availability flags (TAP tables):**
  - `all_source_flags.has_epoch_astrometry`
  - `all_source_rvs.has_epoch_rvs` together with `fraction_double_lined` (0 or null means single-lined)
- **NSS outputs (TAP):**
  - `nss_two_body_orbit`: 5.9 million rows. Solution types include `Orbital`, `OrbitalPoorlyConstrained` (significance < 5), `OrbitalAlternative`, `AstroSpectroSB1/2`, `SB1/2`, `VIMO`, …
  - `nss_multiple_orbits`: 16.7 thousand rows
  - `nss_acceleration_astro`: 3.4 million rows
  - `nss_masses`: 4.3 million rows. It has a field `exoplanet_candidate_hrd_proba`: the probability that the secondary is below 20 M_J, given the primary's position on the HR diagram.
  - `nss_epoch_flags` (DataLink): which epochs the NSS processing used
- **No separate exoplanet table** appears in the content list. DPAC planet candidates are expected inside the NSS tables, plus an "exoplanet list" auxiliary product and the papers *The Exoplanet Candidate Catalogue*, *Astrometric orbit determination for exoplanet discovery* and *Exoplanet orbit processing with radial velocity constraints*. See `dr4_paper_overlap.md`.
- **Scanning law:** published in DR4 as `commanded_scan_law` (TAP, "including Ecliptic Pole Scanning") and `gaia_gap` (excluded intervals). This is useful for Phase 1 injection–recovery. No GOST-type tool is listed yet.

## `epoch_astrometry` (DataLink retrieval type `EPOCH_ASTROMETRY`)

**Row structure.** There is one row per (source_id, transit_id), i.e. per **FoV transit**. Per-CCD quantities are **arrays of length 10**:
- Index 0 is the SkyMapper CCD (SM1/SM2). Indices 1–9 are AF1–AF9.
- Per `gaiasupdate`, SM `centroid_pos_error_al` "is significantly underestimated, especially for bright stars" (data model). `gaiasupdate` filters SM out, and so will we.

**What the data are.** These are measured local-plane coordinates relative to a fixed per-source tangent point (`ra0`, `dec0`). They are **not residuals**. Only the **AL coordinate** is measured in DR4: the AC centroid field and AC usage are empty or false. Model fitting is therefore the 1-D along-scan formulation, as planned in CLAUDE.md.

Fields (37 columns in the pre-release VOTable, release tag "Gaia DR4_RC3"; `[10]` marks a per-CCD array):

| Field | Type / unit | Meaning, and our use |
|---|---|---|
| `solution_id`, `source_id`, `transit_id` | long | IDs. `transit_id` encodes OBMT `= 204800*((tid>>17)&0x3FFFFFFFFFF)` ns, FoV `= 1+((tid>>15)&3)`, CCD row `= (tid>>12)&7` |
| `ra0`, `dec0` | double, deg | Tangent point that w and z are counted from. Constant per source |
| `agis_source_excess_noise` | float, mas | **Unscaled** AGIS source excess noise (differs from the rescaled value in gaia_source) |
| `obs_time_tcb` [10] | long, **ns since J2010.0 (TCB)** | Centre of the exposure per CCD. Absolute time-transformation accuracy is about 150 ns |
| `obs_time_bary_corr` | float, ns | Rømer correction, TCB(barycentre) − TCB(Gaia), computed at AF4. At most about ±370 s |
| `scan_pos_angle` [10] | double, deg | Position angle of scan ψ. sin ψ = ∂w/∂α\*, cos ψ = ∂w/∂δ |
| `zeta` | float, deg | AC field angle of the reference point at AF4. Needed only to move the reference point (resolved pairs) |
| `parallax_factor_al`, `parallax_factor_ac` | float | ∂w/∂ϖ (∂z/∂ϖ), one per transit at AF4 |
| `colour_factor_al` [10] | float, mas·µm | −∂w/∂ν_eff. The pseudo-colour term of the 6-parameter model |
| `colour_factor_ac` [10] | float | Empty in DR4 |
| `nu_eff_used_in_astrometry`, `nu_eff_error` | float, µm⁻¹ | Colour assumed in IPD/AGIS. Default 1.43, error fixed at 0.085 |
| `centroid_pos_al` [10] | double, mas | **The measurement**: AL local-plane coordinate w |
| `centroid_pos_ac` [10], `centroid_pos_error_ac` [10] | double / float, mas | Not available in DR4 |
| `calculated_pos_ac` [10] | double, mas | AC coordinate computed from the source model and attitude, not measured |
| `centroid_pos_error_al` [10] | float, mas | σ² = σ²_IPD + σ²_calibration-excess (does **not** include source excess noise) |
| `ipd_error_al` [10], `ipd_error_ac` [10] | float, mas | IPD-only centroid error (AC is NaN for 1-D windows) |
| `used_by_agis_al` [10] | bool | CCD used in the AL source update. AGIS DR4 downweights are 0 or 1 |
| `used_by_agis_ac` [10] | bool | Always false in DR4 |
| `transit_acq_flags`, `transit_proc_flags` | short | Transit-level bit flags (acquisition, window truncation; processing) |
| `ccd_proc_flags` [10] | short | CCD-level bits: 0x000F IPD problem, 0x0010 reduced window, 0x0020 IPD unavailable, 0x0040 cosmetic, 0x0080 cosmic removed, 0x0100 saturation removed, 0x0200 partial window, 0x0400 no window, 0x0800 InitialCentroid failed, 0x1000 odd background, 0x2000 bias-mitigation fallback, 0x4000 non-target source removed, 0x8000 LSF/PSF clamped. Bit meanings are from `gaiasupdate/utils.py` |
| `multipeak`, `blended` | bool | Multiple peaks in the window; peak also matched to another source. **Impostor and contamination flags for us** |
| `g_mag`, `g_class` | float / byte | Transit G magnitude; window class |
| `gates` [10], `source_dist_to_last_ci` [10], `ac_rate`, `sub_pixel_coord` [10], `mu` [10] | | Instrument context: gating, charge injection, AC rate, AL sub-pixel phase, AC pixel |

**Weights.** To reproduce AGIS, use σ²_tot = `centroid_pos_error_al`² + `agis_source_excess_noise`², restricted to `used_by_agis_al`. **For orbit fitting, drop `agis_source_excess_noise`.** It is the noise AGIS added because the single-star model did not fit, so it soaks up exactly the signal we are looking for. In the pre-release sample it is 0.12 mas for Gaia-4 (the planet host) and 6.6 mas for Gaia BH3, but 0 for most single stars. We should fit our own jitter term instead.

**Known issue (data model).** About 16.4 million CCD values of `centroid_pos_al` fall outside the window, and 208 transits exceed 10⁹ mas, from attitude excursions at data-gap edges. Clip them.

## `rvs_epoch_parameters_single` (DataLink retrieval type `EPOCH_PARAMETERS_RVS_SINGLE`)

- One row per RVS FoV transit for single-lined sources. Double-lined transits are in `rvs_epoch_parameters_double`.
- Epoch RVs exist only where G_RVS ≤ 14 (`external_apparent_grvs` in gaia_source / all_source_rvs).
- Columns (data model §17.2):
  - IDs: `solution_id`, `source_id`, `transit_id`.
  - **`obs_time_rv`**: double, **barycentric JD (TCB) − 2455197.5** (day). Note this is not the ns convention of epoch astrometry.
  - `obs_time_bary_corr` (long, ns).
  - `radial_velocity`, `radial_velocity_error`, `raw_radial_velocity_error`, `rv_systematic_corr`, `rad_vel_bary_corr` (km/s).
  - `radial_velocity_status`.
  - `rv_validity_flag`: only rows with this set are valid.
  - Outlier fields: `rv_outlier_probability`, `rv_outlier_condition`, `outlier_status`.
  - `rv_assumed_sb2`.
  - Emission lines: `emission_lines`, `els_probability`[3].
  - Line broadening: `vbroad`, `vbroad_error`, `vbroad_validity_flag`, `vbroad_invalid_reason`.
  - RVS photometry: `grvs_flux`, `grvs_flux_error`, `grvs_mag`, `grvs_validity_flag`, `grvs_invalid_reason`.
  - `expected_sig_to_noise`.
- **Our use: impostor check.** A 5 M_J planet at 3 AU induces K ≈ 80 m/s. Gaia RVS epoch precision for bright stars is about 1 km/s or worse (to verify in DR4 docs), so epoch RVs cannot confirm a planet. They will readily exclude stellar-mass companions, where K is km/s, and SB2s.

## Access on release day (2 Dec 2026; archive offline from 25 Nov 12:00 CET)

1. **TAP (sync or async ADQL):**
   - `gaia_source` / `all_source_*` for our ~790 source_ids. Map DR3→DR4 IDs; a crossmatch table is expected, so check on release.
   - Read `has_epoch_astrometry` and `has_epoch_rvs`.
   - Read `nss_two_body_orbit`, `nss_acceleration_astro`, `nss_masses` and `nss_multiple_orbits` for every target. DPAC solutions must be cited.
2. **DataLink.** astroquery `Gaia.load_data(ids=[...], data_release='Gaia DR4', retrieval_type='EPOCH_ASTROMETRY', data_structure='RAW', format='votable')`, and the same with `EPOCH_PARAMETERS_RVS_SINGLE`.
   - `gaiasupdate.from_gacs_datalink` already wraps the astrometry call (internal release string `'Gaia DR4_INT4'` for now). The public `data_release` string is unknown until release.
   - Data structures: `DATAMODEL_GAIA` (described above) vs `DATAMODEL_STANDARD`. The pre-release file is the "RAW" Gaia structure.
   - DR3 DataLink allowed **5,000 sources per request**. Our ~800 hosts fit in one or a few batched calls. Expect load on day 1: batch, retry and cache everything under `data/raw/dr4/`.
3. Bulk or CDN download is mentioned by `gaiasupdate` ("via DataLink or bulk download"). It is not needed for a few hundred stars.

## Software status (2026-10-11)

- **`gaiasupdate` is released** (ESA, PyPI 0.1.2, ESA-PL Permissive licence; DOI 10.5281/zenodo.20924676). CLAUDE.md said "announced"; that is out of date.
  - Models: `3p_single_source_without_offsets`, `5p_single_source`, `6p_constrained_colour` (DR4-like, the default) and `6p_perspective_acceleration`.
  - Outputs: parameters, covariances, residuals, outliers, χ², BIC/AIC, F2 and excess noise.
  - **It has no orbit model.** Installed in `.venv` and listed in `requirements.txt`.
  - On Gaia-4 it gives ϖ = 13.624 mas, μ = (−75.552, 17.944) mas/yr from 824 AL CCD measurements (`scripts/plot_prerelease.py`).
- **Other DR4-ready orbit code exists:**
  - ESA `gaia-bhthree` (DR4 pre-release branch): orbit fit for BH3.
  - **`exogaia`** (T. Stolker; PyPI 0.1.3, 2026-09-28; https://github.com/tomasstolker/exogaia): a 12-parameter orbit fit with least squares plus dynesty that already reads the pre-release format.
  - El-Badry's `gaiamock`: forward model and fitting cascade, with GOST-based epochs.
  - Phase 1 should cross-check our fitter against `exogaia` and `gaiamock` as well as `gaiasupdate`.
- ESA says further tools will ship with the data on release day.

## Pre-release sample (`data/raw/dr4_prerelease/epoch_astrometry/`)

- **Contents:** 12 sources, 1,008 FoV-transit rows, astrometry only (**no epoch RVs**). Groups: parallax (3), magnitude (3), QSO (3), and orbit: HD 114762, Gaia-4 and Gaia BH3.
- **Plot:** `figures/prerelease_1457486023639239296.png` shows Gaia-4. Residuals from the 6-parameter fit, binned per transit, have **0.12 mas RMS against 0.05 mas errors**, with a coherent signal of about 1.5 years, as expected from the known giant planet Gaia-4b.
- **Noise-model check:** the median per-CCD `centroid_pos_error_al` over used CCDs is:

  | Source | Per-CCD error (mas) | Holl+2023 curve in `build_targets.py` (mas) |
  |---|---|---|
  | G = 7.2 | 0.125 | 0.17 |
  | G = 9.0 | 0.105 | 0.13 |
  | G = 11.9 | 0.082 | 0.125 |
  | G ≈ 14 | 0.14–0.16 | 0.15 |

  The target-list sensitivity is therefore slightly **conservative** at the bright end. Recalibrate it in Phase 1 against the real DR4 errors.

## Open items to check on release day

- The public `data_release` string for DataLink, and the DR3→DR4 source_id crossmatch.
- CCD-to-CCD correlations within a transit and the attitude-noise floor: the astrometry paper (*Gaia DR4: The astrometric solution*) should cover them. `exogaia` assumes CCDs are uncorrelated; ESA does not say.
- How DPAC planet candidates are labelled in `nss_two_body_orbit` / `nss_masses`, and what is in the "exoplanet list" auxiliary product.
- The RVS epoch precision versus G_RVS (papers *Epoch RVS spectra and radial velocities* and *Properties and validation of the radial velocities*).
