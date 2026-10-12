# Euclid DR1 scoping for brown-dwarf companions to Gaia stars

Written 2026-10-11 for Phase 2 (CLAUDE.md §2b).

**Sources:** ESA pages and documents (URLs below), verified here unless marked otherwise, plus our own archive tests. Agent research was spot-checked against the primary PDFs.

## 1. What is released when

**The November tranche contains what we need.**

**Release dates.** DR1-Foundation is on **12 November 2026**. Full DR1 follows in mid-2027.

**Primary source:** ESA Science Team memo EUCL-EST-ME-8-020, *Content of the ESA Euclid Data Release 1*, issue 1.3, approved 22 Sep 2026 (`https://www.cosmos.esa.int/documents/10647/12245842/DR1-memo_260922.pdf`). Read locally:
- §1: "On 12th November 2026, Euclid will release all raw data and calibrated images, catalogues, spectra, covering a large sky area of about 1900 deg² (DR1-Foundation)."
- §3: "DR1-Foundation will include all products in section 3, except all products in sub-section 3.2.7; DR1 will include all products in section 4, and those in section 3.2.7; there will be no reprocessing of products delivered at DR1-Foundation."

**Section numbers** refer to the archive product list EUCL-ESAC-RP-8-017 v4.2 (30 Jul 2026), checked locally:

| Section | Products |
|---|---|
| 3.2.1 VIS, 3.2.2 NIR | calibrated and stacked frames |
| 3.2.3 SIR | science frames and **combined 1D spectra** |
| 3.2.4 MER | mosaics, segmentation maps and the **final source catalogue** (plus morphology and cutouts) |
| 3.2.5 PHZ | photo-z, including the star SED catalogue |
| 3.2.6 SPE | spectroscopic redshift, lines and models catalogues |
| 3.2.8 SLE | strong-lens detections |
| 3.3 | VMPZ masks |
| **3.2.7 SHE (shear)** | the only section 3 product deferred, with the section 4 Level-3 weak-lensing, clustering and cluster products, to mid-2027 |

**Footprint**
- About 1,900 deg²: roughly 500 deg² in the north and 1,400 deg² in the south (memo §2).
- Wide survey observed 14 Feb 2024 – 21 May 2025. Ice-affected data from 23 Mar – 7 Jun 2024 were replaced by "INFILL" re-observations in Mar–May 2025 (memo §3; per-patch dates in its Appendix Table 1, which gives the first-epoch times we need).
- An expected-coverage HEALPix map exists (`cgv_map_dr1input_o13.fits.gz` on `cosmos.esa.int/web/euclid/euclid-dr1-expected-coverage`, ~2,100 deg² of input). Inclusion is not guaranteed.
- Deep fields in DR1: EDF-N with 8 red + 3 blue grism visits, EDF-S with 10 red, EDF-F with 1 + 1. **Repeated deep-field visits give multi-epoch astrometry in the deep fields** (our inference; not a release product).

**No proper motions in DR1.** No official multi-epoch astrometric product exists (memo and product list); the agent found no statement either. A second wide epoch exists only as a proposal (Bedin 2025, arXiv:2510.23694).

## 2. Access (tested 2026-10-11)

**What works**
- `astroquery.esa.euclid` (astroquery 0.4.11) defaults to the public environment `PDR` and needs **no login**.
- `Euclid.load_tables()` lists 103 tables, including `catalogue.mer_catalogue` (471 columns), `catalogue.mer_morphology`, `catalogue.phz_*`, `catalogue.spectro_*`, `q1.mosaic_product`, `q1.spectra_source` and `q1.vmpz_healpix_*`.
- Relevant `mer_catalogue` columns:
  - positions: `right_ascension`, `declination`, plus `*_psf_fitting`;
  - point-source class: `point_like_prob`, `point_like_flag` (prob > 0.96 + quality), `spurious_flag`, `det_quality_flag`;
  - photometry: `flux_vis_psf`, `flux_{y,j,h}_templfit` and aperture fluxes (µJy; AB = −2.5 log₁₀ F + 23.9).
- Synchronous ADQL cone counts on `catalogue.mer_catalogue` work.

**What doesn't (yet)**
- **TAP table upload with a server-side join failed anonymously.** The astroquery docs say upload and cross-match "will also be available" with DR1. For now, either query per host (sync cones) or log in.

**Other methods (astroquery docs, not yet exercised):** `cone_search`, `get_cutout` (MER mosaics), `get_spectrum` (`SPECTRA_RGS`), `get_product`, `query_sia`.

**The DR1 environment and table names are not yet known.** Q1 uses `catalogue.*` plus `q1.*`; DR1 will likely be similar (UNVERIFIED).

## 3. Depth and distance reach

**Wide Survey median 5σ point-source depth (AB)** (euclid-ec.org science overview): I_E 26.2, Y_E 24.3, J_E 24.5, H_E 24.4. Deep fields are about 2 mag deeper.

**Measured on Q1**
- VIS PSF FWHM 0.16″.
- Absolute astrometric accuracy relative to Gaia: 5 mas (arXiv:2503.15303).
- MER tile offsets: RMS 2–6 mas (arXiv:2503.15305).

**Distance reach for companions** (rough; scoping only):
- d = 10^((m_lim − M)/5 + 1), with M_J (MKO, Vega) from Dupuy & Liu 2012 / Kirkpatrick+2021.
- J_AB ≈ J_Vega + 0.91, and assumes J_E ≈ J_MKO,AB. That carries a ±0.3 mag systematic until calibrated on the Q1 UCDs.
- Late T and Y dwarfs are invisible in VIS (I_E − J_E ≳ 5), so the selection must rest on NISP Y/J/H.

| Type | M_J (Vega) | Euclid J_E 5σ (24.5) | Euclid J_E 10σ (~23.75) | VHS J ≈ 20.2 | CatWISE (W2 ≈ 17.5) |
|---|---|---|---|---|---|
| L5 | ~13.3 | ~1,100 pc | ~800 pc | ~240 pc | — |
| T5 | ~14.9 | ~550 pc | ~390 pc | ~110 pc | ~90 pc |
| T8 | ~16.3 | ~290 pc | ~200 pc | ~60 pc | ~50 pc |
| Y0 | ~19.5 | ~80 pc | ~60 pc | ~15 pc | ~20 pc |

So Euclid alone sees late-T companions out to ~200 pc and Y dwarfs to ~60 pc. That is a factor of about 4–5 in distance (about 100× in volume) beyond VHS and CatWISE.

**Resolution:** 0.16″ VIS and ~0.5″ NISP give access to 1–10″ separations (100–2,000 au at 100–200 pc), where WISE (6″) cannot resolve. The practical inner limit is set by the bright-host halo (pilot design §2).

## 4. Proper motions: earlier epochs and precision

**Euclid epochs:** 2024–2025 (DR1 wide survey; Q1 July–September 2024).

**Euclid centroid precision** ≈ FWHM / (2.35 × S/N):
- VIS: ~7 mas at S/N = 10.
- NISP: ~20 mas for VIS-invisible T dwarfs.

So the first epoch dominates the error budget.

| Earlier survey | Epochs | Depth (relevant band) | Overlap with DR1 | Baseline | Typical μ precision for faint companions* |
|---|---|---|---|---|---|
| 2MASS | 1997–2001 | Ks ≈ 14.3, J ≈ 15.8 (Vega) | all sky | ~25 yr | ~4 mas/yr, but only bright L/early T near the host |
| UKIDSS LAS | 2005–2012 | J ≈ 19.6, K ≈ 18.2 | north equatorial; partial overlap with northern DR1 | 12–19 yr | ~5–8 mas/yr |
| VHS | 2009–2017+ | J ≈ 20.2, Ks ≈ 18.1 | southern sky; large overlap with southern DR1 | 7–15 yr | ~5–15 mas/yr |
| Pan-STARRS1 3π | 2010–2014 | z ≈ 21.5, y ≈ 20.3 (AB) | δ > −30°; partial | ~11–14 yr | L dwarfs only (T too faint in y) |
| Legacy Surveys DR10 (DECaLS/DES/BASS/MzLS) | 2013–2019 | z ≈ 22.5 (AB) | extensive in both hemispheres | 5–11 yr | ~10–20 mas/yr; L and early T |
| CatWISE2020 | 2010–2018 | W1 ≈ 17.7, W2 ≈ 17.5 | all sky | — | its own μ: ~20 mas/yr at W1 = 15, ~100 at W1 = 17; 6″ PSF |
| Euclid deep fields (own repeats) | 2023–2025 visits | 2 mag deeper | EDFs only | 1–2 yr | ~10–20 mas/yr |

\*Precision ≈ first-epoch positional error / baseline. These are scoping estimates from survey specifications and typical faint-source astrometry (~50–150 mas). Measure them in the pilot.

**Host motions:** for hosts within 200 pc, median |μ| ≈ v_tan / (4.74 d) ≈ 30 mas/yr at 200 pc and ≈ 125 mas/yr at 50 pc. So common proper motion is decisive within ~100 pc but only a 2–5σ test at 200 pc. **Companions fainter than about J_Vega ≈ 20–21 have no earlier-epoch detection**; their association rests on photometric distance plus chance-alignment statistics (pilot design §2.7).

## 5. Gaia hosts in the footprint

Source: `scripts/euclid_q1_hosts.py` → `results/euclid_q1_hosts.csv`.

**Method:**
- Q1 VIS MER tiles: 352 tiles in 4 fields (EDF-S, EDF-N, EDF-F, LDN 1641).
- Gaia DR3 stars with parallax > 5 mas and parallax_over_error > 5 in the 30′ tile cores.
- Corrected by the covered fraction (78% of 120 random hosts have MER sources within 30″), giving an effective area of ~69 deg² (published Q1: 63.1 deg² + LDN 1641).

| Distance | Q1 (covered) | Density | DR1 (1,900 deg²) estimate |
|---|---|---|---|
| < 50 pc | ~45 | 0.66 deg⁻² | ~1,250 |
| < 100 pc | ~336 | 4.9 deg⁻² | ~9,300 |
| < 200 pc | ~2,380 | 34.5 deg⁻² | **~65,000** |

**Caveats**
- Q1 fields are at high Galactic latitude; DR1's wide patches span a broader latitude range, so this is approximate (likely a mild underestimate if DR1 includes lower-latitude patches).
- About 75–90% of these hosts are M dwarfs, as for any volume-limited sample.

**Expected yield** (order of magnitude only):
- Marocco+2017 forecast ~24,000 resolved UCD companions to Gaia primaries at |b| > 10°. Scaling by area (1,900 / ~34,000 deg²) gives ~1,300 for DR1, mostly M7–L.
- Wide L/T companion fractions are ~1–4% (dal Ponte+2020).
- So **tens to a few hundred** late-L/T companions to hosts within 200 pc, before chance-alignment cleaning and completeness.

## 6. Answers to the brief

1. **November tranche:** **yes.** It contains the merged MER catalogue, NISP 1D spectra and SPE redshift catalogues, on 12 Nov 2026 (memo §1, §3; product list §3).
2. **Access:** works anonymously on Q1 (ADQL, catalogue columns, tables listed). Table upload needs a login or the DR1 upgrade.
3. **Depth:** late T to ~200 pc and Y0 to ~60 pc at 10σ in J_E; 4–5× the distance of VHS and CatWISE.
4. **Proper motions:** none from Euclid alone. Earlier epochs (VHS south, UKIDSS north, Legacy Surveys, Pan-STARRS, 2MASS, CatWISE) give ~5–20 mas/yr for J_Vega ≲ 20–21. Fainter companions are statistical.
5. **Hosts:** ~65,000 Gaia stars within 200 pc in DR1 (~9,300 within 100 pc).

## 7. Open items before the pilot

- **Done 2026-10-12** (`scripts/euclid_dr1_footprint.py` → `results/euclid_dr1_hosts.csv`):
  - The ESA expected-coverage map (`cgv_map_dr1input_o13`, nside 8192) covers **2,109 deg²** of DR1 input. ESA expects ~1,900 deg² in the release; inclusion is not guaranteed.
  - The footprint's median Galactic latitude is |b| = 43° (10–90%: 33–60°), similar to the Q1 fields, so the §5 host-density scaling is reasonable.
  - **41 of our 788 transiting giant-planet hosts** fall in the footprint (5.2%, consistent with the area fraction). Examples: KELT-24 (97 pc), TOI-201 (112 pc; HGCA 3σ), TOI-2180 (116 pc; 261-d warm Jupiter), TOI-2447, WASP-97.
  - These are the Phase 1 × Phase 2 overlap: Gaia DR4 astrometric limits on close companions plus Euclid limits on wide brown-dwarf companions for the same hosts.
- Calibrate J_E ↔ J_MKO and the absolute-magnitude relations on Q1 UCDs with Gaia parallaxes.
- Measure the real first-epoch astrometric precision (VHS, Legacy Surveys) against Euclid Q1 for UCDs.
- Confirm the DR1 archive environment and table names, and whether table uploads become available, on 12 Nov.
