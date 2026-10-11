# Overlap check: DR4 collaboration papers and other teams

Written 2026-10-11 (Phase 0, step 3).

Our question: outer giant companions to transiting hot/warm-Jupiter hosts, their
occurrence and mutual inclinations, from DR4 epoch astrometry.

## Source

- **Paper list:** https://www.cosmos.esa.int/web/gaia/dr4-papers is an HTML table (Category / Title / ADS / arXiv). It has **86 entries**: 2 Summary, 5 Gaia Collaboration, 15 Performance Verification and 64 Processing.
  - The list was announced on the Gaia news page on 2026-06-29 ("Over 60 papers are expected"). Themes had already been posted on 2026-03-06.
  - No lead authors are listed. Titles were parsed from the raw HTML.
  - The planet-related titles were re-checked against the live page today.

## Overlap ranking

### HIGH: could pre-empt or directly compete with us

| # | Category | Title (verbatim) | Why |
|---|---|---|---|
| 20 | PV | Gaia Data Release 4: The Exoplanet Candidate Catalogue | DPAC's headline astrometric planet list. Any strong orbit around one of our hosts will already be in it, so we must cross-match and cite it |
| 34 | Proc | Gaia Data Release 4: Astrometric orbit determination for exoplanet discovery | The DPAC planet-orbit pipeline. It defines the thresholds we will be compared against, and likely validates on known hosts |
| 75 | Proc | Gaia Data Release 4: Exoplanet orbit processing with radial velocity constraints | Joint astrometry and RV fits for **known RV planet systems**. Hot-Jupiter hosts with RV-detected outer companions are natural targets, so this is the closest to a direct race on individual systems |
| 15 | PV | Gaia Data Release 4: Exploring the Brown Dwarf Desert and the Architectures of Multiple Systems from Stellar Companions down to Giant Planets | "Architectures of multiple systems … down to giant planets". Could include known-planet hosts and their geometry |

### MEDIUM: will include our targets incidentally; cite and use for vetting

| # | Title | Use |
|---|---|---|
| 51 | Gaia Data Release 4: New transiting planet candidates and Gaia's view of known exoplanets | Covers our exact host list. Probably photometric (it is listed with the variability papers), but "Gaia's view of known exoplanets" may include astrometry or RUWE. Read on day 1 |
| 33 | Gaia Data Release 4: Astrometric binaries | NSS `nss_two_body_orbit`. Hosts with stellar or brown-dwarf orbits = impostor check |
| 59 | Gaia Data Release 4: Spectroscopic multiple star processing | SB1/SB2 from epoch RVs = impostor check |
| 36 | Gaia Data Release 4: Spectroscopic characterisation of astrometric binaries by line integral analysis of Gaia RVS spectra | Luminous-companion test for astrometric orbits |
| 13 | Gaia Data Release 4: High accuracy dynamics of the nearest stars | Few of our hosts are this close |
| 14 | Gaia Data Release 4: The mass-luminosity relation | Flux-ratio and photocentre calibration for impostor modelling |

### LOW: methods and data papers we must read and cite

- Mission and release overview:
  - #1 The Gaia Mission
  - #2 Content and properties
- Epoch astrometry and its systematics:
  - #29 The astrometric solution (epoch-astrometry noise model)
  - #24 Modelling of drift-scan-related effects in Gaia's PSF (arXiv:2602.20906)
  - #32 Catalogue validation
  - #28 Parallax-bias recipe
  - #30 The Gaia ESA archive
- Epoch RVs:
  - #39 Epoch RVS spectra and radial velocities
  - #40 Derivation and publication of the RVS data products
  - #46 Properties and validation of the radial velocities
  - #42 Statistical parameters for RV time series
  - #35 RV outliers using machine learning
  - #23 Deblending the RVS spectra
  - #58 SB2 RVs and orbits
- Neighbours and photocentre contamination:
  - #76 Identification and resolution of close source pairs
  - #78 Astrometric processing of resolved pairs
  - #68 Source environment analysis
- Stellar masses:
  - #84 FLAME stellar masses and ages (arXiv:2607.00264), possible host masses
- Worked example:
  - #7 Discovery of a dormant 33 M☉ black hole in pre-release Gaia astrometry, an orbit fit to epoch astrometry. Its code is `esa/gaia-bhthree`.

All other titles (Solar System, variability, galaxies, QSOs, photometry calibration,
Local Group, etc.) do not overlap.

### Assessment of DPAC overlap

No DPAC title targets **outer companions of transiting hot/warm Jupiters** or
**mutual inclinations**. DPAC will, however, publish:
- a planet candidate catalogue (#20);
- known-system orbits with RV constraints (#75).

Together these could already contain orbits for the brightest and nearest of our
hosts. **Our analysis has to be framed as built on #20 and #75**, not in competition
with them. It would be a uniform search plus upper limits over *all* transiting giant
hosts, with an inclination comparison. Whenever DPAC published an orbit, we cite it
and compare.

## Non-DPAC teams (higher risk than DPAC)

Verified on arXiv abstract pages 2026-10-11 unless marked.

1. **Xiao, Teng, Huang, Feng, Lai, Dai & Liu**, "Transiting Planetary Systems with Distant Giant Companions Remain Moderately Coplanar", arXiv:2609.01377 (1 Sep 2026), accepted to ApJL.
   - Method: RV plus Hipparcos–Gaia astrometry for 19 transiting systems with outer giants, analysed with a hierarchical Bayesian model.
   - Result: minimum mutual inclination is Rayleigh with σ ≈ 15.8° (+2.8/−2.6°), preferred over isotropic. Giant-inner-planet systems hint at being more coplanar.
   - They state "Future Gaia DR4 astrometry will enable more robust population-level studies".
   - **This is our science question at population level, from the Feng group, which has the RV plus astrometry machinery ready.**
2. **Do Ó, Thompson, Blakely, … Mawet, Blunt, Bowler et al.**, "Astrometric Evidence for Outer Giant Planets in Known Exoplanet Systems", arXiv:2609.37160 (29 Sep 2026).
   - Method: G23H (Hipparcos plus DR2 plus DR3).
   - Result: 29 outer-companion candidates in known transiting and RV systems.
   - They are explicitly flagged as priority targets "for joint orbit analysis once Gaia DR4 data are available".
   - The Octofitter team will be fitting DR4 epoch data for these systems.
3. **Thompson, Blakely, Xuan et al.**, G23H catalogue paper, arXiv:2602.00235 (Jan 2026): independent astrometric support for 94 of 120 Jovian systems from the NASA Exoplanet Archive.
4. **Ranshaw, Wang, Burt, … S. Wang**, "The Discovery of K2-232c…", arXiv:2609.09077 (Sep 2026).
   - Outer companions are more eccentric in hot-Jupiter systems than in warm-Jupiter systems.
   - A DR4 detection forecast was mentioned in a search snippet. *Not verified in the abstract; check the PDF.*
5. **Shariat, El-Badry, S. Wang et al.**, "Stellar companions sculpt hot Jupiter formation…", arXiv:2609.19249 (Sep 2026): resolved stellar companions to 147 hot Jupiters. **Useful for our known-stellar-companion and impostor list.**
6. **Yang, Lu, Yahalomi & Winn**, arXiv:2606.02705 (Jun 2026): a roughly 34 M_J companion to KELT-20, from acceleration plus TTVs.
7. **Lammers & Winn**, arXiv:2511.04673: DR4 yield forecast. Mock catalogues are at `CalebLammers/GaiaForecasts`.
8. **GEODES I**, Van Zandt, Bowler, Petigura et al., arXiv:2605.23018: a vetting and validation campaign for Gaia giant-planet candidates. Not focused on transiting hosts.
9. **Not verified** (search snippets only):
   - a NExScI 2026 poster (Mazziotti) arguing the DPAC exoplanet catalogue may be conservative because of a RUWE < 1.4 cut;
   - Sozzetti ExoPAG slides on DR4 exoplanet deliverables;
   - an Astrobites post (5 Oct 2026) on vetting DR4 planets.

## Bottom line

- **DPAC overlap: moderate.** Nothing in DPAC answers our question directly. The candidate catalogue (#20) and the RV-constrained orbits (#75) will cover the strongest individual detections, so we cite them and compare.
- **Competition: high.** Two well-equipped external groups (Feng/Xiao; Octofitter/G23H) have publicly said they will apply DR4 to outer companions of known, including transiting, planet hosts. One of them already did the mutual-inclination population study with pre-DR4 data.
- **Room for us:**
  - a fast, uniform, release-day pass over all ~790 transiting giant hosts;
  - the **detection and upper-limit table**, i.e. completeness-corrected non-detections, which the targeted papers will not provide;
  - DR4 epoch-RV impostor vetting.

  This is publishable as an RNAAS if it is out within ~2 weeks and cites the groups above.
