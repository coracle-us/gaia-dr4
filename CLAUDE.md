# Gaia DR4 Companions — project briefing

Hobby research project (Jon McBee + Claude). Goal: be ready on **2 December 2026**, the
day Gaia Data Release 4 goes public, to search Gaia's epoch astrometry for **outer giant
companions to stars that host transiting giant planets**, and publish fast as a
**Research Note of the AAS (RNAAS)**.

Read this file first. It is the source of truth for scope, schedule and rules.

---

## 1. Why now

- **Gaia DR4 releases 2 Dec 2026** (ESA Gaia news, May and Oct 2026 notices). The Gaia
  ESA Archive goes offline 25 Nov 2026 12:00 CET and reopens 2 Dec with DR4.
- DR4 covers the first **5.5 years** of observations and, for the first time, publishes
  **per-observation (epoch) astrometry**: the individual position measurements a planet's
  gravitational tug would wiggle. DR3 only had averaged positions plus a small set of
  pre-fit orbital solutions.
- Scale and span (ESA Image of the Week, 3 Jul 2026): time series for 2.8 billion sources,
  142 billion detections, 66 months of data from **mid-2014 to early 2020**. Published
  time series include positions, G/BP/RP photometry, BP/RP and RVS spectra, and
  **epoch radial velocities for sources with G_RVS ≤ 14**. Epoch RVs are a big deal for
  us: they help unmask impostor binaries and can add a constraint on companion orbits.
- DR4 also introduces "residual images" (preview, Dec 2025), which may reveal dusty
  circumstellar material (debris/protoplanetary disks). Not our focus, but note it.
- Previews index: https://www.cosmos.esa.int/web/gaia/dr4-previews ; also see the DR4
  content and DR4 software tools pages linked from there.
- Forecasts (Lammers et al., arXiv:2511.04673) predict ~7,500 ± 2,100 planet detections
  in DR4, mostly super-Jupiters (3–13 M_J) on 2–5 AU orbits around GKM stars within
  ~500 pc, and warn that unresolved binaries create impostors.
- ESA has already published: a **draft DR4 data model**, a **pre-release sample of epoch
  astrometry** for 12 stars (June 2026), the **released** Python package for fitting
  single-star models to epoch astrometry (`gaiasupdate`, PyPI, June 2026; no orbit model),
  and a list of the **DR4 collaboration paper titles** (86 entries). All on the Gaia news / DR4 pages at
  https://www.cosmos.esa.int/web/gaia/ . Mock planet + impostor-binary catalogs from the
  forecast paper: GitHub `CalebLammers/GaiaForecasts`.
- The first weeks after a big release are when an independent team can contribute
  something new. We use the next seven weeks to be ready.

## 2. The science question (decided 2026-10-11: A, narrowed)

**Do stars hosting transiting hot/warm Jupiters have outer giant companions, and how
often?**

Giant planets form far out; something moves hot and warm Jupiters inward. One leading
idea (high-eccentricity migration) needs an outer perturber to pump the inner planet's
eccentricity. RV surveys ("Friends of Hot Jupiters", Knutson et al. 2014; Bryan et al.
2016) and Hipparcos–Gaia acceleration studies have measured companion rates, but
astrometry gives the *orbit and true mass* (inclination included), and for transiting
hosts the inner planet's inclination is known (~90°), so **mutual inclinations** become
measurable. Aligned vs. misaligned outer companions discriminate between migration
pathways. That's a sharp, timely question that DR4 enables.

Ties to the sister project (`duotransits` repo): warm Jupiters are the population we're
characterizing there.

**Decision (2026-10-11, Jon): option A, narrowed.** Framing for the RNAAS note:
*DR4 constraints on outer giant companions to all ~790 transiting giant-planet hosts.*
- Headline product: one table, per host, of detections, acceleration-only signals, and
  completeness-corrected upper limits (a uniform population result; targeted competitors
  will publish selected detections, not limits).
- Impostor vetting from option C folded in: DR4 epoch RVs (703 hosts), DPAC NSS tables,
  colour–magnitude position, RUWE > 1.4 hosts flagged (52).
- Mutual inclinations: minimum ΔI reported only for detections, framed as building on
  Xiao et al. 2026 (arXiv:2609.01377); cite and compare DPAC's Exoplanet Candidate
  Catalogue and RV-constrained orbits.
- Fallback if DPAC already publishes orbits for most sensitive hosts: limits + vetting.
Rationale (sensitivity, overlap, competition, feasibility): `notes/dr4_paper_overlap.md`,
`notes/literature.md`.

**Alternatives considered in Phase 0** (rejected; kept for the record):
- (A) the question above: transiting giant-planet hosts (default);
- (B) companions to nearby M dwarfs or white dwarfs (high sensitivity, high impostor risk);
- (C) a validation/impostor study: which DR4 candidates are binaries in disguise.
Choose based on (1) sensitivity: enough targets where a Jupiter at 2–5 AU is detectable;
(2) **overlap with the 60+ DR4 collaboration papers**: do not race DPAC on their own
headline analysis; (3) what we can finish and vet within ~2 weeks of release.

## 3. Sensitivity rule of thumb

Astrometric signature: **α [arcsec] = (M_p / M_★) · (a [AU]) / (d [pc])**.
Jupiter (M_p/M_☉ ≈ 9.5e-4) at 3 AU around 1 M_☉ at 100 pc → α ≈ 29 µas.
Gaia per-epoch precision is roughly tens to ~100+ µas for bright stars (G ≈ 8–13),
degrading for fainter ones; with dozens of epochs, α of a few tens of µas is reachable
near the bright end. Verify these numbers against DR4 documentation; don't trust this
file. Use them to rank targets by expected S/N.

## 4. Schedule

### Phase 0 — Scope (now → ~1 Nov)
1. Read the DR4 draft data model: table names, epoch-astrometry fields (times, scan angle,
   along-scan positions, uncertainties, flags), access method (TAP/ADQL,
   DataLink, bulk). Write `notes/dr4_data_model.md`.
2. Download the pre-release epoch-astrometry sample (announced on Gaia news, June 2026;
   it may not be listed on the previews page); load it; plot it. Note how many targets
   have G_RVS ≤ 14 (epoch RVs available).
3. Read the DR4 paper title list; flag every paper overlapping our question
   (`notes/dr4_paper_overlap.md`).
4. Literature: DR3 NSS orbital solutions and Hipparcos–Gaia acceleration (HGCA, Brandt)
   results for transiting-planet hosts; RV companion surveys; Lammers et al. forecasts.
5. Target list: transiting giant planets (M_p or R_p ⇒ giant) from the NASA Exoplanet
   Archive, with Gaia DR3 source_id, distance, G mag, RUWE. Compute expected sensitivity
   to a 1 M_J and 5 M_J companion at 1–5 AU. Save `data/targets.csv`.
6. Decide the question (A/B/C) and write it down here. **Done: A, narrowed (§2).**

### Phase 1 — Build and validate (Nov)
1. Model: 6-parameter single-star astrometry (5 astrometric + pseudo-colour, as in the
   DR4 core solution) + Keplerian orbit (photocenter), using the
   along-scan 1-D formulation Gaia actually measures. Linear Thiele–Innes parameters +
   nonlinear (P, e, T_p) search; periodogram over P ≤ ~4 yr (≤ the DR4 baseline).
2. Use ESA's single-star fitting package (`gaiasupdate`, already released); our orbital
   extension must reproduce its single-star fits exactly when the orbit term is off.
   Note the DR4 core solution is 6-parameter (`6p_constrained_colour`: 5 astrometric +
   pseudo-colour), and the epoch data are measured AL positions (`centroid_pos_al`,
   relative to a per-source reference point `ra0`/`dec0`), **not residuals**: we fit the
   full model to the positions. Weight by `centroid_pos_error_al` only, not
   `agis_source_excess_noise` (which absorbs orbital signal). See `notes/dr4_data_model.md`.
3. **Injection–recovery**: simulate epoch astrometry with realistic scan times/angles and
   noise (use the pre-release sample's sampling, a scanning-law tool, or DR4 docs), inject
   planets and impostor binaries (GaiaForecasts mocks), measure completeness and false
   positives.
4. Cross-check on DR3 NSS orbital solutions with independent confirmation.
5. **Release-day dry run**: scripts that, given only DR4 access, pull epoch data for all
   targets, fit, and write a results table with no manual steps.

### Phase 2 — Release (2 Dec →)
1. Run the pipeline on the target list on release day (expect archive load; be patient,
   use async TAP jobs, cache everything).
2. Read the DR4 collaboration papers relevant to us immediately, compare our detections
   with any published DPAC solutions, cite them.
3. Vet every candidate (§5). Then write the RNAAS note: one table (detections + limits),
   ≤1,500 words.

## 5. Rules of the road

- **Verification over speed.** Speed matters, but a wrong claim is worse than a late one.
  Every claimed companion passes: significance vs. single-star model (ΔBIC or FAP from
  injection tests); RUWE/excess-noise sanity; impostor checks (luminous binary
  photocenter, known stellar companions, DR4 epoch RVs where G_RVS ≤ 14 plus any
  literature RVs); stability of fit
  across reasonable data cuts; independent re-derivation (fresh script or separate agent).
- Upper limits / non-detections are publishable and scientifically useful here.
- Never edit raw downloads (`data/raw/`); derived products go elsewhere.
- Every number in the note traces to a file and a script in this repo.
- Cite, don't claim, anything DPAC or others already published.
- Don't download huge bulk files when a targeted query works; we need hundreds of stars,
  not a billion.

## 6. Publication: Research Notes of the AAS

- Moderated, not peer reviewed; fast; ADS-indexed and citable; no author charges (Oct 2026).
- ≤1,500 words, abstract required, **one figure or one table, not both**. AASTeX v7+.
- Author: Jon McBee (independent; consider Coracle as affiliation). Claude is not an
  author; acknowledge AI assistance and **check AAS journals' current AI policy first.**

## 7. Environment

- Windows machine (Jon). Python venv in `.venv`. Network access to ESA/NASA archives is
  direct from this machine.
- Useful: `astroquery` (Gaia TAP), `pyvo`, `numpy`/`scipy`, `emcee` or `dynesty` for
  posteriors, `matplotlib`. Official DR4 epoch-astrometry package: `gaiasupdate`
  (installed, in `requirements.txt`). Third-party DR4-ready orbit fitter for
  cross-checks: `exogaia` (T. Stolker).

## 8. Status (updated 2026-10-11)

- DR4 date confirmed (2 Dec 2026). Question decided: A, narrowed (§2).
- Phase 0 steps 1–6 done: `notes/dr4_data_model.md`, `notes/dr4_paper_overlap.md`,
  `notes/literature.md`, `data/targets.csv` (`scripts/build_targets.py`), pre-release plot
  (`scripts/plot_prerelease.py`), prior-astrometry cross-match and known-companion
  validation set (`scripts/crossmatch_prior_astrometry.py`, `scripts/validation_systems.py`).
  Step 6 done (A, narrowed). Phase 0 complete; next is Phase 1.
- Sensitivity (edge-on, Lammers & Winn Δχ² > 50, P ≤ 4 yr): 73 hosts reach ≤ 5 M_J,
  370 reach ≤ 13 M_J. Most known outer companions have P > 4 yr, so Phase 1 needs
  acceleration models and a joint astrometry + RV mode (see `notes/literature.md` §0).
