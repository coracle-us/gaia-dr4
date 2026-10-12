# Gaia DR4 release-day runbook

**Timeline**
- **25 Nov 2026, 12:00 CET:** the Gaia archive goes offline.
- **2 Dec 2026, around 12:00 CET:** DR4 is released.

**Goal:** a vetted results table for all ~780 transiting giant-planet hosts within days, then the RNAAS note.

Everything below runs from the repository root in the project venv (`.venv\Scripts\python`). The archive will be under heavy load, so every download step caches and resumes. **Re-running a command never re-downloads what is already in `data/raw/dr4/`.**

---

## Before release (by 24 Nov)

1. **Fast tests:**
   ```bash
   .venv/Scripts/python -m pytest
   ```
2. **Real-data tests:** pre-release Gaia-4, the joint RV fit, and live DataLink.
   ```bash
   .venv/Scripts/python -m pytest -m slow
   ```
3. **Full dry run:** all hosts on simulated DR4 data, with injected planets.
   ```bash
   .venv/Scripts/python scripts/release_day.py --source simulate --inject 25 --out results/dryrun
   ```
   Check:
   - the injected planets are recovered (`results/dryrun/injections.csv` against `results.csv`);
   - false positives ≈ 0;
   - `n_errors` = 0 in `run_summary.json`;
   - the dashboard `results/dryrun/dashboard/index.html` shows the candidates, their plots and the injected planets.
4. **Control stars** for the real-data null calibration. This is already done if `data/control_stars.csv` exists.
   ```bash
   .venv/Scripts/python scripts/calibrate_null_real.py --prepare-dr3 200
   ```
5. **Git:** commit, then tag the code state used on the day.
   ```bash
   git tag pre-dr4
   ```

## Release day: first hour (verify before bulk downloads)

Check the following on the DR4 archive pages and the DR4 documentation, then edit `config/release_day.json`:

| Setting in `config/release_day.json` | Check |
|---|---|
| `data_release` | The exact release string DataLink expects (e.g. `"Gaia DR4"`). Test with one source: `Gaia.load_data(ids=[<id>], data_release=..., retrieval_type="EPOCH_ASTROMETRY", data_structure="RAW")`. |
| `retrieval_type`, `rv_retrieval_type` | `EPOCH_ASTROMETRY`, `EPOCH_PARAMETERS_RVS_SINGLE` (draft data model). |
| `dpac_tables`, `control.gaia_source_table`, `control.nss_tables` | DR4 TAP schema and table names (`gaiadr4.*`?). |
| `rv_summary` | The TAP table and columns holding the DR4 SB2 indicator (draft data model: `all_source_rvs.has_epoch_rvs`, `fraction_double_lined`). This feeds the twin-impostor check and prediction H2 for TIC 147035057. |
| `dr4_source_id` | Do DR3 ids persist? If DR4 ships a DR3→DR4 crossmatch, implement `map_source_ids` in `scripts/release_day.py`. |
| Column names | Spot-check one epoch_astrometry file against `companions/epoch.py` (`obs_time_tcb`, `obs_time_bary_corr`, `scan_pos_angle`, `parallax_factor_al`, `colour_factor_al`, `centroid_pos_al`, `centroid_pos_error_al`, `used_by_agis_al`, `agis_source_excess_noise`). |

**Sanity check on one known system** (Gaia-4, source_id 1457486023639239296):
- Fetch it and run `run_cascade`.
- Expect an orbit with P ≈ 570–580 d and M ≈ 11 M_J.
- Compare against `results/gaia4_validation.csv`.

## Named targets

`data/extra_targets.csv` lists targets that are always run, even if they are not in the planet archive. Columns:
- `p_min_yr` / `p_max_yr`: an optional per-target period window;
- `rv_known_period_d`: a period known from photometry, used to fit Gaia epoch RVs at that period;
- `rv_test_periods_d` (`"P1;P2"`) and `rv_t0_bjd`: every candidate period, and the epoch of a primary eclipse or transit (BJD_TDB). Output `named_rv.csv` has one row per period with:
  - K ± σ and the 3σ upper limit;
  - Δχ² against a constant;
  - the offset of the fitted conjunction from T0 (it should be ≈ 0 if the eclipsing body is the RV companion);
- `note`.

Current entries:
- **TIC 147035057:** the expected 37.3-d photocentre orbit, plus the RV fit at the eclipse period. This settles single star vs SB2.
- **TIC 159614322:** excludes stellar-mass companions at 90 / 180 d (K upper limit at both periods).

The ephemerides (T0, periods) come from the duo-transits repo: `singles/results/vetting/147035057_ephemeris.json` and `singles/writeups/tic159614322/ephemeris.json`. That repo's DR4 checklist is `singles/TODO_gaia_dr4.md`.

Add rows for collaborators' targets: Gaia DR3 id, RA/Dec, G, distance and stellar mass are needed.

## Release day: main run

1. **Fetch and fit every host.**
   ```bash
   .venv/Scripts/python scripts/release_day.py --source datalink --out results/dr4
   ```
   What it does:
   - downloads epoch astrometry (raw zips go to `data/raw/dr4/epoch_astrometry/`);
   - downloads epoch RVs and DPAC NSS rows;
   - runs the cascade on 30 cores, with posteriors for orbit detections.

   Expected time: download depends on archive load; fitting takes about 6 min plus about 1 min per orbit detection.

   **Watch it:** open `results/dr4/dashboard/index.html` in a browser. It reloads every minute and shows the stage, hosts downloaded and fitted, errors, candidates as they appear (with vetting checklist and plots), and the named targets next to their pre-registered predictions. The pipeline rebuilds it as hosts finish; to force a rebuild from another terminal:
   ```bash
   .venv/Scripts/python scripts/dashboard.py results/dr4
   ```

   If downloads fail partway, rerun the same command (it resumes), or refit offline from what is cached:
   ```bash
   .venv/Scripts/python scripts/release_day.py --source cache --out results/dr4
   ```
2. **Real-data null calibration** (the noise prior is the main pre-release assumption).
   ```bash
   .venv/Scripts/python scripts/calibrate_null_real.py --k 5
   ```
   - Compare `results/thresholds_real.csv` with `results/thresholds.csv`.
   - If the real thresholds are higher, set `"thresholds_file": "results/thresholds_real.csv"` and rerun step 1 with `--source cache`.
3. **Population and tables.**
   ```bash
   .venv/Scripts/python scripts/population_analysis.py results/dr4/results.csv --out results/dr4
   ```
   ```bash
   .venv/Scripts/python scripts/make_table.py results/dr4/results.csv --out paper/
   ```
   ```bash
   .venv/Scripts/python scripts/dashboard.py results/dr4
   ```

## Vetting every candidate (CLAUDE.md §5; flags in `results.csv`)

| Flag / check | Action |
|---|---|
| `stability_confirmed` | Strong support (noise passes 3%). |
| `stability_not_confirmed` | Weak evidence against: inspect residuals. |
| `high_e_check` (e > 0.8) | Noise false positives look like this. Inspect transit-level residuals. |
| `marginal` (dlnl2 < 1.25 × threshold) | Report as a candidate only. |
| `impostor_not_excluded` (no literature RVs) | Check the DR4 SB2 table (`rvs_epoch_parameters_double`), epoch `vbroad` variability, and overluminosity (`results/impostor_vetting.csv` → `dM_G`). |
| `stellar_companion_check` | RUWE > 1.4 or a prior HGCA/PMa/DR3-NSS signal. Is the signal a known stellar companion? |
| `gaia_rv_variable` | Epoch RVs are not constant after removing the transiting planet: stellar companion? |
| `fraction_double_lined` > 0 | Double-lined RVS transits: a stellar binary. The dashboard marks the twin check as failed. A value of 0 does not exclude twins split by less than ~25 km/s. |
| `dpac_solutions` | DPAC already has a solution. **Compare and cite.** |
| `brown_dwarf_mass` / `stellar_mass` | Not a planet. Report as a companion of that class. |

**Recording decisions.** Checks done by hand (independent re-derivation, SB2 / overluminosity checks, known-companion searches) go in `results/dr4/manual_vetting.csv` with columns `hostname, check, status, note`. `check` is one of `significance, stability, impostor_twin, stellar_companion, gaia_epoch_rv, dpac, mass, eccentricity, rederivation`; `status` is `pass`, `warn`, `fail` or `pending`. Rebuild the dashboard after editing it. A candidate is ready to claim only when its dashboard status is `pass`.

**For each surviving detection:**
1. Run the joint fit with the literature RVs (`companions.joint.fit_joint_ml`).
2. Check consistency: the astrometry–RV tension statistic, as in `scripts/validate_gaia4_joint.py`.
3. Do an independent re-derivation (fresh script or a separate agent).

**Known pre-release lesson:** on Gaia-4, the DR4 pre-release astrometry and the RVs were in about 3σ tension. Report joint-fit masses and inclinations with that caveat if the tension recurs.

## After the run

- Read the DPAC exoplanet papers on day 1: *Exoplanet Candidate Catalogue*, *Astrometric orbit determination for exoplanet discovery*, and *Exoplanet orbit processing with RV constraints*. Cross-check every overlapping system.
- Recalibrate completeness if the real noise differs. Rerun `scripts/injection_recovery.py` with the noise prior set from the real control stars, then `scripts/target_sensitivity.py` and `scripts/analyze_calibration.py`.
- Every number in the note must trace to a file in `results/` and the script that wrote it.
