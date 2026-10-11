# Phase 1 validation log

## 2026-10-11 — orbit fitter, noise model, Gaia-4b

Code is in the `companions/` package. Tests are in `tests/` (`pytest`, 45 tests). Numbers below come from the scripts named alongside them.

### 1. Single-star model reproduces ESA `gaiasupdate`

- The DR4-like `6p_constrained_colour` fit agrees to 1e-9 relative tolerance on all 12 pre-release sources (`tests/test_linear.py`). The 5p fit also matches.
- Undocumented detail: `gaiasupdate` multiplies the archive's `colour_factor_al` by −1e3, converting −dw/dν_eff in mas·µm to dw/dν_eff in mas·nm. Its colour prior of 0.085e-3 nm⁻¹ (= 0.085 µm⁻¹) refers to the converted column.

### 2. Noise in the pre-release epoch data

**Single stars** (G = 9–20):
- CCD errors are well calibrated: per-CCD χ²/dof is 0.95–1.08.
- Per-transit means are 0.8–1.4, so at most a small per-transit excess (G ≈ 14).

**Bright or orbit stars:** there is excess scatter within transits, which no orbit can produce.

| Star | G | Within-transit χ²/dof |
|---|---|---|
| Gaia-4 | 11.9 | 1.38 |
| BH3 | 11.2 | 3.3 |
| HD 114762 | 7.2 | 9.9 |

**What we did about it**
- `companions.noise` models per-CCD jitter `s_ccd` plus a per-transit correlated jitter `s_tr`. The transit term uses a block-diagonal covariance (Sherman–Morrison), checked against dense algebra.
- Both jitters are fitted by maximum likelihood for each model. The detection statistic is dlnl2 = 2(lnL_orbit − lnL_star).
- Test with pure per-transit noise (s_tr = 0.06 mas, σ = 0.08 mas, 4 seeds):

  | Statistic | Range |
  |---|---|
  | Fixed-noise per-CCD Δχ² | 55–144 (false "detections" at 50) |
  | ML dlnl2 | 9–24 |

  The ML values sit at the pure-noise level of about 20 that Lammers & Winn report.
- For Gaia-4 the excess is per-CCD (s_ccd = 0.047 mas, s_tr ≈ 0.005 mas). `exogaia` gets the same jitter, 0.049 mas.
- **Do not use `agis_source_excess_noise` for orbit fits.** For Gaia-4 it gives χ²/dof = 0.46 with the orbit included, because it has already absorbed the planet.

### 3. Angle conventions

Our Thiele–Innes → Campbell conversion, applied to the DR3 `nss_two_body_orbit` solution for Gaia-4, gives:

| | a₀ (mas) | i (°) | ω (°) | Ω (°) |
|---|---|---|---|---|
| Our conversion | 0.313 | 115.4 | 185.9 | 159.4 |
| nsstools (Stefansson+25, Table 2) | 0.312 | 115.4 | 185 | 159.4 |

These match, and the check is pinned in `tests/test_orbit.py`. Proper-motion signs and sizes match DR3 (pmRA −75.55 vs −75.41 mas/yr).

### 4. Gaia-4b from DR4 pre-release epochs

Sources: `scripts/validate_gaia4.py` → `results/gaia4_validation.csv`; `exogaia` via `scripts/crosscheck_exogaia.py`.

- **Detection:** dlnl2 = 252.
- **Posterior:** emcee, 64k samples, acceptance 0.54.
- **Published values:** Stefansson+25 (arXiv:2410.05654), the joint Gaia DR3 + RV fit.

| | Ours (DR4 pre-release, 6p + s_ccd + s_tr) | `exogaia` 0.1.3 (same data, 5p + jitter) | Published Gaia + RV | DR3 Gaia only |
|---|---|---|---|---|
| P (d) | 578.7 ± 3.6 | 573.1 ± 3.2 | 571.3 ± 1.4 | 564 ± 11 |
| e | 0.43 ± 0.08 | 0.29 ± 0.06 | 0.338 ± 0.025 | 0.51 ± 0.20 |
| i (°) | 120.5 ± 2.8 | 121.8 ± 2.7 | 116.9 ± 4.3 | 115.4 ± 4.6 |
| Ω (°) | 174.3 ± 2.8 | 173.7 ± 2.6 | 158.6 ± 5.1 | 159.4 ± 5.0 |
| ω (°) | 215 ± 10 | 181 ± 13 | 180.3 ± 5 | 185 ± 16 |
| a₀ (mas) | 0.258 ± 0.012 | 0.247 ± 0.008 | — | 0.312 ± 0.040 |
| M_c (M_J), M★ = 0.644 | 11.0 ± 0.5 | 10.6 | 11.8 ± 0.7 | — |

**Reading the table**
- **Inclination, mass and a₀:** consistent with the published values at about 1σ.
- **Node Ω:** DR4 gives 174° and DR3 gave 159°, a 2.6σ difference. `exogaia` independently gets the same 174° from the DR4 data, and our conversion reproduces the DR3 numbers exactly. So the difference is in the data (DR4 versus DR3 solution), not in our code.
- **P and e versus `exogaia`:** with `exogaia`'s exact setup (5p, single jitter), our likelihood still peaks at P = 579 d.
  - The profile likelihood (`results/gaia4_period_profile.csv`) is smooth. P and e are strongly correlated along a ridge, with e = 0.34 at P = 573 d.
  - `exogaia`'s solution lies at −2ΔlnL = 2.9 below the maximum, after only 14 function evaluations. It stopped short of the optimum.
  - The published P = 571.3 d (fixed by the RVs) lies on the same ridge at −2ΔlnL = 4.7, i.e. 2.2σ from the DR4-only peak. Expected value along the ridge: e ≈ 0.32.

**Lessons for Phase 1**
- Report P and e jointly, because they are correlated. Use RV constraints where they exist; this is the joint-RV mode.
- Don't trust a single local optimiser. Our multi-start grid plus refinement found the true maximum where `exogaia`'s did not.
- Ω and ω from DR4 can differ from DR3 at the 2–3σ level for real planets. Treat DR3 orbit angles as approximate when comparing.

### Open

- Joint astrometry + RV fit. Confirm the ω convention (photocentre vs. primary) before combining with RVs.
- Acceleration (7p/9p) models.
- Null and injection calibration of dlnl2 on realistic DR4 sampling (scanning law), including edge-on orbits.

## 2026-10-11 (later): acceleration models and joint astrometry + RV fit

### Acceleration models (`companions/accel.py`)

- 7- and 9-parameter models with the same definition as DR3 `nss_acceleration_astro`: g t²/2 and ġ t³/6 in α* and δ, at J2017.5. Both use ML jitters.
- On pure noise, dlnl2 behaves like χ² with 4 degrees of freedom (max of 5 seeds < 18.5).
- A 14-year orbit (a₀ = 0.6 mas) is detected as acceleration with dlnl2 > 50.

### Joint fit (`companions/joint.py`)

**Parametrisation**
- Nonlinear: P, e, M0, a₀, cos i, ω, Ω, s_ccd, s_tr, and one RV jitter per instrument.
- Profiled linear: single-star astrometry, RV offsets, and the K of inner transiting planets (circular, with fixed P and T_c).
- K is derived as 2π (a₀/ϖ) sin i / (P √(1−e²)).

**Convention**
- Gaia DR3 AstroSpectroSB1 / nsstools: the same ω in the astrometry and in RV = K[cos(ν+ω) + e cos ω].
- Tested as RV = dz/dt, with z = C X + H Y.
- Astrometry alone is degenerate under (ω, Ω) ↔ (ω+π, Ω+π). The fit tries both branches and the RVs choose. On synthetic data the correct branch, K (within 1%) and the offsets are recovered.

### Gaia-4b: DR4 pre-release astrometry + Stefansson+25 RVs

Sources: `scripts/validate_gaia4_joint.py` → `results/gaia4_joint_validation.csv`. The 23 RVs come from AJ 169, 107, Table 5 (machine-readable), saved as `data/raw/literature_rv/`.

**RV-only fit with our RV model**

| | P (d) | e | K (m/s) | ω (°) |
|---|---|---|---|---|
| Ours | 506 | 0.247 | 353.6 | 181.7 |
| Published RV-only | 513 (+21/−16) | 0.263 | 355 | 181.8 |

The RV model and its convention are verified.

**Joint fit**
- The RVs select the Ω = 172.6° branch, matching the published 158.6° branch rather than its +180° twin. This is pinned by a slow test.
- Comparison with the published joint fit, which used the DR3 orbit with the same RVs:

  | Parameter | Our joint | Published joint | Difference |
  |---|---|---|---|
  | P (d) | 575.5 | 571.3 ± 1.4 | 3.1σ |
  | e | 0.32 | 0.338 | 0.6σ |
  | i (°) | 120.3 | 116.9 | 0.8σ |
  | Ω (°) | 172.6 | 158.6 | 2.7σ |
  | ω (°) | 200 | 180.3 | 3.8σ |
  | K (m/s) | 320 | 368 | 2.4σ |
  | M (M_J) | 10.8 | 11.8 | 1.4σ |

**Tension between DR4 astrometry and the RVs**
- 2[lnL_astro-only + lnL_RV-only − lnL_joint] = 17.8, with about 5 shared constraints (P, e, T_p, ω, and K ↔ a₀ sin i). That is roughly 3σ.
- Main drivers:
  - DR4 gives a₀ = 0.25 mas, but the RV K needs a₀ ≈ 0.29 mas at i = 120°.
  - DR4 gives ω = 215 ± 10°, against 182 ± 7° from the RVs.
- **Our code is not the cause.** The RV-only fit matches the paper, the astrometry-only fit matches exogaia, and the conventions are pinned by tests.
- **Possible causes:**
  - pre-release (RC3) processing systematics;
  - correlated noise not captured by our two jitter terms;
  - a second companion or a trend.

**Release-day action.** Re-fit Gaia-4 with the final DR4 data. If astrometry–RV tension of this size recurs on known systems, it will limit our joint-fit claims (masses and mutual inclinations).

## 2026-10-11 (later): DR4 sampling, null calibration, injection–recovery

### Sampling (`companions/scanlaw.py`, `companions/simulate.py`)

**Transit prediction**
- Source: ESA GOST predictions on a HEALPix nside-64 grid. Only the data tables from the gaiascanlaw 0.2.0 wheel are used; its code needs healpy, which doesn't build on Windows. Raw files are in `data/raw/gaiascanlaw/`.
- **Time convention:** the tables' times are Julian years. astropy's calendar `decimalyear` is off by hours.
- **Validation against the 9 pre-release stars:**

  | Check | Result |
  |---|---|
  | Transits matched one-to-one | 37–91% |
  | Timing | ~10–40 s |
  | Scan angle | < 1° |
  | Predicted / real transit count | 0.9–1.4 |

- The grid pixel (0.9°) is comparable to the across-scan FoV width (0.7°), so the predictions are statistically representative but not exact. Release day: use the DR4 `commanded_scan_law` table.

**Parallax factors**
- Computed for Gaia at L2. They match the real `parallax_factor_al` to better than 0.02.

**Simulator settings**
- 9 AF CCDs per transit; 10% of transits rejected.
- Holl+2023 per-CCD σ (`companions/precision.py`, shared with `scripts/build_targets.py`).
- Optional s_ccd and s_tr jitter, and an injected orbit.

### Null calibration

Sources: `scripts/calibrate_null.py` → `results/calibration_null.csv` (10,000 realisations), `scripts/analyze_calibration.py` → `results/thresholds.csv`.

**Setup:** single stars at real target positions and G magnitudes, with excess noise s_ccd, s_tr ~ U(0, 0.05) mas, run through the full pipeline statistics.

| Statistic | Median | 99.9% | Threshold at FAP 1e-3 | Threshold at FAP 1e-4 | Max of 10⁴ |
|---|---|---|---|---|---|
| Orbit dlnl2 | 14.8 | 35.1 | 35.1 | 42.6 (exponential tail) | 50.7 |
| 7p acceleration | 1.5 | 16.4 | 16.4 | 21.7 | — |
| 9p acceleration | 3.6 | 20.8 | 20.8 | 27.1 | — |

- The orbit threshold is stable across magnitude: 39–45 for G < 11, 11–13 and > 13.
- Lammers & Winn's Δχ² > 50 is slightly more conservative than our FAP-1e-4 threshold.
- **False-positive mode:** the extreme nulls are e > 0.9 fits to a few noisy transits, often where the single-star ML fit set s_tr to 0. Vetting rule: any detection with e > 0.8 needs a leave-transits-out (jackknife) check.

### Injection–recovery

Source: `scripts/injection_recovery.py` → `results/injection_recovery.csv` (6,000 injections).

**Setup**
- P log-uniform over 0.2–8 yr, e ~ U(0, 0.6).
- Half isotropic, half aligned (i = 90 ± 3°).
- Amplitude set by the expected full-orbit Δχ², log-uniform over 3–3000.

**Detection at FAP 1e-4 per star** (`results/completeness.csv`, `figures/completeness.png`)
- **P ≤ 4 yr, aligned:** logistic fits give 50% completeness at an expected Δχ² of 74–84, and 90% at 126–166. Isotropic orbits give the same curve once the geometry factor is included. The period is right (within 10%) for over 90% of detections at P < 3 yr.
- **P = 4–5.5 yr:** orbit detections still work, but the period is right in only 55–90% of them. The acceleration models recover these well.
- **P > 5.5 yr:** strong signals trigger the orbit statistic with a wrong period (about 0% correct). Use the acceleration models, which reach 60–96% at high signal.

### Calibrated per-host sensitivity

Source: `scripts/target_sensitivity.py` → `results/target_sensitivity.csv`.

Aligned (edge-on) companions, detected at FAP 1e-4 per star. Entries are the number of hosts whose detectable mass at that completeness is at most the stated limit:

| P | 50% complete ≤ 1 M_J | ≤ 2 M_J | ≤ 5 M_J | ≤ 13 M_J | 90% complete ≤ 5 M_J | 90% complete ≤ 13 M_J |
|---|---|---|---|---|---|---|
| 1 yr | 0 | 1 | 6 | 68 | 1 | 23 |
| 2 yr | 0 | 1 | 17 | 189 | 12 | 119 |
| 3 yr | 1 | 2 | 39 | 266 | 17 | 178 |
| 4 yr | 1 | 5 | 56 | 326 | 19 | 200 |

- Calibrated 50% masses are a median 1.13× the analytic estimate in `targets.csv` (`mmin_mj_edge_P4yr`).
- **Best host:** HD 189733, with m50 = 0.82 M_J at P = 3 yr.
- **The 39 hosts with m50 ≤ 5 M_J at P = 3 yr:**
  - all have epoch RVs;
  - 32 are hot-Jupiter hosts;
  - 3 have RUWE > 1.4.

### Caveats

- **Noise prior:** the excess-noise prior (U(0, 0.05) mas) is a guess informed by 12 pre-release stars. On release day, recalibrate the null on real DR4 residuals of single stars, e.g. RUWE < 1.1 stars near our targets.
- **Sampling:** this is GOST-grid sampling, not exact. Repeat the calibration with the DR4 `commanded_scan_law` and `gaia_gap` tables.
- **Simulation model:** 5p only (no colour term) and dark companions only. Luminous binaries (impostors) are not simulated yet; the GaiaForecasts impostor mocks are the next step.
