# Pre-registered Gaia DR4 predictions for named targets

Written **2026-10-12**, before Gaia DR4 (2 Dec 2026); TIC 238606535, TIC 184895645 and TIC 301328232 added the same day. The git commit and push timestamps of this file are the record. **Do not edit after release**; add outcomes in a separate section below.

**Sources**
- Numbers: `scripts/dr4_predictions.py` → `results/dr4_predictions.csv`.
- Stellar parameters, periods and DR3 values: McBee 2026, RNAAS (submitted), *Revisiting Stale TESS Single-transit Community Candidates*, and its per-target write-ups.
- Pipeline used on release day: `scripts/release_day.py`. These two stars are named targets in `data/extra_targets.csv`.

**Detection criteria** (calibrated before release)
- Astrometric orbit detection: dlnl2 > 42.6 (false-alarm probability 1e-4 per star; `results/thresholds.csv`).
- 50% completeness at an expected full-orbit Δχ² ≈ 80 (`results/completeness_models.csv`).
- Expected Δχ² = N_obs a0² (1 + cos² i) / (4 σ_eff²), with σ_eff² = σ_ccd(G)² + (0.025 mas)².

## TIC 147035057 (CTOI 147035057.01)

**System:** G = 8.83, 166 pc, M1 = 1.44 M☉, P = 37.31497 d, grazing eclipse (b ≈ 1.28).
- DR3 RV robust amplitude is 49.8 km/s, so K ≈ 25 km/s and the mass function is f(M) = 0.060 M☉.
- DR3 RUWE is 1.12.

| Hypothesis | M2 | Photocentre a0 | Expected Δχ² (astrometry) | Epoch RVs |
|---|---|---|---|---|
| H1: single-lined star | ≥ 0.64 M☉ (sin i = 1) | 0.49 mas | ≈ 4,000 | single-lined, K1 ≈ 25 km/s at 37.315 d |
| H2: near-equal SB2, q = 0.97 | ~1.2 M☉ each | 0.077 mas | ≈ 100 | double-lined, separation up to ≈ 85 km/s |
| H2: q = 0.99 | | 0.034 mas | ≈ 20 | double-lined |
| H2: q = 1.0 | | 0 | 0 | double-lined |

**If H1 (single-lined stellar companion), DR4 shows:**
1. An astrometric orbit with P = 37.31 ± 0.05 d, photocentre a0 = 0.4–0.6 mas, and inclination i = 87.6° or 92.4° (±3°), consistent with the grazing eclipse geometry. The detection is unmistakable (dlnl2 ≫ 42.6).
2. Single-lined epoch RVs following a 37.315-d orbit with K1 ≈ 25 km/s. For a circular orbit, the RV passes through the systemic velocity, decreasing, at the eclipse times (T0 = BTJD 1515.0661). The lack of a secondary eclipse suggests e > 0, which would shift this.

**If H2 (near-equal pair, q ≳ 0.99), DR4 shows:**
1. No significant astrometric orbit at 37.3 d (dlnl2 < 42.6).
2. The star in the double-lined products (`rvs_epoch_parameters_double`, fraction_double_lined > 0), with component separations up to ~85 km/s. Single-template epoch RVs scatter rather than following one clean K ≈ 25 km/s sinusoid.

**Our expectation: H2, or another configuration with little photocentre motion.**
- H1 predicts a DR3 RUWE of about 5 (a 0.49 mas orbit against σ_FoV ≈ 0.05 mas at G = 8.8), but DR3 measured 1.12.
- The modest DR3 astrometric_excess_noise (0.12 mas, 28σ) and the duplicated_source flag fit a small, non-zero motion.
- A strong 0.5 mas DR4 orbit would mean the DR3 RUWE was anomalous.
- Intermediate mass ratios (q ≈ 0.97) give a ~0.08 mas orbit (Δχ² ≈ 100), a marginal detection, together with SB2 flags.

## TIC 159614322 (CTOI 159614322.01)

**System:** G = 11.63, 401 pc, M1 = 1.11 M☉; P = 89.965 d or 179.930 d.
- Transit depth implies a companion of ≈ 0.9 R_J.
- DR3: RV-stable (rv_chisq_pvalue = 0.88, 22 epochs; per-epoch scatter ≈ 1.44 km/s); RUWE 1.51.

| Companion | P = 89.965 d: a0 | Δχ² | K1 | P = 179.930 d: a0 | Δχ² | K1 |
|---|---|---|---|---|---|---|
| 1 M_J planet | 0.0009 mas | 0.01 | 0.04 km/s | 0.0014 mas | 0.02 | 0.03 km/s |
| 13 M_J | 0.011 mas | 1.3 | 0.55 km/s | 0.018 mas | 3.2 | 0.43 km/s |
| 50 M_J brown dwarf | 0.042 mas | 18 | 2.1 km/s | 0.067 mas | 46 | 1.6 km/s |
| 0.15 M☉ M dwarf | 0.126 mas | 160 | 6.1 km/s | 0.20 mas | 403 | 4.8 km/s |

**DR4 sensitivity at these periods**

| Method | P = 89.965 d | P = 179.930 d |
|---|---|---|
| Astrometry (50% completeness) | 108 M_J | 67 M_J |
| Epoch RVs (3σ) | 22 M_J | 28 M_J |

The RV limits assume the DR3 per-epoch precision and twice the DR3 number of epochs. They are for an edge-on, circular orbit.

**Predictions**
1. **If the companion is a planet (≤ 13 M_J):** no astrometric orbit (Δχ² ≤ 3) and no RV variability at either period above ~0.9 km/s. DR4 then gives a combined upper limit of roughly 22–28 M_J (RVs). DR4 cannot detect the planet itself.
2. **An RV signal at 89.965 or 179.930 d** with K ≳ 1 km/s would mean a brown dwarf or low-mass star.
3. **An astrometric orbit** would mean a companion of ≳ 0.1 M☉ at 90 d (≳ 67 M_J at 180 d). DR3's RV stability already makes a ~0.15 M☉ companion unlikely (K1 ≈ 5–6 km/s against 1.44 km/s scatter).
4. **Our expectation: a null result in both astrometry and RVs.** The photometric transit predicted for 2026-12-28 06:42 UTC (±45 min) remains the decisive test of the period.
   - *Erratum, 2026-10-12, still before DR4:* the source time omitted the 1.2-min barycentric light-travel correction. The transit is BJD_TDB 2461402.7800, i.e. geocentric UTC 2026-12-28 06:41 (±45 min); this conversion was recomputed independently here with astropy. The prediction is otherwise unchanged.

## Addendum, 2026-10-12 (same day, still before DR4): TIC 147035057

Expert feedback on this target raised two points that the hypotheses above did not cover.

**1. Gaia DR3 RV amplitudes are often not confirmed by precise RVs.**
- H1's photocentre size (0.49 mas) assumed the DR3 amplitude is real (K ≈ 25 km/s).
- Add **H0: the DR3 RV amplitude is spurious.** DR4 epoch RVs then show no 37.315-d variation above ~1–2 km/s (the per-epoch scatter).
- The eclipsing companion's mass is then unconstrained by RVs. Any astrometric signal at 37.3 d then measures the mass directly.

**2. DR3 non_single_star = 0 weighs against a double-lined pair (H2), and the field is crowded.**
- Add **H3: a blended or hierarchical eclipsing binary**, where the eclipse arises on a fainter star within ~1–2″, bound or not.
- Its prediction: no significant 37.3-d photocentre orbit of the bright star.
- Possible DR4 signatures of the close neighbour:
  - a resolved or close pair (`nss_resolved_pair`, `duplicated_source`);
  - an acceleration solution if the pair is bound and wide;
  - an elevated excess noise.
- If the DR3 RV variability is real but belongs to the bright star, it would not phase on 37.315 d.

**Revised expectation.** H1 (a real 0.5 mas orbit) remains disfavoured by DR3 RUWE = 1.12. We no longer favour H2 over H0 and H3.

**Discriminating DR4 observables:**

| DR4 observable | Points to |
|---|---|
| 37.3-d astrometric orbit with a0 ~ 0.5 mas | H1 |
| Double-lined epoch spectra / SB2 products | H2 |
| Flat epoch RVs at 37.315 d | H0 (or H3 with a quiet primary) |
| RV variability not phased on 37.315 d | H3, or variability unrelated to the eclipses |

## TIC 238606535 (TOI-5560.01 = EPIC 211924561): added 2026-10-12, before DR4

**Source:** the duo-transits repo, `k2singles/README.md` ("EPIC 211924561: blind re-analysis").

**System**
- Gaia DR3 660640844292494336: G = 10.70, parallax 2.625 mas (381 pc), RUWE 1.19, astrometric_excess_noise 0.21 mas (sig 58), non_single_star = 0.
- FLAME M1 = 1.424 M☉; Teff 6680 K; v sin i 15 km/s (TRES, ExoFOP).
- DR3 values re-queried independently on 2026-10-12 and matched those reported.
- Four transits (K2 C5, K2 C16, TESS S46, TESS S72): P = 140.49401 d or its alias 70.24701 d, T0 = BJD_TDB 2457179.89774.
  - The C16 and S46 epochs fall at cycles 7 and 17 of the 140.494-d period, as computed here.
  - The transit is V-shaped with b ≳ 0.9. Duration vs stellar density needs e ≳ 0.5 at 140.5 d (≳ 0.3 at 70.2 d).
- DR3 RVs: constant (22 transits, p = 0.082; error 0.60 km/s, so ≈ 2.2 km/s per epoch).
- The nearest Gaia neighbour is 8.8″ away at G = 20.5, irrelevant to the photocentre.

| Companion | P = 70.247 d: a0 | Δχ² | K1 | P = 140.494 d: a0 | Δχ² | K1 |
|---|---|---|---|---|---|---|
| 1 M_J planet | 0.0007 mas | 0.003 | 0.04 km/s | 0.0010 mas | 0.008 | 0.03 km/s |
| 13 M_J | 0.0085 mas | 0.6 | 0.50 km/s | 0.0135 mas | 1.4 | 0.40 km/s |
| 50 M_J brown dwarf | 0.032 mas | 8 | 1.9 km/s | 0.051 mas | 20 | 1.5 km/s |
| 0.1 M☉ M dwarf | 0.066 mas | 33 | 3.9 km/s | 0.105 mas | 84 | 3.1 km/s |
| 0.2 M☉ M dwarf | 0.126 mas | 123 | 7.5 km/s | 0.200 mas | 309 | 5.9 km/s |

**DR4 sensitivity at these periods** (edge-on, circular)

| Method | P = 70.247 d | P = 140.494 d |
|---|---|---|
| Astrometry (50% completeness) | 166 M_J | 102 M_J |
| Epoch RVs (3σ; K_lim = 1.41 km/s) | 37 M_J | 47 M_J |

- Assumptions are as for TIC 159614322: DR3 per-epoch precision, twice the DR3 number of epochs.
- An eccentric orbit (e ≈ 0.5) raises K by ~15% at fixed mass, so the RV limits are slightly conservative.
- Numbers come from `scripts/dr4_predictions.py` and `results/dr4_predictions.csv`.

**Predictions**
1. **If the transiting companion is a planet or a low-mass brown dwarf (≲ 40 M_J):**
   - no astrometric orbit at either period (Δχ² ≲ 10);
   - no epoch-RV signal at either period;
   - DR4 then excludes ≳ 37–47 M_J (RVs). It cannot detect a planet.
2. **An RV signal with K ≳ 1.5 km/s at 70.247 or 140.494 d,** phased so that the fitted conjunction falls at the transit epochs (`named_rv.csv` t0_offset_d ≈ 0), would mean a brown dwarf or low-mass star. Its period would also break the alias.
3. **An astrometric orbit** would mean ≳ 0.1 M☉ at 140 d (≳ 0.15 M☉ at 70 d). DR3 RV constancy already disfavours ≳ 0.2 M☉ (K1 ≈ 6–7.5 km/s against ≈ 2.2 km/s per-epoch scatter).
4. **On the DR3 excess noise:** astrometric_excess_noise = 0.21 mas at G = 10.7 is not by itself evidence of a companion (RUWE is only 1.19). If DR4 shows a ~0.2 mas photocentre orbit at one of these periods, it would point to a ~0.2 M☉ companion, in tension with the DR3 RVs.
5. **Our expectation:** a null result in both astrometry and RVs, giving the RV upper limit above. The period alias then remains for photometry to settle.

## Eclipsing binaries from the K2 work: added 2026-10-12, before DR4

**Source:** the duo-transits repo, `k2singles/results/full2_triage.md`; ephemerides in `k2singles/results/ephemerides/<EPIC>.json`. Each binary was vetted there by an adversarial review and a blind re-derivation.

**Independent checks here (2026-10-12)**
- Gaia DR3 values re-queried.
- The DR3 SB1 orbit propagated to eclipse times.
- The HD 220967 identification confirmed.

Numbers come from `scripts/dr4_predictions.py` and `results/dr4_predictions.csv`. The astrometric Δχ² uses the circular-orbit expression, which is an approximation for these eccentric orbits.

### TIC 184895645 (EPIC 211892898): red-giant eclipsing SB1, a pipeline validation case

**System**
- Gaia DR3 660920219027943552: G = 11.81, G_RVS = 10.85, parallax 0.878 mas (1.14 kpc), RUWE 1.60, non_single_star = 2.
- DR3 SB1 solution: P = 72.534 ± 0.029 d, e = 0.480 ± 0.019, ω = 201.5 ± 1.8°, T_p = BJD 2457418.905 ± 0.25, K1 = 22.75 ± 0.85 km/s, so f(M) = 0.060 M☉.
- K2 eclipses (C5, C16, C18): P = 72.50887 ± 0.0014 d, T0 = BJD_TDB 2457189.3272 ± 0.019; depth 6–7 ppt; T14 51–56 h.
- Propagated here, the DR3 elements put the primary conjunction (ν + ω = 90°) 1.0 d from the K2 eclipse ephemeris (≈ 1.5σ in DR3 T_p and P). The secondary conjunction would be 17 d off.
- With M1 = 1.11 M☉ (APO-K2 asteroseismology, as reported by the duo-transits reviewer; not re-checked here) and sin i ≈ 1: M2 = 0.55 M☉.
- The companion contributes ~0.1% of the G flux, so the photocentre follows the giant.

**Predictions**
1. **Epoch RVs:** single-lined (`fraction_double_lined` = 0 or null). The fixed-period Keplerian fit (`named_rv.csv`, kep_*) at 72.509 d should give:
   - K1 = 22.7 ± ~1 km/s and e = 0.48 ± ~0.03;
   - ω ≈ 200°;
   - a primary conjunction within ~1 d of the K2 eclipses (kep_t0_offset_d ≈ 0; |offset| ≲ 1 d).
   The circular-fit offset is not meaningful at this eccentricity.
2. **Astrometry:**
   - photocentre orbit a0 ≈ 0.12 mas at P ≈ 72.5 d, expected Δχ² ≈ 145 (above the 50%-completeness level of ~80), so a detection in the 0.189–0.208 yr window;
   - its inclination should be near 90° (eclipsing), and its period consistent with 72.51 d.
   - A non-detection would most likely mean extra astrometric noise from the giant, not a wrong orbit.
3. **DPAC:** an SB1 or AstroSpectroSB1 solution in the DR4 NSS tables, consistent with DR3.
4. **Use:** if 1–3 hold, the system validates the named-target RV and astrometry chain on a known binary.

### TIC 301328232 (EPIC 246194159 = HD 220967): eccentric EB of two F stars

**System**
- Gaia DR3 2633217515438236544: G = 9.39, parallax 5.624 mas (178 pc), RUWE 1.24, non_single_star = 0.
- DR3 published **no RVS output** for this star (no G_RVS, no RV) despite G = 9.4.
- Identification: the 1.8″ offset from the J2000 position equals 16 yr of its 110 mas/yr proper motion, and it is the only Gaia source within 36″.
- Ephemeris (K2 C12, TESS S70, S92): P = 73.245178 ± 0.000007 d, primary T0 = BJD_TDB 2457747.38501 ± 0.0003.
  - The secondary eclipse falls at phase 0.867, so e cos ω ≈ 0.577 and e ≥ 0.58 (blind estimate e ≈ 0.66).
  - Depths are 37% and 30%, suggesting two similar F stars.

| Mass ratio (M1 = 1.25 M☉) | G flux ratio | a0 | Expected Δχ² | K1 (e = 0.66) | Max line separation |
|---|---|---|---|---|---|
| q = 0.8 | 0.37 | 0.44 mas | ≈ 1,500 | 39 km/s | ≈ 140 km/s |
| q = 0.9 | 0.60 | 0.25 mas | ≈ 500 | 43 km/s | ≈ 143 km/s |
| q = 1.0 | 1.00 | 0 | 0 | 46 km/s | ≈ 145 km/s |

**Predictions**
1. **Epoch RVs:** DR4 places this star in the double-lined products (`fraction_double_lined` > 0). Component separations reach ~140 km/s near periastron, far above the RVS resolution (~25 km/s).
   - Any single-template epoch RVs will scatter strongly.
   - If DR4 publishes double-lined epoch velocities, they give a benchmark dynamical mass for both stars, with e ≈ 0.6 at P = 73.245 d.
2. **Astrometry:**
   - DR3 RUWE = 1.24 at G = 9.4 favours a small photocentre orbit, i.e. a near-equal pair (q ≳ 0.95, a0 ≲ 0.1 mas). An a0 ≳ 0.25 mas would conflict with it.
   - Our expectation: no significant photocentre orbit, or a marginal one, in the 0.190–0.211 yr window.
   - A strong orbit would point to q ≲ 0.9 and a fainter secondary.
3. **Phase check:**
   - If the RV fit is single-lined (blended), K and ω are not physical; the double-lined products decide.
   - The Keplerian fit's e should be ≥ 0.58 if the RVs track the primary.

## Outcomes (to be added after 2 Dec 2026)

*(empty until DR4)*
