"""Pre-registered Gaia DR4 predictions for the two named TESS targets (written before 2 Dec 2026).

For each hypothesis: the photocentre semi-major axis a0, the expected full-orbit Delta-chi2 of the
DR4 epoch astrometry (dchi2_exp = N_obs a0^2 (1 + cos^2 i) / (4 sigma_eff^2), sigma_eff^2 =
sigma_ccd(G)^2 + 0.025^2; N_obs = 0.9 x 9 x predicted DR4 transits; the calibrated detection
threshold and 50% completeness are dlnl2 > 42.6 and dchi2_exp ~ 80, results/thresholds.csv and
results/completeness_models.csv), and the primary's RV semi-amplitude K1 compared with the Gaia
per-epoch RV scatter inferred from DR3 (sqrt((err^2 - 0.11^2) 2N/pi), DR3 documentation).
Stellar parameters, periods and DR3 values: McBee 2026 RNAAS (submitted) and its write-ups.
Output: results/dr4_predictions.csv (printed). The narrative is notes/dr4_predictions.md.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from companions.impostor import MainSequence, photocentre_factor
from companions.precision import sigma_ccd_holl2023
from companions.scanlaw import ScanLaw
from companions.simulate import N_CCD

ROOT = Path(__file__).resolve().parents[1]
MJUP = 9.5458e-4
AU_YR_KMS = 4.740470463
THRESHOLD_DLNL2, D50 = 42.6, 80.0
S_EXCESS = 0.025


def n_obs(sl, ra, dec):
    return 0.9 * N_CCD * len(sl.transits(ra, dec)[0])


def orbit_terms(m1, m2, p_days, plx, flux_ratio=0.0):
    p = p_days / 365.25
    a_rel = ((m1 + m2) * p**2) ** (1 / 3)
    q = m2 / m1
    a0 = a_rel * abs(photocentre_factor(q, flux_ratio)) * plx
    k1 = 2 * np.pi * a_rel * (m2 / (m1 + m2)) / p * AU_YR_KMS          # km/s, sin i = 1
    return a_rel, a0, k1


def dchi2(n, a0, g, inc_deg=90.0):
    s2 = sigma_ccd_holl2023(g) ** 2 + S_EXCESS**2
    return n * a0**2 * (1 + np.cos(np.radians(inc_deg)) ** 2) / (4 * s2)


def min_mass_function_m2(fm, m1):
    """M2 for mass function fm with sin i = 1."""
    lo, hi = 1e-4, 10.0
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if mid**3 / (m1 + mid) ** 2 < fm else (lo, mid)
    return 0.5 * (lo + hi)


def mass_at_d50(n, g, m1, p_days, plx):
    """Dark-companion mass [M_J] with dchi2_exp = D50 (50% DR4 completeness), edge-on."""
    s2 = sigma_ccd_holl2023(g) ** 2 + S_EXCESS**2
    a0_req = np.sqrt(4 * s2 * D50 / n)
    lo, hi = 1e-3, 5000.0
    for _ in range(100):
        mid = np.sqrt(lo * hi)
        lo, hi = (mid, hi) if orbit_terms(m1, mid * MJUP, p_days, plx)[1] < a0_req else (lo, mid)
    return np.sqrt(lo * hi)


def main() -> None:
    sl, ms = ScanLaw.load(), MainSequence.load()
    rows = []

    # --- TIC 147035057: G = 8.83, 166 pc (plx 6.023), M1 = 1.44, P = 37.31497 d, grazing eclipse
    g, plx, m1, p = 8.83, 6.023, 1.44, 37.31497
    n = n_obs(sl, 104.27440, -16.47848)
    k_obs = 49.8 / 2                                                     # DR3 rv_amplitude_robust / 2
    fm = (k_obs / AU_YR_KMS) ** 3 * (p / 365.25) / (2 * np.pi) ** 3      # mass function [M_sun]
    m2 = min_mass_function_m2(fm, m1)
    f_g = ms.flux_ratio(m1, m2 / m1, "G")
    a_rel, a0, k1 = orbit_terms(m1, m2, p, plx, f_g)
    inc = 90 - np.degrees(np.arcsin(min(1.0, 1.28 * 2.0 * 0.00465 / a_rel)))   # b = 1.28, R1 = 2 R_sun
    rows.append(dict(target="TIC 147035057", hypothesis="H1 single-lined star (from DR3 RV amplitude)",
                     period_d=p, m2_msun=m2, flux_ratio_g=f_g, a_rel_au=a_rel, a0_mas=a0,
                     inc_deg=inc, n_obs=n, dchi2_exp=dchi2(n, a0, g, inc), k1_kms=k1))
    m2_h1, fm_h1 = m2, fm
    m_twin = 1.2                                                         # near-equal pair
    for q in (0.97, 0.99, 1.0):
        f_t = ms.flux_ratio(m_twin, q, "G")
        a_rel_t, a0_t, k1_t = orbit_terms(m_twin, q * m_twin, p, plx, f_t)
        rows.append(dict(target="TIC 147035057", hypothesis=f"H2 near-equal SB2 pair (q = {q})",
                         period_d=p, m2_msun=q * m_twin, flux_ratio_g=f_t, a_rel_au=a_rel_t,
                         a0_mas=a0_t, inc_deg=90.0, n_obs=n, dchi2_exp=dchi2(n, a0_t, g),
                         k1_kms=k1_t, line_split_kms=k1_t * (1 + q) / q))

    # --- TIC 159614322: G = 11.63, 401 pc (plx 2.495), M1 = 1.11; P = 89.965 or 179.930 d
    g, plx, m1 = 11.63, 2.495, 1.11
    n = n_obs(sl, 220.00977, 65.48404)
    sig_rv = np.sqrt((0.4**2 - 0.11**2) * 2 * 22 / np.pi)               # DR3: err 0.4 km/s, 22 transits
    for p in (89.9651, 179.9302):
        for label, m2 in (("H1 planet, 1 M_J", 1 * MJUP), ("H1 planet, 13 M_J", 13 * MJUP),
                          ("H2 brown dwarf, 50 M_J", 50 * MJUP), ("H3 M dwarf, 0.15 M_sun", 0.15)):
            f2 = ms.flux_ratio(m1, m2 / m1, "G") if m2 > 0.08 else 0.0
            a_rel, a0, k1 = orbit_terms(m1, m2, p, plx, f2)
            rows.append(dict(target="TIC 159614322", hypothesis=label, period_d=p, m2_msun=m2,
                             flux_ratio_g=f2, a_rel_au=a_rel, a0_mas=a0, inc_deg=90.0, n_obs=n,
                             dchi2_exp=dchi2(n, a0, g), k1_kms=k1, gaia_rv_sigma_epoch_kms=sig_rv))
        rows.append(dict(target="TIC 159614322", hypothesis="DR4 astrometry 50%-completeness mass",
                         period_d=p, m2_msun=mass_at_d50(n, g, m1, p, plx) * MJUP, n_obs=n))
        # DR4 epoch-RV limit: assume DR3 per-epoch scatter and twice the DR3 epochs (66 vs 34
        # months), 3-sigma on K for a circular, edge-on orbit: K_lim = 3 sigma sqrt(2 / N)
        k_lim = 3 * sig_rv * np.sqrt(2 / 44)
        lo, hi = 1e-3, 1000.0
        for _ in range(100):
            mid = np.sqrt(lo * hi)
            lo, hi = (mid, hi) if orbit_terms(m1, mid * MJUP, p, plx)[2] < k_lim else (lo, mid)
        rows.append(dict(target="TIC 159614322", hypothesis="DR4 epoch-RV 3-sigma mass limit",
                         period_d=p, m2_msun=np.sqrt(lo * hi) * MJUP, k1_kms=k_lim,
                         gaia_rv_sigma_epoch_kms=sig_rv))
    df = pd.DataFrame(rows)
    df["m2_mjup"] = df.m2_msun / MJUP
    df["detectable_astrometry"] = df.dchi2_exp > D50
    df.to_csv(ROOT / "results/dr4_predictions.csv", index=False, float_format="%.5g")
    cols = ["target", "hypothesis", "period_d", "m2_mjup", "a0_mas", "inc_deg", "dchi2_exp",
            "k1_kms", "line_split_kms"]
    with pd.option_context("display.width", 220):
        print(df[cols].round(4).to_string(index=False))
    print(f"TIC 147035057 mass function f(M) = {fm_h1:.4f} M_sun -> M2 >= {m2_h1:.2f} M_sun (sin i = 1)")
    print(f"TIC 159614322 Gaia per-epoch RV scatter (from DR3) ~ {sig_rv:.2f} km/s")


if __name__ == "__main__":
    main()
