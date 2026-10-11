"""Simulated DR4 epoch astrometry for injection-recovery and null calibration.

Sampling: predicted DR4 astrometric transits (companions.scanlaw), a random fraction `reject`
dropped (gaiamock and Lammers & Winn use 10%), 9 AF CCDs per transit spaced 4.85 s, scan angle
constant within a transit, parallax factor from Gaia's barycentric position.
Noise: per-CCD sigma from the Holl+2023 curve (or sigma_override), plus optional per-CCD jitter
s_ccd and per-transit correlated jitter s_tr. Signal: 5-parameter single star (+ Keplerian
photocentre orbit, Gaia DR3 Thiele-Innes convention).
"""
from __future__ import annotations

import numpy as np

from companions.epoch import DAYS_PER_YEAR, EpochData
from companions.linear import single_star_design
from companions.orbit import orbit_design, thiele_innes
from companions.precision import sigma_ccd_holl2023
from companions.scanlaw import parallax_factor_al

N_CCD = 9
CCD_SPACING_S = 4.85


def simulate_epochs(scanlaw, ra, dec, g_mag, star_params, rng, orbit=None, reject=0.1,
                    s_ccd=0.0, s_tr=0.0, sigma_override=None, noise=True,
                    source_id=0) -> EpochData:
    t_tr, psi_tr = scanlaw.transits(ra, dec)
    keep = rng.uniform(size=len(t_tr)) >= reject
    t_tr, psi_tr = t_tr[keep], psi_tr[keep]
    n_tr = len(t_tr)
    offsets = (np.arange(N_CCD) - (N_CCD - 1) / 2) * CCD_SPACING_S / 86400 / DAYS_PER_YEAR
    transit = np.repeat(np.arange(n_tr), N_CCD)
    t = (t_tr[:, None] + offsets[None, :]).ravel()
    psi = np.repeat(psi_tr, N_CCD)
    plx = np.repeat(parallax_factor_al(t_tr, psi_tr, ra, dec), N_CCD)
    sigma = np.full(len(t), sigma_ccd_holl2023(g_mag) if sigma_override is None else sigma_override)
    ep = EpochData(source_id=source_id, t=t, psi=psi, plx_factor=plx,
                   colour_factor=np.zeros(len(t)), w=np.zeros(len(t)), sigma=sigma,
                   transit=transit, agis_source_excess_noise=0.0)
    w = single_star_design(ep, "5p") @ np.asarray(star_params, float)
    if orbit is not None:
        ti = np.array(thiele_innes(orbit["a0"], orbit["inc"], orbit["omega"], orbit["big_omega"]))
        w = w + orbit_design(ep, orbit["period"], orbit["e"], orbit["tp"]) @ ti
    if noise:
        w = (w + np.sqrt(sigma**2 + s_ccd**2) * rng.standard_normal(len(t))
             + s_tr * rng.standard_normal(n_tr)[transit])
    ep.w = w
    return ep
