"""Gaia per-CCD along-scan precision model.

Holl et al. 2023 (A&A 674, A10; arXiv:2206.05439) Fig. 3, "EDR3 adjusted" curve, per CCD, as
digitised in El-Badry's gaiamock (al_uncertainty_per_ccd_interp). Used by Lammers & Winn 2025.
The DR4 pre-release centroid errors are equal or smaller at G < 12 (notes/dr4_data_model.md),
so this is conservative at the bright end.
"""
import numpy as np

_G_GRID = [4, 5, 6, 7, 8.2, 8.4, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
_SIG_CCD = [0.4, 0.35, 0.15, 0.17, 0.23, 0.13, 0.13, 0.135, 0.125, 0.13, 0.15, 0.23,
            0.36, 0.63, 1.05, 2.05, 4.1]


def sigma_ccd_holl2023(g_mag):
    """Per-CCD AL uncertainty [mas] at G magnitude g_mag."""
    return np.interp(g_mag, _G_GRID, _SIG_CCD)
