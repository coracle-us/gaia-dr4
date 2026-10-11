"""Gaia DR4 epoch radial velocities (rvs_epoch_parameters_single) for impostor vetting.

Field names from the DR4 draft data model, section 17.2: obs_time_rv [BJD(TCB) - 2455197.5, d],
radial_velocity, radial_velocity_error [km/s], rv_validity_flag. Epoch RVs exist for
G_RVS <= 14. At ~1 km/s per epoch they cannot see planets or twin impostors (blended amplitude
~0.1 km/s; notes/phase1_validation.md) but do reveal non-twin stellar companions.
rv_constancy: chi-square test of a constant RV after subtracting known inner (transiting)
planets on circular orbits, v_b = -K_b sin(2 pi (t - T_c) / P_b).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.stats import chi2 as chi2_dist

from companions.epoch import DAYS_PER_YEAR, JD_REF_EPOCH, JD_TIME_ORIGIN
from companions.joint import RVData


def rvdata_from_gaia_epoch(table, source_id: int) -> RVData:
    rows = table[np.asarray(table["source_id"]) == source_id]
    rv = np.ma.filled(np.ma.asarray(rows["radial_velocity"], dtype=float), np.nan)
    err = np.ma.filled(np.ma.asarray(rows["radial_velocity_error"], dtype=float), np.nan)
    t_d = np.ma.filled(np.ma.asarray(rows["obs_time_rv"], dtype=float), np.nan)
    ok = np.isfinite(rv) & np.isfinite(err) & np.isfinite(t_d)
    if "rv_validity_flag" in rows.colnames:
        ok &= np.ma.filled(np.ma.asarray(rows["rv_validity_flag"]), False).astype(bool)
    t = (t_d[ok] + JD_TIME_ORIGIN - JD_REF_EPOCH) / DAYS_PER_YEAR
    return RVData(t=t, rv=rv[ok] * 1000.0, err=err[ok] * 1000.0,
                  inst=np.zeros(ok.sum(), int), inst_names=["Gaia RVS"])


@dataclass
class ConstancyTest:
    chi2: float
    dof: int
    pvalue: float
    rms_ms: float


def rv_constancy(rv: RVData, inner=(), floor_ms: float = 0.0) -> ConstancyTest:
    """inner: [(P_b [yr], T_c [yr from J2017.5], K_b [m/s]), ...] subtracted before the test."""
    y = rv.rv.copy()
    for p_b, tc_b, k_b in inner:
        y += k_b * np.sin(2 * np.pi * (rv.t - tc_b) / p_b)
    var = rv.err**2 + floor_ms**2
    mean = np.sum(y / var) / np.sum(1 / var)
    chi2 = float(np.sum((y - mean) ** 2 / var))
    dof = len(y) - 1
    return ConstancyTest(chi2, dof, float(chi2_dist.sf(chi2, dof)) if dof > 0 else np.nan,
                         float(np.std(y - mean)))
