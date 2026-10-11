"""Per-CCD along-scan epoch astrometry for one source, from the DR4 `epoch_astrometry` table.

Conventions follow ESA's gaiasupdate (and the DR4 data model, section 3.1):
* t: Julian years from J2017.5 (TCB), using the barycentric-corrected observation time.
  obs_time_tcb counts ns from JD 2455197.5 (TCB), called "J2010.0" in the data model.
* Only CCDs with used_by_agis_al are kept (this also drops the SkyMapper CCD).
* w is centroid_pos_al [mas], relative to the per-source reference point (ra0, dec0).
* sigma is centroid_pos_error_al [mas] (IPD + calibration noise, *no* source excess noise).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

NS_PER_DAY = 86400e9
JD_TIME_ORIGIN = 2455197.5        # obs_time_tcb zero point (TCB)
JD_REF_EPOCH = 2457936.875        # J2017.5 (TCB), DR4 reference epoch
DAYS_PER_YEAR = 365.25


def _flat(column, fill, dtype):
    """Concatenate a per-CCD array column (one masked array per FoV-transit row)."""
    return np.concatenate([np.ma.filled(np.ma.asarray(x).astype(dtype), fill) for x in column])


@dataclass
class EpochData:
    source_id: int
    t: np.ndarray            # yr from J2017.5 (TCB, barycentric)
    psi: np.ndarray          # scan position angle [rad]
    plx_factor: np.ndarray   # parallax factor AL
    colour_factor: np.ndarray
    w: np.ndarray            # centroid_pos_al [mas]
    sigma: np.ndarray        # centroid_pos_error_al [mas]
    transit: np.ndarray      # 0..n_transit-1, groups CCDs of one FoV transit
    agis_source_excess_noise: float  # [mas], unscaled AGIS value

    @classmethod
    def from_table(cls, table, source_id: int) -> "EpochData":
        rows = table[table["source_id"] == source_id]
        if len(rows) == 0:
            raise ValueError(f"source_id {source_id} not in table")
        n_ccd = np.array([len(np.ma.asarray(x)) for x in rows["obs_time_tcb"]])
        per_transit = lambda col: np.repeat(np.ma.filled(rows[col].astype(float), np.nan), n_ccd)

        t_ns = _flat(rows["obs_time_tcb"], np.nan, float) + per_transit("obs_time_bary_corr")
        jd = JD_TIME_ORIGIN + t_ns / NS_PER_DAY
        used = _flat(rows["used_by_agis_al"], False, bool)
        transit_all = np.repeat(np.arange(len(rows)), n_ccd)

        keep = used & np.isfinite(jd)
        _, transit = np.unique(transit_all[keep], return_inverse=True)
        return cls(
            source_id=int(source_id),
            t=((jd - JD_REF_EPOCH) / DAYS_PER_YEAR)[keep],
            psi=np.deg2rad(_flat(rows["scan_pos_angle"], np.nan, float))[keep],
            plx_factor=per_transit("parallax_factor_al")[keep],
            colour_factor=_flat(rows["colour_factor_al"], np.nan, float)[keep],
            w=_flat(rows["centroid_pos_al"], np.nan, float)[keep],
            sigma=_flat(rows["centroid_pos_error_al"], np.nan, float)[keep],
            transit=transit,
            agis_source_excess_noise=float(rows["agis_source_excess_noise"][0]),
        )
