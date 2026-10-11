"""Write epoch astrometry in the DR4 `epoch_astrometry` (DataLink RAW) layout.

Inverse of companions.epoch.EpochData.from_table: one row per FoV transit, per-CCD arrays of
length 10 (slot 0 = SkyMapper, 1-9 = AF1-AF9). CCDs present in the EpochData go to AF slots with
used_by_agis_al = True; empty slots carry the transit's values with used_by_agis_al = False.
obs_time_tcb is written in ns from JD 2455197.5 (TCB) with obs_time_bary_corr = 0, i.e. the
stored time is already barycentric. Used by the simulated archive source and for caching.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
from astropy.io.votable import from_table
from astropy.table import Table

from companions.epoch import DAYS_PER_YEAR, JD_REF_EPOCH, JD_TIME_ORIGIN, NS_PER_DAY

N_SLOTS = 10


def epochs_to_dr4_table(ep, ra0: float, dec0: float) -> Table:
    n_tr = int(ep.transit.max()) + 1
    t_ns = np.rint((ep.t * DAYS_PER_YEAR + (JD_REF_EPOCH - JD_TIME_ORIGIN)) * NS_PER_DAY).astype(np.int64)
    shape = (n_tr, N_SLOTS)
    obs = np.zeros(shape, np.int64)
    psi = np.zeros(shape)
    colour = np.zeros(shape, np.float32)
    w = np.zeros(shape)
    err = np.ones(shape, np.float32)
    used = np.zeros(shape, bool)
    plx = np.zeros(n_tr, np.float32)
    for k in range(n_tr):
        idx = np.where(ep.transit == k)[0]
        if len(idx) > N_SLOTS - 1:
            raise ValueError(f"transit {k} has {len(idx)} CCDs (max {N_SLOTS - 1})")
        slots = np.arange(1, len(idx) + 1)
        obs[k, :] = t_ns[idx[0]]
        psi[k, :] = np.rad2deg(ep.psi[idx[0]])
        obs[k, slots], psi[k, slots] = t_ns[idx], np.rad2deg(ep.psi[idx])
        colour[k, slots], w[k, slots], err[k, slots] = ep.colour_factor[idx], ep.w[idx], ep.sigma[idx]
        used[k, slots] = True
        plx[k] = ep.plx_factor[idx[0]]
    return Table({
        "source_id": np.full(n_tr, ep.source_id, np.int64),
        "transit_id": np.arange(n_tr, dtype=np.int64),
        "ra0": np.full(n_tr, ra0), "dec0": np.full(n_tr, dec0),
        "agis_source_excess_noise": np.full(n_tr, ep.agis_source_excess_noise, np.float32),
        "obs_time_tcb": obs, "obs_time_bary_corr": np.zeros(n_tr, np.float32),
        "scan_pos_angle": psi, "parallax_factor_al": plx, "colour_factor_al": colour,
        "centroid_pos_al": w, "centroid_pos_error_al": err, "used_by_agis_al": used,
    })


def write_votable(table: Table, path: Path) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    from_table(table).to_xml(str(path), tabledata_format="binary2")
