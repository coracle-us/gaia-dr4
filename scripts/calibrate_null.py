"""Null calibration of the detection statistics on simulated DR4 single stars.

For N realisations: pick a target from data/targets.csv (its RA, Dec, G), simulate a single star
(no companion) with the DR4 sampling (companions.simulate), add realistic excess noise, and run
the full pipeline statistics:
  dlnl2_orbit  (search_orbit_ml, P = 0.1-4 yr)   dlnl2_acc7 / dlnl2_acc9 (acceleration models)
Excess noise: per-CCD s_ccd ~ U(0, 0.05) mas and per-transit s_tr ~ U(0, 0.05) mas, for all G
(gaiamock adds U(0, 0.04) mas for G < 13; the DR4 pre-release shows per-CCD excess ~0.05 mas at
G ~ 12 and a mild per-transit excess at G ~ 14; notes/phase1_validation.md).
Output: results/calibration_null.csv (one row per realisation; seed = index).
Usage:  python scripts/calibrate_null.py --n 10000 --procs 30
"""
from __future__ import annotations

import argparse
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "calibration_null.csv"
_SL = None
_TARGETS = None


def _init():
    global _SL, _TARGETS
    from companions.scanlaw import ScanLaw
    _SL = ScanLaw.load()
    t = pd.read_csv(ROOT / "data" / "targets.csv").drop_duplicates("hostname")
    _TARGETS = t.dropna(subset=["phot_g_mean_mag", "ra", "dec", "dist_pc"]).reset_index(drop=True)


def run_one(i: int) -> dict:
    from companions.accel import fit_acceleration_ml
    from companions.search import search_orbit_ml
    from companions.simulate import simulate_epochs
    rng = np.random.default_rng(i)
    row = _TARGETS.iloc[rng.integers(len(_TARGETS))]
    s_ccd, s_tr = rng.uniform(0, 0.05, 2)
    star = [0.0, 0.0, 1000.0 / row["dist_pc"], rng.normal(0, 30), rng.normal(0, 30)]
    ep = simulate_epochs(_SL, row["ra"], row["dec"], row["phot_g_mean_mag"], star, rng,
                         s_ccd=s_ccd, s_tr=s_tr)
    orb = search_orbit_ml(ep, model="5p")
    a7 = fit_acceleration_ml(ep, order=7, model="5p")
    a9 = fit_acceleration_ml(ep, order=9, model="5p")
    return dict(seed=i, hostname=row["hostname"], g_mag=row["phot_g_mean_mag"],
                n_obs=len(ep.w), n_transit=int(ep.transit.max() + 1),
                s_ccd_true=s_ccd, s_tr_true=s_tr, s_ccd_fit=orb.star.s_ccd, s_tr_fit=orb.star.s_tr,
                dlnl2_orbit=orb.dlnl2, period=orb.period, e=orb.e,
                dlnl2_acc7=a7.dlnl2, dlnl2_acc9=a9.dlnl2)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=10000)
    ap.add_argument("--procs", type=int, default=30)
    ap.add_argument("--start", type=int, default=0)
    args = ap.parse_args()
    t0 = time.time()
    with Pool(args.procs, initializer=_init) as pool:
        rows = pool.map(run_one, range(args.start, args.start + args.n), chunksize=8)
    df = pd.DataFrame(rows)
    OUT.parent.mkdir(exist_ok=True)
    df.to_csv(OUT, index=False, float_format="%.6g")
    print(f"{len(df)} nulls in {time.time() - t0:.0f} s -> {OUT}")
    for col in ("dlnl2_orbit", "dlnl2_acc7", "dlnl2_acc9"):
        q = np.percentile(df[col], [50, 90, 99, 99.9])
        print(f"{col:12s} median {q[0]:6.1f}  90% {q[1]:6.1f}  99% {q[2]:6.1f}  99.9% {q[3]:6.1f}  "
              f"max {df[col].max():6.1f}")


if __name__ == "__main__":
    main()
