"""Injection-recovery on simulated DR4 epoch astrometry of the target hosts.

Each injection: a random host from data/targets.csv (RA, Dec, G, distance, M_star), a dark
companion on a Keplerian orbit with
  P ~ log-uniform(0.2, 8) yr, e ~ U(0, 0.6), omega, Omega, T_p uniform,
  inclination: half isotropic (cos i ~ U(-1, 1)), half aligned with a transiting planet
               (i ~ N(90, 3) deg -- the science case, and the hardest geometry),
  amplitude set by a target expected Delta-chi2  dchi2_exp ~ log-uniform(3, 3000), where
  dchi2_exp = N_obs a0^2 (1 + cos^2 i) / (4 sigma_eff^2),  sigma_eff^2 = sigma_ccd^2 + s_ccd^2
  (i.e. the full-orbit signal before absorption by the 5-parameter fit),
plus the same excess-noise prior as the null calibration (s_ccd, s_tr ~ U(0, 0.05) mas).
The full pipeline statistics are recorded (orbit search, 7p/9p accelerations) with the fitted
orbit. Output: results/injection_recovery.csv. Usage: python scripts/injection_recovery.py --n 6000
"""
from __future__ import annotations

import argparse
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "injection_recovery.csv"
MJUP_MSUN = 9.5458e-4
SEED_OFFSET = 1_000_000          # keep injection seeds disjoint from the null calibration
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
    from companions.orbit import campbell_from_thiele_innes
    from companions.posterior import companion_mass_mj
    from companions.precision import sigma_ccd_holl2023
    from companions.search import search_orbit_ml
    from companions.simulate import N_CCD, simulate_epochs

    rng = np.random.default_rng(SEED_OFFSET + i)
    row = _TARGETS.iloc[rng.integers(len(_TARGETS))]
    s_ccd, s_tr = rng.uniform(0, 0.05, 2)
    aligned = rng.uniform() < 0.5
    inc = np.deg2rad(rng.normal(90, 3)) if aligned else np.arccos(rng.uniform(-1, 1))
    period = np.exp(rng.uniform(np.log(0.2), np.log(8.0)))
    e = rng.uniform(0, 0.6)
    dchi2_exp = np.exp(rng.uniform(np.log(3), np.log(3000)))
    g = row["phot_g_mean_mag"]
    sig_eff2 = sigma_ccd_holl2023(g) ** 2 + s_ccd**2
    n_obs_pred = 0.9 * N_CCD * len(_SL.transits(row["ra"], row["dec"])[0])
    a0 = np.sqrt(4 * sig_eff2 * dchi2_exp / (n_obs_pred * (1 + np.cos(inc) ** 2)))
    orbit = dict(period=period, e=e, tp=rng.uniform(0, period), a0=a0, inc=inc,
                 omega=rng.uniform(0, 2 * np.pi), big_omega=rng.uniform(0, np.pi))
    plx = 1000.0 / row["dist_pc"]
    star = [0.0, 0.0, plx, rng.normal(0, 30), rng.normal(0, 30)]
    ep = simulate_epochs(_SL, row["ra"], row["dec"], g, star, rng, orbit=orbit,
                         s_ccd=s_ccd, s_tr=s_tr)
    orb = search_orbit_ml(ep, model="5p")
    a7 = fit_acceleration_ml(ep, order=7, model="5p")
    a9 = fit_acceleration_ml(ep, order=9, model="5p")
    a0_fit, inc_fit, _, _ = campbell_from_thiele_innes(*orb.fit.params[-4:])
    mstar = row["mstar_msun"] if np.isfinite(row["mstar_msun"]) else 1.0
    return dict(seed=SEED_OFFSET + i, hostname=row["hostname"], g_mag=g, dist_pc=row["dist_pc"],
                mstar=mstar, n_obs=len(ep.w), s_ccd_true=s_ccd, s_tr_true=s_tr,
                aligned=aligned, inc_true=inc, period_true=period, e_true=e, a0_true=a0,
                mass_true_mj=float(companion_mass_mj(a0, plx, period, mstar)),
                dchi2_exp=dchi2_exp,
                dchi2_exp_actual=len(ep.w) * a0**2 * (1 + np.cos(inc) ** 2) / (4 * sig_eff2),
                dlnl2_orbit=orb.dlnl2, period_fit=orb.period, e_fit=orb.e, a0_fit=a0_fit,
                inc_fit=inc_fit, dlnl2_acc7=a7.dlnl2, dlnl2_acc9=a9.dlnl2)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=6000)
    ap.add_argument("--procs", type=int, default=30)
    args = ap.parse_args()
    t0 = time.time()
    with Pool(args.procs, initializer=_init) as pool:
        rows = pool.map(run_one, range(args.n), chunksize=8)
    df = pd.DataFrame(rows)
    OUT.parent.mkdir(exist_ok=True)
    df.to_csv(OUT, index=False, float_format="%.6g")
    print(f"{len(df)} injections in {time.time() - t0:.0f} s -> {OUT}")


if __name__ == "__main__":
    main()
