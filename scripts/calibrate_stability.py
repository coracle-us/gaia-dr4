"""Calibrate the drop-one-fifth stability check on genuine vs spurious orbit detections.

Genuine: injections from results/injection_recovery.csv with dlnl2_orbit in [threshold, 120],
P_true < 3.8 yr and correctly recovered period (regenerated from their seeds).
Spurious: the highest-dlnl2 null realisations from results/calibration_null.csv (regenerated
from their seeds); these are the noise configurations most likely to pass detection.
For each: the full orbit search and the five searches with one time-ordered fifth of the
transits removed -> min(dlnl2_sub / dlnl2_full) and max |P_sub / P_full - 1|.
Output: results/stability_calibration.csv; a criterion summary is printed.
"""
from __future__ import annotations

import argparse
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
_SL = None
_TARGETS = None


def _init():
    global _SL, _TARGETS
    import calibrate_null
    calibrate_null._init()
    _SL, _TARGETS = calibrate_null._SL, calibrate_null._TARGETS


def _null_epochs(seed):
    from companions.simulate import simulate_epochs
    rng = np.random.default_rng(seed)                     # same draws as calibrate_null.run_one
    row = _TARGETS.iloc[rng.integers(len(_TARGETS))]
    s_ccd, s_tr = rng.uniform(0, 0.05, 2)
    star = [0.0, 0.0, 1000.0 / row["dist_pc"], rng.normal(0, 30), rng.normal(0, 30)]
    return simulate_epochs(_SL, row["ra"], row["dec"], row["phot_g_mean_mag"], star, rng,
                           s_ccd=s_ccd, s_tr=s_tr)


def _injection_epochs(seed):
    """Same draws as injection_recovery.run_one (seed includes SEED_OFFSET)."""
    from companions.precision import sigma_ccd_holl2023
    from companions.simulate import N_CCD, simulate_epochs
    rng = np.random.default_rng(seed)
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
    star = [0.0, 0.0, 1000.0 / row["dist_pc"], rng.normal(0, 30), rng.normal(0, 30)]
    return simulate_epochs(_SL, row["ra"], row["dec"], g, star, rng, orbit=orbit,
                           s_ccd=s_ccd, s_tr=s_tr)


def _run(job):
    from companions.cascade import N_STABILITY_BLOCKS, subset_epochs
    from companions.search import search_orbit_ml
    kind, seed = job
    ep = _null_epochs(seed) if kind == "null" else _injection_epochs(seed)
    full = search_orbit_ml(ep, model="5p")
    tt = np.array([ep.t[ep.transit == k].mean() for k in range(ep.transit.max() + 1)])
    edges = np.quantile(tt, np.linspace(0, 1, N_STABILITY_BLOCKS + 1))
    blk = np.clip(np.searchsorted(edges, tt, side="right") - 1, 0, N_STABILITY_BLOCKS - 1)
    ratios, dps = [], []
    for b in range(N_STABILITY_BLOCKS):
        s = search_orbit_ml(subset_epochs(ep, blk[ep.transit] != b), model="5p")
        ratios.append(s.dlnl2 / full.dlnl2)
        dps.append(abs(s.period / full.period - 1))
    return dict(kind=kind, seed=seed, dlnl2=full.dlnl2, period=full.period, e=full.e,
                min_ratio=min(ratios), second_min_ratio=sorted(ratios)[1], max_dp=max(dps),
                n_dp_gt10=int(np.sum(np.array(dps) > 0.1)), n_dp_gt25=int(np.sum(np.array(dps) > 0.25)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--procs", type=int, default=30)
    ap.add_argument("--n-null", type=int, default=60)
    ap.add_argument("--n-inj", type=int, default=120)
    args = ap.parse_args()
    thr = pd.read_csv(ROOT / "results/thresholds.csv")
    t_orb = float(thr[(thr.statistic == "dlnl2_orbit") & (thr.fap_per_star == 1e-4)].threshold.iloc[0])
    null = pd.read_csv(ROOT / "results/calibration_null.csv").nlargest(args.n_null, "dlnl2_orbit")
    inj = pd.read_csv(ROOT / "results/injection_recovery.csv")
    inj = inj[(inj.dlnl2_orbit > t_orb) & (inj.dlnl2_orbit < 120) & (inj.period_true < 3.8)
              & ((inj.period_fit / inj.period_true - 1).abs() < 0.1)].head(args.n_inj)
    jobs = [("null", int(s)) for s in null.seed] + [("inj", int(s)) for s in inj.seed]
    with Pool(args.procs, initializer=_init) as pool:
        df = pd.DataFrame(pool.map(_run, jobs, chunksize=1))
    df.to_csv(ROOT / "results/stability_calibration.csv", index=False, float_format="%.5g")
    for name, ok in [
        ("current: min_ratio>0.5 & max_dp<0.1", (df.min_ratio > 0.5) & (df.max_dp < 0.1)),
        ("min_ratio>0.4 & n_dp>25% == 0", (df.min_ratio > 0.4) & (df.n_dp_gt25 == 0)),
        ("min_ratio>0.4 & n_dp>10% <= 1", (df.min_ratio > 0.4) & (df.n_dp_gt10 <= 1)),
        ("second_min_ratio>0.5 & n_dp>10% <= 1", (df.second_min_ratio > 0.5) & (df.n_dp_gt10 <= 1)),
        ("n_dp>10% <= 1", df.n_dp_gt10 <= 1),
    ]:
        print(f"{name:40s}  genuine pass {ok[df.kind == 'inj'].mean():.2f}   "
              f"null pass {ok[df.kind == 'null'].mean():.2f}")


if __name__ == "__main__":
    main()
