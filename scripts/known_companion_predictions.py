"""Pre-registered Gaia DR4 detectability of KNOWN outer companions to transiting giant-planet hosts.

For every companion in data/prediction_companions.csv (literature RV / RV+HG companions), simulate
DR4 epoch astrometry of its host n times (companions.simulate: GOST scan law, Holl+2023 per-CCD
noise, 10% transits rejected, excess noise s_ccd, s_tr ~ U(0, 0.05) mas as in the calibrations),
inject the companion, and run the release-day statistics (orbit search, 7p/9p accelerations).
Detection = statistic above its FAP 1e-4 per-star threshold (results/thresholds.csv).

Orbit draws per realization:
  * period, e from the table (e ~ U(0, 0.5) when not given); omega, T_p (phase) uniform: the
    published periastron epochs are not used yet, so phase is marginalised (matters for P > 5.5 yr);
  * mass_type msini: inclination from a Rayleigh(16 deg) mutual inclination about the transiting
    planet (Xiao+2026 population width) with a uniform node, i.e. cos i = sin(dI) cos(phi);
    true mass M = msini / sin i;
  * mass_type true: M fixed; sin i = msini/M when msini is given (sense random), else the same
    aligned draw as above.
Output: results/known_companion_predictions.csv (summary) and ..._draws.csv.
Usage: python scripts/known_companion_predictions.py --n 200 --procs 4
"""
from __future__ import annotations

import argparse
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
MJUP_MSUN = 9.5458e-4
SIGMA_DI_DEG = 16.0
SEED = 7_000_000
COMPANIONS = "prediction_companions_v2.csv"
_SL = None
_ROWS = None


def _init(rows):
    global _SL, _ROWS
    from companions.scanlaw import ScanLaw
    _SL = ScanLaw.load()
    _ROWS = rows


def load_rows():
    c = pd.read_csv(ROOT / "data" / COMPANIONS)
    for col in ("omega_deg", "tp_bjd", "tp_err_d"):
        if col not in c:
            c[col] = np.nan
    t = pd.read_csv(ROOT / "data/targets.csv").drop_duplicates("hostname")
    df = c.merge(t[["hostname", "source_id", "ra", "dec", "phot_g_mean_mag", "dist_pc", "mstar_msun",
                    "ruwe", "epoch_rv_expected"]], on="hostname", how="left")
    miss = df[df.ra.isna()].hostname.tolist()
    if miss:
        print("not in targets.csv:", miss)
    df = df[df.ra.notna()].reset_index(drop=True)
    m = df.mstar_msun.fillna(1.0)
    # fill P from a (Kepler III, companion mass neglected) and a from P
    p_yr = df.period_d / 365.25
    a = df.a_au.copy()
    a[a.isna()] = (m * p_yr**2) ** (1 / 3)
    p_yr[p_yr.isna()] = np.sqrt(df.a_au**3 / m)
    df["period_yr"] = p_yr
    df["a_rel_au"] = a
    return df


def run_one(job):
    from companions.accel import fit_acceleration_ml
    from companions.epoch import DAYS_PER_YEAR, JD_REF_EPOCH
    from companions.precision import sigma_ccd_holl2023
    from companions.search import search_orbit_ml
    from companions.simulate import N_CCD, simulate_epochs

    k, j, use_phase = job
    r = _ROWS.iloc[k]
    rng = np.random.default_rng(SEED + 1000 * k + j)   # same draws with and without phase
    mstar = r.mstar_msun if np.isfinite(r.mstar_msun) else 1.0
    s_ccd, s_tr = rng.uniform(0, 0.05, 2)
    if r.mass_type == "true" and np.isfinite(r.msini_mj):
        sini = min(r.msini_mj / r.mass_mj, 1.0)
        inc = np.arcsin(sini) if rng.uniform() < 0.5 else np.pi - np.arcsin(sini)
        m_mj = r.mass_mj
    else:
        d_i = np.deg2rad(SIGMA_DI_DEG) * np.sqrt(-2 * np.log(rng.uniform()))
        inc = np.arccos(np.clip(np.sin(d_i) * np.cos(rng.uniform(0, 2 * np.pi)), -1, 1))
        m_mj = r.mass_mj if r.mass_type == "true" else r.mass_mj / np.sin(inc)
    period = float(r.period_yr)
    e = float(r.ecc) if np.isfinite(r.ecc) else rng.uniform(0, 0.5)
    q = m_mj * MJUP_MSUN / mstar
    a_rel = ((mstar + m_mj * MJUP_MSUN) * period**2) ** (1 / 3)
    plx = 1000.0 / r.dist_pc
    a0 = a_rel * q / (1 + q) * plx                       # photocentre semi-major axis, mas (dark companion)
    tp_rand, om_rand, big_om = rng.uniform(0, period), rng.uniform(0, 2 * np.pi), rng.uniform(0, np.pi)
    tp_shift = rng.standard_normal()
    if use_phase:   # published periastron time (BJD ~ TCB JD at this level) and stellar omega
        tp = (r.tp_bjd - JD_REF_EPOCH) / DAYS_PER_YEAR
        if np.isfinite(r.tp_err_d):
            tp += tp_shift * r.tp_err_d / DAYS_PER_YEAR
        om = np.deg2rad(r.omega_deg) if np.isfinite(r.omega_deg) else om_rand
    else:
        tp, om = tp_rand, om_rand
    orbit = dict(period=period, e=e, tp=tp, a0=a0, inc=inc, omega=om, big_omega=big_om)
    star = [0.0, 0.0, plx, rng.normal(0, 30), rng.normal(0, 30)]
    ep = simulate_epochs(_SL, r.ra, r.dec, r.phot_g_mean_mag, star, rng, orbit=orbit,
                         s_ccd=s_ccd, s_tr=s_tr)
    sig_eff2 = sigma_ccd_holl2023(r.phot_g_mean_mag) ** 2 + s_ccd**2
    orb = search_orbit_ml(ep, model="5p")
    a7 = fit_acceleration_ml(ep, order=7, model="5p")
    a9 = fit_acceleration_ml(ep, order=9, model="5p")
    return dict(k=k, draw=j, phase_known=use_phase, hostname=r.hostname, companion=r.companion, inc_deg=np.rad2deg(inc),
                mass_mj=m_mj, e=e, a0_mas=a0, n_obs=len(ep.w),
                dchi2_exp=len(ep.w) * a0**2 * (1 + np.cos(inc) ** 2) / (4 * sig_eff2),
                dlnl2_orbit=orb.dlnl2, period_fit_yr=orb.period, dlnl2_acc7=a7.dlnl2,
                dlnl2_acc9=a9.dlnl2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=200)
    ap.add_argument("--procs", type=int, default=4)
    args = ap.parse_args()
    from companions.cascade import Thresholds
    thr = Thresholds.from_csv()
    rows = load_rows()
    jobs = [(k, j, False) for k in range(len(rows)) for j in range(args.n)]
    jobs += [(k, j, True) for k in range(len(rows)) if np.isfinite(rows.tp_bjd.iloc[k])
             for j in range(args.n)]
    t0 = time.time()
    with Pool(args.procs, initializer=_init, initargs=(rows,)) as pool:
        out = pd.DataFrame(pool.map(run_one, jobs, chunksize=4))
    out["det_orbit"] = out.dlnl2_orbit > thr.orbit
    out["det_orbit_period"] = out.det_orbit & (abs(out.period_fit_yr / out.k.map(rows.period_yr) - 1) < 0.2)
    out["det_acc"] = (out.dlnl2_acc7 > thr.acc7) | (out.dlnl2_acc9 > thr.acc9)
    out["det_any"] = out.det_orbit | out.det_acc
    out.to_csv(ROOT / "results/known_companion_predictions_draws.csv", index=False, float_format="%.5g")
    g = out[~out.phase_known].groupby("k")
    gp = out[out.phase_known].groupby("k")
    summ = rows.assign(
        p_any_phase=gp.det_any.mean(), p_orbit_phase=gp.det_orbit.mean(), p_acc_phase=gp.det_acc.mean(),
        a0_mas_med=g.a0_mas.median(), alpha_uas_med=1000 * g.a0_mas.median(),
        dchi2_exp_med=g.dchi2_exp.median(), n_obs_med=g.n_obs.median(),
        p_orbit=g.det_orbit.mean(), p_orbit_period=g.det_orbit_period.mean(),
        p_acc=g.det_acc.mean(), p_any=g.det_any.mean(),
        dlnl2_orbit_med=g.dlnl2_orbit.median(), dlnl2_acc9_med=g.dlnl2_acc9.median(),
        n_draws=g.size())
    summ.to_csv(ROOT / "results/known_companion_predictions.csv", index=False, float_format="%.4g")
    cols = ["hostname", "companion", "mass_mj", "mass_type", "period_yr", "phot_g_mean_mag", "dist_pc",
            "alpha_uas_med", "dchi2_exp_med", "p_orbit", "p_orbit_period", "p_acc", "p_any", "p_any_phase"]
    print(summ[cols].sort_values("p_any", ascending=False).to_string(index=False, float_format="%.3g"))
    print(f"{len(out)} draws in {time.time() - t0:.0f} s; thresholds {thr}")


if __name__ == "__main__":
    main()
