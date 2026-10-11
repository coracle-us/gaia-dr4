"""Validate the orbit fitter on Gaia-4b using ESA's DR4 pre-release epoch astrometry.

Input : data/raw/dr4_prerelease/epoch_astrometry/GAIA_DR4_PRERELEASE_EPOCH_ASTROMETRY_RAW.xml
Output: results/gaia4_validation.csv (posterior summary vs. published values)
        results/gaia4_period_profile.csv (profile likelihood in P; 5p model, single jitter,
        i.e. exogaia's setup, for the cross-check in scripts/crosscheck_exogaia.py)
Reference values: Stefansson et al. 2025, AJ 169, 107 (arXiv:2410.05654), Table 2:
  Gaia+RV adopted solution, and the Gaia DR3 orbit; host mass 0.644 M_sun (abstract).
"""
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io.votable import parse_single_table

from scipy.optimize import minimize

from companions.epoch import EpochData
from companions.orbit import orbit_gls
from companions.posterior import companion_mass_mj, sample_orbit_campbell
from companions.search import search_orbit_ml

ROOT = Path(__file__).resolve().parents[1]
XML = ROOT / "data/raw/dr4_prerelease/epoch_astrometry/GAIA_DR4_PRERELEASE_EPOCH_ASTROMETRY_RAW.xml"
OUT = ROOT / "results" / "gaia4_validation.csv"
OUT_PROFILE = ROOT / "results" / "gaia4_period_profile.csv"
GAIA4, MSTAR = 1457486023639239296, 0.644

# (median, +err, -err) from Stefansson+2025 Table 2; DR3 column symmetric errors
PUBLISHED = {
    "period_d": ((571.3, 1.4, 1.3), (564, 11, 11)),
    "e": ((0.338, 0.026, 0.023), (0.51, 0.20, 0.20)),
    "inc_deg": ((116.9, 4.2, 4.4), (115.4, 4.6, 4.6)),
    "Omega_deg": ((158.6, 5.1, 5.1), (159.4, 5.0, 5.0)),
    "omega_deg": ((180.3, 5.5, 4.8), (185, 16, 16)),
    "a0_mas": ((np.nan, np.nan, np.nan), (0.312, 0.040, 0.040)),
    "mass_mj": ((11.8, 0.7, 0.7), (np.nan, np.nan, np.nan)),
}


def period_profile(ep, periods_d, model="5p"):
    """max lnL over (e, T_p, s_ccd) at fixed P, s_tr = 0."""
    rows = []
    for pd_ in periods_d:
        period, best = pd_ / 365.25, None
        for e0 in (0.25, 0.4):
            for tp0 in np.linspace(0, period, 8, endpoint=False):
                f = lambda x: -orbit_gls(ep, period, x[0], x[1], x[2], 0.0, model).lnlike
                s = minimize(f, [e0, tp0, 0.048], method="Powell",
                             bounds=[(0, 0.95), (tp0 - 2, tp0 + 2), (0, 1)],
                             options={"xtol": 1e-9, "ftol": 1e-12})
                if best is None or s.fun < best.fun:
                    best = s
        rows.append(dict(period_d=pd_, lnlike=-best.fun, e=best.x[0], s_ccd=best.x[2]))
    df = pd.DataFrame(rows)
    df["m2dlnl"] = 2 * (df["lnlike"].max() - df["lnlike"])
    return df


def main(seed: int = 2) -> None:
    ep = EpochData.from_table(parse_single_table(XML).to_table(), GAIA4)
    ml = search_orbit_ml(ep, model="6p")
    post = sample_orbit_campbell(ep, ml, model="6p", n_walkers=40, n_steps=3000, burn=1200, seed=seed)
    draws = {
        "period_d": post.period * 365.25, "e": post.e, "inc_deg": np.rad2deg(post.inc),
        "Omega_deg": np.rad2deg(post.big_omega), "omega_deg": np.rad2deg(post.omega),
        "a0_mas": post.a0,
        "mass_mj": companion_mass_mj(post.a0, post.linear[:, 2], post.period, MSTAR),
    }
    rows = []
    for k, x in draws.items():
        lo, md, hi = np.percentile(x, [16, 50, 84])
        sd = (hi - lo) / 2
        (pv, pp, pm), (dv, dp, dm) = PUBLISHED[k]
        rows.append(dict(param=k, dr4_median=md, dr4_plus=hi - md, dr4_minus=md - lo,
                         pub_gaia_rv=pv, pub_gaia_rv_err=(pp + pm) / 2,
                         z_vs_gaia_rv=(md - pv) / np.hypot(sd, (pp + pm) / 2),
                         dr3=dv, dr3_err=(dp + dm) / 2,
                         z_vs_dr3=(md - dv) / np.hypot(sd, (dp + dm) / 2)))
    df = pd.DataFrame(rows)
    OUT.parent.mkdir(exist_ok=True)
    df.to_csv(OUT, index=False, float_format="%.5g")
    print(f"ML: dlnl2 = {ml.dlnl2:.0f}; orbit-model jitter s_ccd = {ml.s_ccd:.3f} mas, "
          f"s_tr = {ml.s_tr:.3f} mas; acceptance {post.acceptance:.2f}")
    print(df.round(3).to_string(index=False))
    prof = period_profile(ep, [565, 569, 571.3, 573, 575, 577, 579, 581, 583, 587])
    prof.to_csv(OUT_PROFILE, index=False, float_format="%.6g")
    print(prof.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
