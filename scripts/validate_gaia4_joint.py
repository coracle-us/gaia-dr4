"""Joint DR4-pre-release astrometry + literature RV fit of Gaia-4b.

Inputs : data/raw/dr4_prerelease/epoch_astrometry/GAIA_DR4_PRERELEASE_EPOCH_ASTROMETRY_RAW.xml
         data/raw/literature_rv/stefansson2025_gaia4_rv_t5_mrt.txt
           (AJ 169, 107, Table 5 machine-readable; 23 RVs, HPF/NEID/FIES, BJD_TDB)
Output : results/gaia4_joint_validation.csv
Reference: Stefansson+2025 Table 2 "Gaia+RV (adopted)" (DR3 Gaia orbit + same RVs).
"""
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io.votable import parse_single_table
from scipy.optimize import minimize

from companions.epoch import EpochData
from companions.joint import RVData, fit_joint_ml, rv_kepler
from companions.noise import gls
from companions.posterior import companion_mass_mj
from companions.search import search_orbit_ml

ROOT = Path(__file__).resolve().parents[1]
XML = ROOT / "data/raw/dr4_prerelease/epoch_astrometry/GAIA_DR4_PRERELEASE_EPOCH_ASTROMETRY_RAW.xml"
RV_FILE = ROOT / "data/raw/literature_rv/stefansson2025_gaia4_rv_t5_mrt.txt"
OUT = ROOT / "results" / "gaia4_joint_validation.csv"
GAIA4, MSTAR = 1457486023639239296, 0.644
PUBLISHED = {"period_d": (571.3, 1.35), "e": (0.338, 0.025), "inc_deg": (116.9, 4.3),
             "Omega_deg": (158.6, 5.1), "omega_deg": (180.3, 5.2), "K_ms": (368.0, 20.0),
             "mass_mj": (11.8, 0.7)}


def read_mrt_rv(path: Path) -> pd.DataFrame:
    rows = []
    for line in path.read_text(encoding="utf-8").splitlines():
        parts = line.split()
        if len(parts) == 4 and parts[0].replace(".", "", 1).isdigit() and parts[0].startswith("24"):
            rows.append((float(parts[0]), float(parts[1]), float(parts[2]), parts[3]))
    return pd.DataFrame(rows, columns=["bjd_tdb", "rv", "err", "inst"])


def fit_rv_only(rv: RVData):
    """Keplerian + offsets + jitters to the RVs alone (multi-start ML). Returns (x, lnL),
    x = (P [yr], e, M0, K [m/s], omega, jitters...)."""
    def nll(x):
        period, e, m0, k, om = x[:5]
        r = rv.rv - rv_kepler(rv.t, period, e, -m0 * period / (2 * np.pi), k, om)
        return -gls(rv.linear_design(), r, rv.err**2 + x[5:][rv.inst] ** 2,
                    np.arange(len(r)), 0.0).lnlike
    best = None
    bounds = [(1.0, 2.5), (0, 0.9), (-20, 20), (50, 800), (-10, 10)] + [(0, 200)] * rv.n_inst
    for p0 in np.linspace(1.3, 1.7, 9):
        for om0 in np.linspace(0, 2 * np.pi, 6, endpoint=False):
            for m00 in np.linspace(0, 2 * np.pi, 6, endpoint=False):
                s = minimize(nll, [p0, 0.3, m00, 350, om0] + [10] * rv.n_inst, method="Powell",
                             bounds=bounds, options={"maxfev": 20000})
                if best is None or s.fun < best.fun:
                    best = s
    return best.x, -best.fun


def main() -> None:
    ep = EpochData.from_table(parse_single_table(XML).to_table(), GAIA4)
    df = read_mrt_rv(RV_FILE)
    rv = RVData.from_bjd_tdb(df["bjd_tdb"], df["rv"], df["err"], list(df["inst"]))
    astro_only = search_orbit_ml(ep, model="6p")
    res = fit_joint_ml(ep, rv, model="6p", start=astro_only)
    plx = res.astro_fit.params[2]
    ours = {"period_d": res.period * 365.25, "e": res.e, "inc_deg": np.rad2deg(res.inc),
            "Omega_deg": np.rad2deg(res.big_omega), "omega_deg": np.rad2deg(res.omega),
            "K_ms": res.k, "mass_mj": float(companion_mass_mj(res.a0, plx, res.period, MSTAR))}
    out = pd.DataFrame([dict(param=k, joint_ml=v, published=PUBLISHED[k][0],
                             published_err=PUBLISHED[k][1],
                             diff_in_pub_sigma=(v - PUBLISHED[k][0]) / PUBLISHED[k][1])
                        for k, v in ours.items()])
    OUT.parent.mkdir(exist_ok=True)
    out.to_csv(OUT, index=False, float_format="%.5g")
    print(f"{len(df)} RVs ({', '.join(rv.inst_names)}); a0 = {res.a0:.3f} mas, plx = {plx:.3f} mas")
    print(f"RV jitters {np.round(res.rv_jitter, 1)} m/s, offsets {np.round(res.rv_offsets, 1)} m/s; "
          f"astro jitter s_ccd = {res.s_ccd:.3f}, s_tr = {res.s_tr:.3f} mas")
    print(f"lnL joint = {res.lnlike:.2f} (astro {res.lnlike_astro:.2f}, rv {res.lnlike_rv:.2f}); "
          f"astro-only lnL = {astro_only.fit.lnlike:.2f}")
    print(out.round(3).to_string(index=False))
    x_rv, lnl_rv_only = fit_rv_only(rv)
    print(f"RV-only ML: P = {x_rv[0] * 365.25:.1f} d, e = {x_rv[1]:.3f}, K = {x_rv[3]:.1f} m/s, "
          f"omega = {np.rad2deg(x_rv[4]) % 360:.1f} deg (published RV-only: 513 d, 0.263, 355, 181.8)")
    tension = 2 * (astro_only.fit.lnlike + lnl_rv_only - res.lnlike)
    print(f"tension 2*[lnL_astro_only + lnL_rv_only - lnL_joint] = {tension:.1f} "
          f"(shared: P, e, T_p, omega, and K <-> a0 sin i)")


if __name__ == "__main__":
    main()
