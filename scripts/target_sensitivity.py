"""Per-host DR4 sensitivity from the injection-recovery calibration.

Inputs : results/injection_recovery.csv, results/thresholds.csv, data/targets.csv,
         gaiascanlaw grid (companions.scanlaw) for each host's predicted DR4 transits
Output : results/target_sensitivity.csv; summary printed.
Method : for each period bin and geometry, fit a logistic detection probability
         p(detect) = 1 / (1 + exp(-k (ln dchi2_exp - ln d50)))
to the injections (detection: orbit dlnl2 > threshold at FAP 1e-4 per star). For each host,
dchi2_exp = N_obs a0^2 (1 + cos^2 i) / (4 sigma_eff^2), N_obs = 0.9 * 9 * N_transits(predicted),
sigma_eff^2 = sigma_ccd(G)^2 + 0.025^2 (mean injected excess). Inverting for a0 at p = 0.5, 0.9 and
converting to companion mass (dark companion, M_star from the archive) gives m50 / m90 [M_J]
for aligned (i = 90 deg) companions at P = 1, 2, 3, 4 yr.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize

from companions.posterior import MJUP_MSUN
from companions.precision import sigma_ccd_holl2023
from companions.scanlaw import ScanLaw
from companions.simulate import N_CCD

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
P_EVAL = (1.0, 2.0, 3.0, 4.0)
P_BINS = [(0.5, 1.5), (1.5, 2.5), (2.5, 3.5), (3.5, 4.5)]
S_EXCESS = 0.025


def fit_logistic(x_ln, det):
    def nll(p):
        k, ln_d50 = p
        z = np.clip(k * (x_ln - ln_d50), -50, 50)
        pr = np.clip(1 / (1 + np.exp(-z)), 1e-9, 1 - 1e-9)
        return -np.sum(det * np.log(pr) + (1 - det) * np.log(1 - pr))
    return minimize(nll, [3.0, np.log(50)], method="Nelder-Mead").x


def a0_for_mass(m_mj, period, mstar, plx):
    mc = m_mj * MJUP_MSUN
    a_rel = ((mstar + mc) * period**2) ** (1 / 3)
    return a_rel * mc / (mstar + mc) * plx


def mass_for_a0(a0, period, mstar, plx):
    lo, hi = 1e-3, 500.0
    for _ in range(80):
        mid = np.sqrt(lo * hi)
        lo, hi = (mid, hi) if a0_for_mass(mid, period, mstar, plx) < a0 else (lo, mid)
    return np.sqrt(lo * hi)


def main() -> None:
    thr = pd.read_csv(RES / "thresholds.csv")
    t_orb = float(thr[(thr.statistic == "dlnl2_orbit") & (thr.fap_per_star == 1e-4)].threshold.iloc[0])
    inj = pd.read_csv(RES / "injection_recovery.csv")
    inj["det"] = (inj.dlnl2_orbit > t_orb).astype(float)
    models = {}
    for plo, phi in P_BINS:
        s = inj[(inj.aligned) & (inj.period_true >= plo) & (inj.period_true < phi)]
        k, ln_d50 = fit_logistic(np.log(s.dchi2_exp_actual.to_numpy()), s.det.to_numpy())
        models[(plo, phi)] = (k, ln_d50)
        print(f"aligned, P {plo}-{phi} yr: n = {len(s)}, d50 = {np.exp(ln_d50):.0f}, "
              f"d90 = {np.exp(ln_d50 + np.log(9) / k):.0f}, slope k = {k:.2f}")

    sl = ScanLaw.load()
    hosts = (pd.read_csv(ROOT / "data" / "targets.csv")
             .dropna(subset=["ra", "dec", "phot_g_mean_mag", "dist_pc"])
             .drop_duplicates("hostname").reset_index(drop=True))
    rows = []
    for _, h in hosts.iterrows():
        n_obs = 0.9 * N_CCD * len(sl.transits(h.ra, h.dec)[0])
        sig2 = sigma_ccd_holl2023(h.phot_g_mean_mag) ** 2 + S_EXCESS**2
        plx, mstar = 1000.0 / h.dist_pc, (h.st_mass if np.isfinite(h.st_mass) else 1.0)
        out = dict(hostname=h.hostname, source_id=h.source_id, g_mag=h.phot_g_mean_mag,
                   dist_pc=h.dist_pc, mstar=mstar, n_obs_pred=n_obs,
                   epoch_rv_expected=h.epoch_rv_expected, ruwe=h.ruwe, pl_orbper=h.pl_orbper)
        for p_eval, key in zip(P_EVAL, P_BINS):
            k, ln_d50 = models[key]
            for prob, tag in ((0.5, "m50"), (0.9, "m90")):
                d_req = np.exp(ln_d50 + np.log(prob / (1 - prob)) / k)
                a0_req = np.sqrt(4 * sig2 * d_req / n_obs)        # edge-on: 1 + cos^2 i = 1
                out[f"{tag}_mj_P{p_eval:g}"] = mass_for_a0(a0_req, p_eval, mstar, plx)
        rows.append(out)
    df = pd.DataFrame(rows)
    df.to_csv(RES / "target_sensitivity.csv", index=False, float_format="%.4g")
    print(f"\n{len(df)} hosts; aligned companions, detection at FAP 1e-4/star (dlnl2 > {t_orb:.1f})")
    for p in P_EVAL:
        m = df[f"m50_mj_P{p:g}"]
        m9 = df[f"m90_mj_P{p:g}"]
        print(f"  P = {p:g} yr: hosts with m50 <= 1/2/5/13 MJ: "
              f"{(m <= 1).sum()}/{(m <= 2).sum()}/{(m <= 5).sum()}/{(m <= 13).sum()}   "
              f"m90 <= 5/13 MJ: {(m9 <= 5).sum()}/{(m9 <= 13).sum()}")


if __name__ == "__main__":
    main()
