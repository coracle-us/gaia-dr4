"""Per-host model-selection cascade with vetting flags.

Steps (all with ML per-CCD and per-transit jitter):
  1. single star; 2. 7p / 9p accelerations; 3. orbit search (P = 0.1-4 yr).
Classification (thresholds from the null calibration, results/thresholds.csv, FAP per star):
  orbit         dlnl2_orbit > T_orbit, and if an acceleration is also significant, the orbit must
                beat it by > 10 in dlnl2 with fitted P < 3.8 yr (calibrated on
                results/injection_recovery.csv: 93% correct when both fire, 1% of P > 4 yr
                companions called orbits; misclassified short periods all have P = 3.5-4 yr);
  acceleration  otherwise, if any statistic is significant (long-period companion);
  none          nothing significant -> upper limits (results/target_sensitivity.csv).
Orbit characterisation: ML Thiele-Innes -> a0, i; companion mass for a dark companion; minimum
mutual inclination with the transiting planet, min(|i - i_b|, |i - (180 - i_b)|) (the transit
does not fix the sense of i_b, and astrometry gives no relative node).
Optional: posterior sampling (companions.posterior) and joint fit with literature RVs.
Vetting flags (notes/phase1_validation.md):
  high_e_check (e > 0.8; null false positives are e > 0.9 fits), unstable_subsets (orbit not
  recovered when any one fifth of the transits, in time order, is dropped), marginal
  (dlnl2 < 1.25 T), twin_excluded_lit_rv / impostor_not_excluded, stellar_companion_check
  (RUWE > 1.4 or prior HGCA / PMa / DR3 NSS evidence), brown_dwarf_mass (13-80 M_J) /
  stellar_mass (> 80 M_J; dark-companion mass, so a luminous companion is heavier still).
"""
from __future__ import annotations

import dataclasses
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from companions.accel import fit_acceleration_ml
from companions.orbit import campbell_from_thiele_innes
from companions.posterior import companion_mass_mj
from companions.search import search_orbit_ml

ROOT = Path(__file__).resolve().parents[1]
ORBIT_OVER_ACCEL_DLNL2 = 10.0
ORBIT_MAX_FIT_PERIOD_YR = 3.8
N_STABILITY_BLOCKS = 5


@dataclass
class Thresholds:
    orbit: float
    acc7: float
    acc9: float

    @classmethod
    def from_csv(cls, path: Path = ROOT / "results/thresholds.csv", fap: float = 1e-4):
        t = pd.read_csv(path)
        get = lambda s: float(t[(t.statistic == s) & np.isclose(t.fap_per_star, fap)].threshold.iloc[0])
        return cls(get("dlnl2_orbit"), get("dlnl2_acc7"), get("dlnl2_acc9"))


@dataclass
class HostInfo:
    hostname: str
    mstar: float
    inner_inc_deg: float = 90.0
    lit_rv_followup: bool = False
    ruwe: float = np.nan
    prior_flags: str = ""


@dataclass
class CascadeResult:
    host: HostInfo
    classification: str
    dlnl2_orbit: float
    dlnl2_acc7: float
    dlnl2_acc9: float
    s_ccd: float
    s_tr: float
    period: float = np.nan
    e: float = np.nan
    a0: float = np.nan
    inc: float = np.nan
    parallax: float = np.nan
    mass_mj: float = np.nan
    min_mutual_inc_deg: float = np.nan
    stable: bool | None = None
    flags: list = field(default_factory=list)
    posterior: dict = field(default_factory=dict)
    joint: dict = field(default_factory=dict)

    def row(self) -> dict:
        return dict(hostname=self.host.hostname, classification=self.classification,
                    dlnl2_orbit=self.dlnl2_orbit, dlnl2_acc7=self.dlnl2_acc7,
                    dlnl2_acc9=self.dlnl2_acc9, s_ccd=self.s_ccd, s_tr=self.s_tr,
                    period_yr=self.period, e=self.e, a0_mas=self.a0,
                    inc_deg=np.rad2deg(self.inc), parallax_mas=self.parallax,
                    mass_mj=self.mass_mj, min_mutual_inc_deg=self.min_mutual_inc_deg,
                    stable=self.stable, flags=";".join(self.flags),
                    **{f"post_{k}": v for k, v in self.posterior.items()},
                    **{f"joint_{k}": v for k, v in self.joint.items()})


def subset_epochs(ep, keep: np.ndarray):
    _, transit = np.unique(ep.transit[keep], return_inverse=True)
    arrays = {f.name: getattr(ep, f.name)[keep] for f in dataclasses.fields(ep)
              if isinstance(getattr(ep, f.name), np.ndarray) and f.name != "transit"}
    return dataclasses.replace(ep, transit=transit, **arrays)


def min_mutual_inclination_deg(inc_deg, inner_inc_deg):
    return float(min(abs(inc_deg - inner_inc_deg), abs(inc_deg - (180 - inner_inc_deg))))


def stability_check(ep, period, dlnl2, model) -> bool:
    t_transit = np.array([ep.t[ep.transit == k].mean() for k in range(ep.transit.max() + 1)])
    edges = np.quantile(t_transit, np.linspace(0, 1, N_STABILITY_BLOCKS + 1))
    block = np.clip(np.searchsorted(edges, t_transit, side="right") - 1, 0, N_STABILITY_BLOCKS - 1)
    for b in range(N_STABILITY_BLOCKS):
        sub = subset_epochs(ep, block[ep.transit] != b)
        r = search_orbit_ml(sub, model=model)
        if abs(r.period / period - 1) > 0.1 or r.dlnl2 < 0.5 * dlnl2:
            return False
    return True


def run_cascade(ep, host: HostInfo, thr: Thresholds, model="6p", rv=None,
                do_stability=True, do_posterior=False, posterior_kw=None) -> CascadeResult:
    orb = search_orbit_ml(ep, model=model)
    a7 = fit_acceleration_ml(ep, order=7, model=model)
    a9 = fit_acceleration_ml(ep, order=9, model=model)
    orb_sig = orb.dlnl2 > thr.orbit
    acc_sig = a7.dlnl2 > thr.acc7 or a9.dlnl2 > thr.acc9
    beats_acc = (orb.dlnl2 - max(a7.dlnl2, a9.dlnl2) > ORBIT_OVER_ACCEL_DLNL2
                 and orb.period < ORBIT_MAX_FIT_PERIOD_YR)
    if orb_sig and (not acc_sig or beats_acc):
        cls = "orbit"
    elif orb_sig or acc_sig:
        cls = "acceleration"
    else:
        cls = "none"
    res = CascadeResult(host, cls, orb.dlnl2, a7.dlnl2, a9.dlnl2, orb.star.s_ccd, orb.star.s_tr)

    if cls != "none":
        res.flags.append("twin_excluded_lit_rv" if host.lit_rv_followup else "impostor_not_excluded")
        prior = host.prior_flags or ""
        if (np.isfinite(host.ruwe) and host.ruwe > 1.4) or any(
                k in prior for k in ("hgca", "pma", "dr3_nss")):
            res.flags.append("stellar_companion_check")
    if cls == "orbit":
        a0, inc, _, _ = campbell_from_thiele_innes(*orb.fit.params[-4:])
        plx = float(orb.fit.params[2])
        res.period, res.e, res.a0, res.inc, res.parallax = orb.period, orb.e, a0, inc, plx
        res.mass_mj = float(companion_mass_mj(a0, plx, orb.period, host.mstar))
        res.min_mutual_inc_deg = min_mutual_inclination_deg(np.rad2deg(inc), host.inner_inc_deg)
        if res.mass_mj > 80:
            res.flags.append("stellar_mass")
        elif res.mass_mj > 13:
            res.flags.append("brown_dwarf_mass")
        if orb.e > 0.8:
            res.flags.append("high_e_check")
        if orb.dlnl2 < 1.25 * thr.orbit:
            res.flags.append("marginal")
        if do_stability:
            res.stable = stability_check(ep, orb.period, orb.dlnl2, model)
            if not res.stable:
                res.flags.append("unstable_subsets")
        if do_posterior:
            from companions.posterior import sample_orbit
            post = sample_orbit(ep, orb, model=model, **(posterior_kw or {}))
            mass = companion_mass_mj(post.a0, post.linear[:, 2], post.period, host.mstar)
            dinc = np.array([min_mutual_inclination_deg(x, host.inner_inc_deg)
                             for x in np.rad2deg(post.inc)])
            for name, x in (("period_yr", post.period), ("e", post.e), ("mass_mj", mass),
                            ("inc_deg", np.rad2deg(post.inc)), ("min_mutual_inc_deg", dinc)):
                lo, md, hi = np.percentile(x, [16, 50, 84])
                res.posterior.update({f"{name}_p16": lo, f"{name}_p50": md, f"{name}_p84": hi})
        if rv is not None:
            from companions.joint import fit_joint_ml
            j = fit_joint_ml(ep, rv, model=model, start=orb)
            plx_j = float(j.astro_fit.params[2])
            res.joint = dict(period_yr=j.period, e=j.e, inc_deg=np.rad2deg(j.inc),
                             mass_mj=float(companion_mass_mj(j.a0, plx_j, j.period, host.mstar)),
                             k_ms=j.k, lnlike=j.lnlike)
    return res
