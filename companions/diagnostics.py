"""Plot-ready diagnostics for one host's cascade fits (for the release-day dashboard).

Everything is reduced to one point per FoV transit (CCD values averaged with weights
1 / (sigma^2 + s_ccd^2); error sqrt(1 / sum(w) + s_tr^2)), as plain lists so the dict can be
written to JSON:
  transits      t [yr from J2017.5], phase (on the orbit fit's period), err [mas],
                orbit_resid = data - single-star part of the ML orbit fit, orbit_model = its
                Keplerian part (AL projection); accel_resid / accel_model the same for the
                more significant of the 7p / 9p acceleration fits;
  period_profile  Delta-chi2 vs period of the grid search (at the single-star ML noise);
  sky_orbit     the fitted photocentre orbit on the sky (Delta alpha*, Delta delta) [mas],
                from the Thiele-Innes constants (Delta alpha* = B X + G Y, Delta delta = A X + F Y).
"""
from __future__ import annotations

import numpy as np

from companions.accel import acceleration_columns
from companions.orbit import campbell_from_thiele_innes, orbit_design, solve_kepler


def _per_transit(values, ep, s_ccd):
    w = 1.0 / (ep.sigma**2 + s_ccd**2)
    n = ep.transit.max() + 1
    sw = np.bincount(ep.transit, weights=w, minlength=n)
    return np.bincount(ep.transit, weights=w * values, minlength=n) / sw, sw


def _clean(x):
    return [None if not np.isfinite(v) else float(v) for v in np.asarray(x, float)]


def cascade_diagnostics(ep, orb, accels, thresholds) -> dict:
    """orb: companions.search.MLOrbitResult; accels: AccelResults (7p, 9p); thresholds: dict
    with orbit / acc7 / acc9."""
    orbit_part = orbit_design(ep, orb.period, orb.e, orb.tp) @ orb.fit.params[-4:]
    orbit_resid, sw = _per_transit(orb.fit.residuals + orbit_part, ep, orb.s_ccd)
    orbit_model, _ = _per_transit(orbit_part, ep, orb.s_ccd)
    err = np.sqrt(1.0 / sw + orb.s_tr**2)
    t, _ = _per_transit(ep.t, ep, orb.s_ccd)

    acc = max(accels, key=lambda a: a.dlnl2 / thresholds[f"acc{a.order}"])
    k = acc.order - 5
    accel_part = acceleration_columns(ep, acc.order) @ acc.fit.params[-k:]
    accel_resid, _ = _per_transit(acc.fit.residuals + accel_part, ep, acc.s_ccd)
    accel_model, _ = _per_transit(accel_part, ep, acc.s_ccd)

    big_a, big_b, big_f, big_g = orb.fit.params[-4:]
    a0, inc, _, _ = campbell_from_thiele_innes(big_a, big_b, big_f, big_g)
    ecc_anom = solve_kepler(np.linspace(0, 2 * np.pi, 181), orb.e)
    x, y = np.cos(ecc_anom) - orb.e, np.sqrt(1 - orb.e**2) * np.sin(ecc_anom)
    grid = orb.grid
    return dict(
        transits=dict(t=_clean(t), phase=_clean(np.mod((t - orb.tp) / orb.period, 1.0)),
                      err=_clean(err), orbit_resid=_clean(orbit_resid),
                      orbit_model=_clean(orbit_model), accel_resid=_clean(accel_resid),
                      accel_model=_clean(accel_model)),
        orbit=dict(period_yr=orb.period, e=orb.e, tp=orb.tp, dlnl2=orb.dlnl2,
                   s_ccd=orb.s_ccd, s_tr=orb.s_tr, a0_mas=float(a0), inc_deg=float(np.rad2deg(inc))),
        accel=dict(order=acc.order, dlnl2=acc.dlnl2),
        period_profile=dict(period_yr=_clean(grid.periods), dchi2=_clean(grid.dchi2),
                            best_period_yr=grid.best_period) if grid is not None else {},
        sky_orbit=dict(ra_mas=_clean(big_b * x + big_g * y), dec_mas=_clean(big_a * x + big_f * y)),
        thresholds={k: float(v) for k, v in thresholds.items()},
    )
