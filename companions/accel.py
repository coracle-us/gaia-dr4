"""Acceleration models for long-period companions (P >> 5.5 yr baseline).

Same definition as Gaia DR3 nss_acceleration_astro (Halbwachs+2023), at the DR4 reference epoch:
  7-parameter: da* += g_a* t^2 / 2,  dd += g_d t^2 / 2                          [mas/yr^2]
  9-parameter: additionally  da* += gdot_a* t^3 / 6,  dd += gdot_d t^3 / 6        [mas/yr^3]
Fitted with ML per-CCD and per-transit jitter; dlnl2 is relative to the ML single-star fit.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from companions.linear import single_star_design, single_star_prior
from companions.noise import NoiseFit, fit_noise


def acceleration_columns(ep, order: int = 7) -> np.ndarray:
    s, c, t = np.sin(ep.psi), np.cos(ep.psi), ep.t
    cols = [s * t**2 / 2, c * t**2 / 2]
    if order == 9:
        cols += [s * t**3 / 6, c * t**3 / 6]
    elif order != 7:
        raise ValueError("order must be 7 or 9")
    return np.column_stack(cols)


@dataclass
class AccelResult:
    order: int
    fit: object               # GLSFit; params = [single star..., accelerations]
    s_ccd: float
    s_tr: float
    star: NoiseFit

    @property
    def dlnl2(self) -> float:
        return 2 * (self.fit.lnlike - self.star.fit.lnlike)


def fit_acceleration_ml(ep, order: int = 7, model: str = "6p") -> AccelResult:
    star_design = single_star_design(ep, model)
    star_prior = single_star_prior(model)
    star = fit_noise(star_design, ep.w, ep.sigma, ep.transit, star_prior)
    design = np.hstack([star_design, acceleration_columns(ep, order)])
    prior = None if star_prior is None else np.concatenate([star_prior, np.full(order - 5, np.nan)])
    res = fit_noise(design, ep.w, ep.sigma, ep.transit, prior)
    return AccelResult(order, res.fit, res.s_ccd, res.s_tr, star)
