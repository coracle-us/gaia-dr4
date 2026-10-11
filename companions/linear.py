"""Single-star along-scan models and the weighted linear least-squares solver.

The 6p model reproduces ESA gaiasupdate's DR4-like `6p_constrained_colour` solution:
  w = da* sin(psi) + dd cos(psi) + plx f_plx + (mua* sin(psi) + mud cos(psi)) t + c f_colour
with a Gaussian prior on the colour parameter. Parameter order: [da*, dd, plx, mua*, mud, (c)]
in mas, mas, mas, mas/yr, mas/yr, (pseudo-colour offset in 1/nm).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

# Archive colour_factor_al is -dw/dnu_eff in mas*um; gaiasupdate converts it to dw/dnu_eff in
# mas*nm (x -1e3) and puts a prior of 0.085 um^-1 = 0.085e-3 nm^-1 on the pseudo-colour offset
# (data model: nu_eff_error fixed at 0.085 um^-1).
COLOUR_FACTOR_TO_PER_NM = -1e3
COLOUR_PRIOR_SIGMA = 0.085e-3


def single_star_design(ep, model: str = "6p") -> np.ndarray:
    s, c = np.sin(ep.psi), np.cos(ep.psi)
    cols = [s, c, ep.plx_factor, s * ep.t, c * ep.t]
    if model == "6p":
        cols.append(COLOUR_FACTOR_TO_PER_NM * ep.colour_factor)
    elif model != "5p":
        raise ValueError(f"unknown model {model!r}")
    return np.column_stack(cols)


def single_star_prior(model: str = "6p"):
    """prior_sigma vector for single_star_design columns (NaN = no prior), or None."""
    return np.array([np.nan] * 5 + [COLOUR_PRIOR_SIGMA]) if model == "6p" else None


@dataclass
class LinearFit:
    params: np.ndarray
    cov: np.ndarray
    residuals: np.ndarray
    variance: np.ndarray
    chi2: float


def wls(design: np.ndarray, y: np.ndarray, variance: np.ndarray,
        prior_sigma: np.ndarray | None = None) -> LinearFit:
    """Weighted least squares; prior_sigma (NaN = none) adds zero-mean Gaussian priors."""
    sw = 1.0 / np.sqrt(variance)
    a = design * sw[:, None]
    normal = a.T @ a
    if prior_sigma is not None:
        ps = np.asarray(prior_sigma, dtype=float)
        has = np.isfinite(ps)
        normal[np.where(has)[0], np.where(has)[0]] += 1.0 / ps[has] ** 2
    rhs = a.T @ (y * sw)
    params = np.linalg.solve(normal, rhs)
    resid = y - design @ params
    return LinearFit(params=params, cov=np.linalg.inv(normal), residuals=resid,
                     variance=variance, chi2=float(np.sum(resid**2 / variance)))


def noise_variance(ep, noise) -> np.ndarray:
    """noise = 'agis' (AGIS source excess noise, reproduces DR4) or a jitter in mas."""
    jitter = ep.agis_source_excess_noise if noise == "agis" else float(noise)
    return ep.sigma**2 + jitter**2


def fit_single_star(ep, model: str = "6p", noise="agis") -> LinearFit:
    return wls(single_star_design(ep, model), ep.w, noise_variance(ep, noise),
               single_star_prior(model))
