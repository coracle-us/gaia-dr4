"""Keplerian photocentre orbit in the along-scan formulation (Thiele-Innes, linear at fixed P, e, T_p).

Convention (as Gaia DR3 nss_two_body_orbit and gaiamock): photocentre offsets
    da* = B X + G Y,   dd = A X + F Y,   X = cos E - e,   Y = sqrt(1 - e^2) sin E,
with mean anomaly M = 2 pi (t - T_p) / P, t and T_p in yr from J2017.5, A..G in mas.
AL displacement: dw = da* sin(psi) + dd cos(psi).
NB: the sign of the photocentre-vs-relative orbit (omega vs omega + pi) matters for joint RV fits;
check it against the DR4 documentation before combining with RVs. It does not affect detection,
the semi-major axis, or the inclination.
"""
from __future__ import annotations

import numpy as np

from companions.linear import COLOUR_PRIOR_SIGMA, LinearFit, noise_variance, single_star_design, wls
from companions.noise import GLSFit, gls


def solve_kepler(m, e: float, tol: float = 1e-14, max_iter: int = 50) -> np.ndarray:
    """Eccentric anomaly E with E - e sin E = M (Newton with Danby starting value)."""
    m = np.asarray(m, dtype=float)
    mw = np.mod(m + np.pi, 2 * np.pi) - np.pi
    big_e = mw + 0.85 * e * np.sign(np.sin(mw))
    for _ in range(max_iter):
        f = big_e - e * np.sin(big_e) - mw
        step = f / (1 - e * np.cos(big_e))
        big_e -= step
        if np.max(np.abs(step)) < tol:
            break
    return big_e + (m - mw)


def thiele_innes(a, inc, omega, big_omega):
    """(A, B, F, G) from semi-major axis a, inclination, argument of periastron, node [rad]."""
    co, so, cO, sO, ci = np.cos(omega), np.sin(omega), np.cos(big_omega), np.sin(big_omega), np.cos(inc)
    return (a * (co * cO - so * sO * ci), a * (co * sO + so * cO * ci),
            -a * (so * cO + co * sO * ci), -a * (so * sO - co * cO * ci))


def campbell_from_thiele_innes(A, B, F, G):
    """Inverse of thiele_innes: (a, inc, omega, big_omega), big_omega in [0, pi). Vectorised."""
    A, B, F, G = (np.asarray(v, dtype=float) for v in (A, B, F, G))
    k = (A**2 + B**2 + F**2 + G**2) / 2
    m = A * G - B * F
    a = np.sqrt(k + np.sqrt(np.maximum(k**2 - m**2, 0.0)))
    inc = np.arccos(np.clip(m / a**2, -1, 1))
    w_plus_o = np.arctan2(B - F, A + G)
    w_minus_o = np.arctan2(-B - F, A - G)
    omega, big_omega = (w_plus_o + w_minus_o) / 2, (w_plus_o - w_minus_o) / 2
    shift = np.where(big_omega < 0, np.pi, np.where(big_omega >= np.pi, -np.pi, 0.0))
    omega, big_omega = omega + shift, big_omega + shift
    out = a, inc, np.mod(omega, 2 * np.pi), big_omega
    return tuple(float(v) for v in out) if np.ndim(a) == 0 else out


def orbit_design(ep, period: float, e: float, tp: float) -> np.ndarray:
    """Columns multiplying (A, B, F, G) in the AL equation."""
    big_e = solve_kepler(2 * np.pi * (ep.t - tp) / period, e)
    x, y = np.cos(big_e) - e, np.sqrt(1 - e**2) * np.sin(big_e)
    s, c = np.sin(ep.psi), np.cos(ep.psi)
    return np.column_stack([c * x, s * x, c * y, s * y])


def star_orbit_design(ep, period: float, e: float, tp: float, model: str = "6p"):
    """Design matrix [single star..., A, B, F, G] and matching prior_sigma vector."""
    design = np.hstack([single_star_design(ep, model), orbit_design(ep, period, e, tp)])
    prior = np.full(design.shape[1], np.nan)
    if model == "6p":
        prior[5] = COLOUR_PRIOR_SIGMA
    return design, prior


def fit_orbit(ep, period: float, e: float, tp: float, noise=0.0, model: str = "6p") -> LinearFit:
    """Single-star model + Keplerian orbit at fixed (P, e, T_p). Params: [star..., A, B, F, G]."""
    design, prior = star_orbit_design(ep, period, e, tp, model)
    return wls(design, ep.w, noise_variance(ep, noise), prior)


def orbit_gls(ep, period: float, e: float, tp: float, s_ccd: float = 0.0, s_tr: float = 0.0,
              model: str = "6p") -> GLSFit:
    """As fit_orbit, with per-CCD jitter s_ccd and per-transit correlated jitter s_tr [mas]."""
    design, prior = star_orbit_design(ep, period, e, tp, model)
    return gls(design, ep.w, ep.sigma**2 + s_ccd**2, ep.transit, s_tr**2, prior)
