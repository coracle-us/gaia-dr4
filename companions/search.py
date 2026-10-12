"""Orbit search: grid in (frequency, e, T_p) with linear Thiele-Innes + single-star parameters
profiled out, then nonlinear refinement of (P, e, T_p) from the best grid minima.

Delta-chi2 is chi2(single star) - chi2(single star + orbit), same model and noise for both.
search_orbit_ml instead fits per-CCD and per-transit jitter by maximum likelihood for each model
and reports dlnl2 = 2 (lnL_orbit - lnL_star), which is robust to transit-correlated noise.
Default period range 0.1-4 yr: Lammers & Winn (2025) find DR4 periods unreliable beyond ~4 yr.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import least_squares, minimize

from companions.linear import LinearFit, fit_single_star, single_star_design, single_star_prior
from companions.noise import GLSFit, NoiseFit, fit_noise, gls
from companions.orbit import fit_orbit, orbit_gls


@dataclass
class GridResult:
    periods: np.ndarray       # yr
    chi2: np.ndarray          # min over (e, T_p) at each period
    e: np.ndarray             # best e at each period
    tp: np.ndarray            # best T_p at each period
    chi2_star: float

    @property
    def best_period(self) -> float:
        return float(self.periods[np.argmin(self.chi2)])

    @property
    def dchi2(self) -> np.ndarray:
        return self.chi2_star - self.chi2


@dataclass
class OrbitResult:
    period: float
    e: float
    tp: float
    fit: LinearFit
    chi2_star: float

    @property
    def dchi2(self) -> float:
        return self.chi2_star - self.fit.chi2


def period_search(ep, model="6p", noise=0.0, p_min=0.1, p_max=4.0, oversample=5,
                  e_grid=(0.0, 0.25, 0.5, 0.75), n_tp=8, s_tr=0.0) -> GridResult:
    """noise: per-CCD jitter (float or 'agis'); s_tr > 0 adds per-transit correlated jitter."""
    baseline = np.ptp(ep.t)
    freqs = np.arange(1 / p_max, 1 / p_min, 1 / (oversample * baseline))
    periods = 1 / freqs
    chi2 = np.full(len(periods), np.inf)
    best_e, best_tp = np.zeros(len(periods)), np.zeros(len(periods))
    for k, period in enumerate(periods):
        for e in e_grid:
            # circular orbits: all T_p span the same linear space
            for tp in ([0.0] if e == 0 else period * np.arange(n_tp) / n_tp):
                c2 = (fit_orbit(ep, period, e, tp, noise, model).chi2 if s_tr == 0 else
                      orbit_gls(ep, period, e, tp, float(noise), s_tr, model).chi2)
                if c2 < chi2[k]:
                    chi2[k], best_e[k], best_tp[k] = c2, e, tp
    if s_tr == 0:
        chi2_star = fit_single_star(ep, model, noise).chi2
    else:
        chi2_star = _star_noise_gls(ep, model, float(noise), s_tr).chi2
    return GridResult(periods, chi2, best_e, best_tp, chi2_star)


def _local_minima(chi2: np.ndarray, n: int) -> np.ndarray:
    idx = [i for i in range(len(chi2))
           if (i == 0 or chi2[i] <= chi2[i - 1]) and (i == len(chi2) - 1 or chi2[i] <= chi2[i + 1])]
    return np.array(sorted(idx, key=lambda i: chi2[i])[:n], dtype=int)


def refine_orbit(ep, grid: GridResult, model="6p", noise=0.0, n_starts=5,
                 e_max=0.95) -> OrbitResult:
    def whitened(x):
        f = fit_orbit(ep, x[0], x[1], x[2], noise, model)
        return f.residuals / np.sqrt(f.variance)

    best = None
    for i in _local_minima(grid.chi2, n_starts):
        p0 = grid.periods[i]
        x0 = [p0, min(grid.e[i], e_max - 1e-3), grid.tp[i]]
        lo, hi = [0.8 * p0, 0.0, x0[2] - 2 * p0], [1.25 * p0, e_max, x0[2] + 2 * p0]
        sol = least_squares(whitened, x0, bounds=(lo, hi), x_scale=[0.01 * p0, 0.1, 0.1 * p0],
                            xtol=1e-12, ftol=1e-12, gtol=1e-12)
        period, e, tp = sol.x
        fit = fit_orbit(ep, period, e, tp, noise, model)
        if best is None or fit.chi2 < best.fit.chi2:
            best = OrbitResult(float(period), float(e), float(tp), fit, grid.chi2_star)
    return best


def _star_noise_gls(ep, model, s_ccd, s_tr):
    return gls(single_star_design(ep, model), ep.w, ep.sigma**2 + s_ccd**2, ep.transit,
               s_tr**2, single_star_prior(model))


@dataclass
class MLOrbitResult:
    period: float
    e: float
    tp: float
    s_ccd: float
    s_tr: float
    fit: GLSFit
    star: NoiseFit            # ML single-star fit (its own jitters)
    grid: GridResult | None = None   # period profile at the single-star ML noise

    @property
    def dlnl2(self) -> float:
        return 2 * (self.fit.lnlike - self.star.fit.lnlike)


def search_orbit_ml(ep, model="6p", n_starts=5, e_max=0.95, **grid_kw) -> MLOrbitResult:
    """ML single-star noise fit -> period grid at that noise -> joint ML refinement of
    (P, e, T_p, s_ccd, s_tr) from the best grid minima."""
    star = fit_noise(single_star_design(ep, model), ep.w, ep.sigma, ep.transit, single_star_prior(model))
    grid = period_search(ep, model=model, noise=star.s_ccd, s_tr=star.s_tr, **grid_kw)

    def nll(x):
        return -orbit_gls(ep, x[0], x[1], x[2], x[3], x[4], model).lnlike

    best = None
    for i in _local_minima(grid.chi2, n_starts):
        p0 = grid.periods[i]
        x0 = [p0, min(grid.e[i], e_max - 1e-3), grid.tp[i], star.s_ccd, star.s_tr]
        bounds = [(0.8 * p0, 1.25 * p0), (0.0, e_max), (x0[2] - 2 * p0, x0[2] + 2 * p0),
                  (0.0, 5.0), (0.0, 5.0)]
        sol = minimize(nll, x0, method="Powell", bounds=bounds,
                       options={"xtol": 1e-8, "ftol": 1e-10, "maxfev": 20000})
        if best is None or sol.fun < best.fun:
            best = sol
    period, e, tp, s_ccd, s_tr = (float(v) for v in best.x)
    return MLOrbitResult(period, e, tp, s_ccd, s_tr,
                         orbit_gls(ep, period, e, tp, s_ccd, s_tr, model), star, grid)
