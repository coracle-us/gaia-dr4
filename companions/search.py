"""Orbit search: grid in (frequency, e, T_p) with linear Thiele-Innes + single-star parameters
profiled out, then nonlinear refinement of (P, e, T_p) from the best grid minima.

Delta-chi2 is chi2(single star) - chi2(single star + orbit), same model and noise for both.
Default period range 0.1-4 yr: Lammers & Winn (2025) find DR4 periods unreliable beyond ~4 yr.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import least_squares

from companions.linear import LinearFit, fit_single_star
from companions.orbit import fit_orbit


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
                  e_grid=(0.0, 0.25, 0.5, 0.75), n_tp=8) -> GridResult:
    baseline = np.ptp(ep.t)
    freqs = np.arange(1 / p_max, 1 / p_min, 1 / (oversample * baseline))
    periods = 1 / freqs
    chi2 = np.full(len(periods), np.inf)
    best_e, best_tp = np.zeros(len(periods)), np.zeros(len(periods))
    for k, period in enumerate(periods):
        for e in e_grid:
            # circular orbits: all T_p span the same linear space
            for tp in ([0.0] if e == 0 else period * np.arange(n_tp) / n_tp):
                c2 = fit_orbit(ep, period, e, tp, noise, model).chi2
                if c2 < chi2[k]:
                    chi2[k], best_e[k], best_tp[k] = c2, e, tp
    return GridResult(periods, chi2, best_e, best_tp,
                      fit_single_star(ep, model, noise).chi2)


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
