"""Generalised least squares with per-CCD white noise plus per-transit correlated noise.

Covariance of the CCD measurements in FoV transit k:
    C_k = diag(var_ccd) + s_tr^2 * 1 1^T
(the per-transit term models noise shared by all CCDs of one transit, e.g. attitude errors).
C^-1 and log det C follow per block from the Sherman-Morrison formula, so cost is O(n p^2).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class GLSFit:
    params: np.ndarray
    cov: np.ndarray           # (A^T C^-1 A + priors)^-1
    residuals: np.ndarray
    chi2: float               # r^T C^-1 r (prior penalty not included)
    logdet: float             # log det C
    n: int

    @property
    def lnlike(self) -> float:
        return -0.5 * (self.chi2 + self.logdet + self.n * np.log(2 * np.pi))


def _block_sums(values: np.ndarray, transit: np.ndarray, n_transit: int) -> np.ndarray:
    out = np.zeros((n_transit,) + values.shape[1:])
    np.add.at(out, transit, values)
    return out


def gls(design: np.ndarray, y: np.ndarray, var_ccd: np.ndarray, transit: np.ndarray,
        s_tr2: float, prior_sigma: np.ndarray | None = None) -> GLSFit:
    n_tr = int(transit.max()) + 1
    dinv = 1.0 / var_ccd
    a_d = design * dinv[:, None]
    s_k = _block_sums(dinv, transit, n_tr)
    c_k = s_tr2 / (1.0 + s_tr2 * s_k)
    u = _block_sums(a_d, transit, n_tr)                  # sum_i D^-1 A_i per transit
    v = _block_sums(dinv * y, transit, n_tr)
    normal = design.T @ a_d - (u * c_k[:, None]).T @ u
    rhs = a_d.T @ y - u.T @ (c_k * v)
    if prior_sigma is not None:
        ps = np.asarray(prior_sigma, dtype=float)
        idx = np.where(np.isfinite(ps))[0]
        normal[idx, idx] += 1.0 / ps[idx] ** 2
    params = np.linalg.solve(normal, rhs)
    r = y - design @ params
    rk = _block_sums(dinv * r, transit, n_tr)
    chi2 = float(np.sum(dinv * r**2) - np.sum(c_k * rk**2))
    logdet = float(np.sum(np.log(var_ccd)) + np.sum(np.log1p(s_tr2 * s_k)))
    return GLSFit(params, np.linalg.inv(normal), r, chi2, logdet, len(y))


@dataclass
class NoiseFit:
    s_ccd: float              # per-CCD jitter [mas], added in quadrature to sigma
    s_tr: float               # per-transit correlated jitter [mas]
    fit: GLSFit


def fit_noise(design, y, sigma, transit, prior_sigma=None, s_max=5.0) -> NoiseFit:
    """Maximum-likelihood (s_ccd, s_tr) with the linear parameters profiled out."""
    from scipy.optimize import minimize

    def nll(s):
        return -gls(design, y, sigma**2 + s[0] ** 2, transit, s[1] ** 2, prior_sigma).lnlike

    best = None
    for s0 in ([0.0, 0.0], [0.05, 0.05], [0.3, 0.3]):
        sol = minimize(nll, s0, method="L-BFGS-B", bounds=[(0, s_max), (0, s_max)])
        if best is None or sol.fun < best.fun:
            best = sol
    s_ccd, s_tr = (float(v) for v in best.x)
    return NoiseFit(s_ccd, s_tr, gls(design, y, sigma**2 + s_ccd**2, transit, s_tr**2,
                                     prior_sigma))
