"""Posterior sampling of an orbit with the linear parameters marginalised analytically.

Nonlinear parameters theta = (P [yr], e, M0, s_ccd, s_tr), with M0 the mean anomaly at t = 0
(J2017.5), so T_p = -M0 P / (2 pi). Flat priors: P within +-p_width of the ML period (a local
posterior around one mode; aliases are a separate question), 0 <= e < e_max, s >= 0.
With flat priors on the linear parameters (single star + Thiele-Innes), integrating them out gives
    ln p(theta) = lnL(theta, best linear) + 0.5 ln det Cov_linear(theta) + const.
For each theta sample, linear parameters are drawn from their conditional Gaussian.
"""
from __future__ import annotations

from dataclasses import dataclass

import emcee
import numpy as np

from companions.orbit import campbell_from_thiele_innes, orbit_gls

MJUP_MSUN = 9.5458e-4


def companion_mass_mj(a0_mas, parallax_mas, period_yr, mstar_msun, n_iter=60):
    """Dark-companion mass [M_J] from photocentre (= primary) semi-major axis a0.
    Solves a_star^3 / P^2 = M_c^3 / (M_star + M_c)^2 (AU, yr, M_sun)."""
    q = (np.asarray(a0_mas, float) / parallax_mas) ** 3 / np.asarray(period_yr, float) ** 2
    mc = np.cbrt(q * mstar_msun**2)
    for _ in range(n_iter):
        f = mc**3 / (mstar_msun + mc) ** 2 - q
        fp = mc**2 * (mc + 3 * mstar_msun) / (mstar_msun + mc) ** 3
        mc = np.maximum(mc - f / fp, 1e-12)
    return mc / MJUP_MSUN


@dataclass
class OrbitPosterior:
    period: np.ndarray       # yr
    e: np.ndarray
    tp: np.ndarray           # yr from J2017.5
    s_ccd: np.ndarray
    s_tr: np.ndarray
    linear: np.ndarray       # (n, n_linear): single star..., A, B, F, G
    a0: np.ndarray           # photocentre semi-major axis [mas]
    inc: np.ndarray          # [rad]
    omega: np.ndarray
    big_omega: np.ndarray
    acceptance: float


def sample_orbit(ep, ml, model="6p", n_walkers=32, n_steps=3000, burn=1000, seed=0,
                 p_width=0.3, e_max=0.95) -> OrbitPosterior:
    p_lo, p_hi = ml.period * (1 - p_width), ml.period * (1 + p_width)

    def unpack(x):
        period, e, m0, s_ccd, s_tr = x
        return period, e, -m0 * period / (2 * np.pi), s_ccd, s_tr

    def log_prob(x):
        period, e, tp, s_ccd, s_tr = unpack(x)
        if not (p_lo < period < p_hi and 0 <= e < e_max and s_ccd >= 0 and s_tr >= 0):
            return -np.inf
        fit = orbit_gls(ep, period, e, tp, s_ccd, s_tr, model)
        return fit.lnlike + 0.5 * np.linalg.slogdet(fit.cov)[1]

    rng = np.random.default_rng(seed)
    x0 = np.array([ml.period, ml.e, -2 * np.pi * ml.tp / ml.period, ml.s_ccd, ml.s_tr])
    scale = np.array([1e-3 * ml.period, 0.01, 0.01, 1e-3, 1e-3])
    p0 = x0 + scale * rng.standard_normal((n_walkers, 5))
    p0[:, 1] = np.clip(np.abs(p0[:, 1]), 0, e_max - 1e-3)
    p0[:, 3:] = np.abs(p0[:, 3:])
    sampler = emcee.EnsembleSampler(n_walkers, 5, log_prob)
    sampler.random_state = np.random.RandomState(seed).get_state()
    sampler.run_mcmc(p0, n_steps, progress=False)
    chain = sampler.get_chain(discard=burn, flat=True)

    lin = []
    for x in chain:
        period, e, tp, s_ccd, s_tr = unpack(x)
        fit = orbit_gls(ep, period, e, tp, s_ccd, s_tr, model)
        lin.append(rng.multivariate_normal(fit.params, fit.cov))
    lin = np.array(lin)
    a0, inc, om, big_om = campbell_from_thiele_innes(*lin[:, -4:].T)
    period, e, tp = chain[:, 0], chain[:, 1], -chain[:, 2] * chain[:, 0] / (2 * np.pi)
    return OrbitPosterior(period, e, tp, chain[:, 3], chain[:, 4], lin, a0, inc, om, big_om,
                          float(np.mean(sampler.acceptance_fraction)))
