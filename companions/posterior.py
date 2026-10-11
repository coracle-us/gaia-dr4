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


def orbit_log_prob(x, ep, model, p_lo, p_hi, e_max) -> float:
    """Marginal log posterior of theta = (P, e, M0, s_ccd, s_tr); M0 is wrapped to [-pi, pi)
    (walkers can drift by many cycles) and numerical failures count as zero probability."""
    period, e, m0, s_ccd, s_tr = x
    if not (p_lo < period < p_hi and 0 <= e < e_max and s_ccd >= 0 and s_tr >= 0
            and np.isfinite(m0)):
        return -np.inf
    m0 = np.mod(m0 + np.pi, 2 * np.pi) - np.pi
    try:
        fit = orbit_gls(ep, period, e, -m0 * period / (2 * np.pi), s_ccd, s_tr, model)
        lp = fit.lnlike + 0.5 * np.linalg.slogdet(fit.cov)[1]
    except np.linalg.LinAlgError:
        return -np.inf
    return float(lp) if np.isfinite(lp) else -np.inf


def sample_orbit(ep, ml, model="6p", n_walkers=32, n_steps=3000, burn=1000, seed=0,
                 p_width=0.3, e_max=0.95) -> OrbitPosterior:
    p_lo, p_hi = ml.period * (1 - p_width), ml.period * (1 + p_width)

    def unpack(x):
        period, e, m0, s_ccd, s_tr = x
        m0 = np.mod(m0 + np.pi, 2 * np.pi) - np.pi
        return period, e, -m0 * period / (2 * np.pi), s_ccd, s_tr

    def log_prob(x):
        return orbit_log_prob(x, ep, model, p_lo, p_hi, e_max)

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
    m0 = np.mod(chain[:, 2] + np.pi, 2 * np.pi) - np.pi
    period, e, tp = chain[:, 0], chain[:, 1], -m0 * chain[:, 0] / (2 * np.pi)
    return OrbitPosterior(period, e, tp, chain[:, 3], chain[:, 4], lin, a0, inc, om, big_om,
                          float(np.mean(sampler.acceptance_fraction)))


def thiele_innes_jacobian(a, inc, omega, big_omega, h=1e-6) -> float:
    """|d(A, B, F, G) / d(a, i, omega, Omega)| (numerical); equals a^3 sin^3 i."""
    from companions.orbit import thiele_innes
    x0 = np.array([a, inc, omega, big_omega], float)
    jac = np.zeros((4, 4))
    for k in range(4):
        dx = np.zeros(4)
        dx[k] = h
        jac[:, k] = (np.array(thiele_innes(*(x0 + dx))) - np.array(thiele_innes(*(x0 - dx)))) / (2 * h)
    return float(abs(np.linalg.det(jac)))


def campbell_log_prob(x, ep, model, p_lo, p_hi, e_max, ln_a_lo, ln_a_hi) -> float:
    """theta = (P, e, M0, ln a0, cos i, omega, Omega, s_ccd, s_tr). Priors: P uniform (local),
    e uniform, M0 / omega / Omega uniform (periodic), ln a0 uniform, cos i uniform (isotropic),
    jitters flat >= 0. The single-star parameters are marginalised analytically (flat priors)."""
    from companions.linear import single_star_design, single_star_prior
    from companions.noise import gls
    from companions.orbit import orbit_design, thiele_innes
    period, e, m0, ln_a, cos_i, om, big_om, s_ccd, s_tr = x
    if not (p_lo < period < p_hi and 0 <= e < e_max and ln_a_lo < ln_a < ln_a_hi
            and -1 <= cos_i <= 1 and s_ccd >= 0 and s_tr >= 0 and np.all(np.isfinite(x))):
        return -np.inf
    m0 = np.mod(m0 + np.pi, 2 * np.pi) - np.pi
    ti = np.array(thiele_innes(np.exp(ln_a), np.arccos(cos_i), om, big_om))
    try:
        y = ep.w - orbit_design(ep, period, e, -m0 * period / (2 * np.pi)) @ ti
        fit = gls(single_star_design(ep, model), y, ep.sigma**2 + s_ccd**2, ep.transit, s_tr**2,
                  single_star_prior(model))
        lp = fit.lnlike + 0.5 * np.linalg.slogdet(fit.cov)[1]
    except np.linalg.LinAlgError:
        return -np.inf
    return float(lp) if np.isfinite(lp) else -np.inf


def sample_orbit_campbell(ep, ml, model="6p", n_walkers=40, n_steps=3000, burn=1200, seed=0,
                          p_width=0.3, e_max=0.95) -> OrbitPosterior:
    """Posterior in Campbell elements with isotropic orientation and log-uniform a0 priors.
    Use this for masses and inclinations: flat Thiele-Innes priors (sample_orbit) imply
    p(a0, i) ~ a0^3 sin^3 i, biasing low-S/N masses high and inclinations towards edge-on."""
    from companions.linear import single_star_design, single_star_prior
    from companions.noise import gls
    from companions.orbit import orbit_design, thiele_innes
    a0, inc, om, big_om = campbell_from_thiele_innes(*ml.fit.params[-4:])
    p_lo, p_hi = ml.period * (1 - p_width), ml.period * (1 + p_width)
    ln_a_lo, ln_a_hi = np.log(a0) - 5, np.log(a0) + 3
    x0 = np.array([ml.period, min(ml.e, e_max - 1e-3), -2 * np.pi * ml.tp / ml.period, np.log(a0),
                   np.cos(inc), om, big_om, ml.s_ccd, ml.s_tr])
    scale = np.array([1e-3 * ml.period, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 1e-3, 1e-3])
    rng = np.random.default_rng(seed)
    p0 = x0 + scale * rng.standard_normal((n_walkers, len(x0)))
    p0[:, 1] = np.clip(np.abs(p0[:, 1]), 0, e_max - 1e-3)
    p0[:, 4] = np.clip(p0[:, 4], -1, 1)
    p0[:, 7:] = np.abs(p0[:, 7:])
    sampler = emcee.EnsembleSampler(n_walkers, len(x0), campbell_log_prob,
                                    args=(ep, model, p_lo, p_hi, e_max, ln_a_lo, ln_a_hi))
    sampler.random_state = np.random.RandomState(seed).get_state()
    sampler.run_mcmc(p0, n_steps, progress=False)
    chain = sampler.get_chain(discard=burn, flat=True)

    period, e = chain[:, 0], chain[:, 1]
    m0 = np.mod(chain[:, 2] + np.pi, 2 * np.pi) - np.pi
    tp = -m0 * period / (2 * np.pi)
    a0s, incs = np.exp(chain[:, 3]), np.arccos(chain[:, 4])
    om_s, big_om_s = chain[:, 5], chain[:, 6]
    shift = np.floor(big_om_s / np.pi)                     # fold (omega, Omega) ~ (omega+pi, Omega+pi)
    big_om_s = big_om_s - shift * np.pi
    om_s = np.mod(om_s - shift * np.pi, 2 * np.pi)
    lin = []
    for k in range(len(chain)):
        ti = np.array(thiele_innes(a0s[k], incs[k], om_s[k], big_om_s[k]))
        y = ep.w - orbit_design(ep, period[k], e[k], tp[k]) @ ti
        fit = gls(single_star_design(ep, model), y, ep.sigma**2 + chain[k, 7] ** 2, ep.transit,
                  chain[k, 8] ** 2, single_star_prior(model))
        lin.append(np.concatenate([rng.multivariate_normal(fit.params, fit.cov), ti]))
    return OrbitPosterior(period, e, tp, chain[:, 7], chain[:, 8], np.array(lin), a0s, incs,
                          om_s, big_om_s, float(np.mean(sampler.acceptance_fraction)))
