"""Joint fit of Gaia epoch astrometry and literature radial velocities.

Convention (Gaia DR3 AstroSpectroSB1, Halbwachs+2023; nsstools): with the photocentre (= primary,
dark companion) Thiele-Innes elements A, B, F, G built from (a0, i, omega, Omega), the primary's
RV is  v = gamma + K [cos(nu + omega) + e cos(omega)],  K = 2 pi a_star sin(i) / (P sqrt(1-e^2)),
a_star = a0 / parallax (AU). Equivalently v = dz/dt with z = C X + H Y, C = a sin(w) sin(i),
H = a cos(w) sin(i), z positive receding (tests/test_joint.py). The same omega appears in both:
checked empirically on Gaia-4 (DR3 nsstools omega = 185 +- 16 deg, RV-only 181.8 +- 6.5 deg;
Stefansson+2025). Astrometry alone is invariant under (omega, Omega) -> (omega + pi, Omega + pi);
the RVs select the branch, so Omega is reported in [0, 2 pi).

Nonlinear parameters: P, e, M0 (mean anomaly at J2017.5), a0, cos i, omega, Omega, s_ccd, s_tr,
and one RV jitter per instrument. Profiled linear parameters: single-star astrometry, RV offset
per instrument, and the RV semi-amplitudes of known inner (transiting) planets on circular
orbits with fixed (P_b, T_c,b). The RV semi-amplitude of the outer orbit uses the parallax from
the profiled astrometric solution.
Times: astrometry in TCB, RVs in BJD_TDB; TCB - TDB ~ 20 s around 2017-2024, negligible here.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import minimize

from companions.epoch import DAYS_PER_YEAR, JD_REF_EPOCH
from companions.linear import single_star_design, single_star_prior
from companions.noise import GLSFit, gls
from companions.orbit import campbell_from_thiele_innes, orbit_design, solve_kepler, thiele_innes

AU_PER_YR_MS = 4740.470463          # 1 AU/yr in m/s


def semi_amplitude_ms(a0_mas, parallax_mas, period_yr, e, inc):
    a_star_au = np.asarray(a0_mas) / parallax_mas
    return 2 * np.pi * a_star_au * np.sin(inc) / (period_yr * np.sqrt(1 - e**2)) * AU_PER_YR_MS


def rv_kepler(t, period, e, tp, k, omega):
    big_e = solve_kepler(2 * np.pi * (np.asarray(t) - tp) / period, e)
    nu = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(big_e / 2), np.sqrt(1 - e) * np.cos(big_e / 2))
    return k * (np.cos(nu + omega) + e * np.cos(omega))


@dataclass
class RVData:
    t: np.ndarray            # yr from J2017.5
    rv: np.ndarray           # m/s
    err: np.ndarray          # m/s
    inst: np.ndarray         # int instrument index 0..n-1
    inst_names: list = field(default_factory=list)
    inner: list = field(default_factory=list)   # [(P_b [yr], T_c,b [yr from J2017.5]), ...]

    @classmethod
    def from_bjd_tdb(cls, bjd, rv, err, inst_labels, inner_bjd=()):
        names = list(dict.fromkeys(inst_labels))
        inst = np.array([names.index(x) for x in inst_labels])
        to_yr = lambda jd: (np.asarray(jd, float) - JD_REF_EPOCH) / DAYS_PER_YEAR
        inner = [(p_days / DAYS_PER_YEAR, float(to_yr(tc))) for p_days, tc in inner_bjd]
        return cls(to_yr(bjd), np.asarray(rv, float), np.asarray(err, float), inst, names, inner)

    @property
    def n_inst(self) -> int:
        return int(self.inst.max()) + 1

    def linear_design(self) -> np.ndarray:
        cols = [(self.inst == k).astype(float) for k in range(self.n_inst)]
        cols += [-np.sin(2 * np.pi * (self.t - tc) / p) for p, tc in self.inner]
        return np.column_stack(cols)


@dataclass
class JointResult:
    period: float
    e: float
    tp: float
    a0: float
    inc: float
    omega: float
    big_omega: float
    s_ccd: float
    s_tr: float
    rv_jitter: np.ndarray
    k: float
    rv_offsets: np.ndarray
    inner_k: np.ndarray
    lnlike_astro: float
    lnlike_rv: float
    astro_fit: GLSFit
    theta: np.ndarray

    @property
    def lnlike(self) -> float:
        return self.lnlike_astro + self.lnlike_rv


def _evaluate(theta, ep, rv, model):
    period, e, m0, a0, cosi, om, big_om, s_ccd, s_tr = theta[:9]
    jit = np.asarray(theta[9:])
    tp = -m0 * period / (2 * np.pi)
    inc = np.arccos(np.clip(cosi, -1, 1))
    ti = np.array(thiele_innes(a0, inc, om, big_om))
    y = ep.w - orbit_design(ep, period, e, tp) @ ti
    astro = gls(single_star_design(ep, model), y, ep.sigma**2 + s_ccd**2, ep.transit,
                s_tr**2, single_star_prior(model))
    k = semi_amplitude_ms(a0, astro.params[2], period, e, inc)
    r = rv.rv - rv_kepler(rv.t, period, e, tp, k, om)
    rvfit = gls(rv.linear_design(), r, rv.err**2 + jit[rv.inst] ** 2,
                np.arange(len(r)), 0.0)
    return astro, rvfit, k, tp, inc


def joint_lnlike(theta, ep, rv, model="6p") -> float:
    astro, rvfit, *_ = _evaluate(theta, ep, rv, model)
    return astro.lnlike + rvfit.lnlike


def fit_joint_ml(ep, rv: RVData, model="6p", start=None, n_restarts=2) -> JointResult:
    """start: an MLOrbitResult (astrometry-only); found with search_orbit_ml if None."""
    if start is None:
        from companions.search import search_orbit_ml
        start = search_orbit_ml(ep, model=model)
    a0, inc, om, big_om = campbell_from_thiele_innes(*start.fit.params[-4:])
    m0 = -2 * np.pi * start.tp / start.period
    jit0 = [10.0] * rv.n_inst
    p = start.period
    bounds = ([(0.8 * p, 1.25 * p), (0, 0.95), (m0 - 10, m0 + 10), (0, 10 * a0 + 1),
               (-1, 1), (-4 * np.pi, 4 * np.pi), (-4 * np.pi, 4 * np.pi), (0, 5), (0, 5)]
              + [(0, 1000)] * rv.n_inst)
    nll = lambda x: -joint_lnlike(x, ep, rv, model)

    best = None
    for branch in (0.0, np.pi):                     # (omega, Omega) and (omega+pi, Omega+pi)
        x = np.array([p, min(start.e, 0.9), m0, a0, np.cos(inc), om + branch, big_om + branch,
                      start.s_ccd, start.s_tr] + jit0)
        for _ in range(n_restarts):                 # Powell, restarted from its own optimum
            sol = minimize(nll, x, method="Powell", bounds=bounds,
                           options={"xtol": 1e-8, "ftol": 1e-11, "maxfev": 40000})
            x = sol.x
        if best is None or sol.fun < best.fun:
            best = sol
    theta = best.x
    astro, rvfit, k, tp, inc = _evaluate(theta, ep, rv, model)
    n = rv.n_inst
    return JointResult(period=theta[0], e=theta[1], tp=tp, a0=theta[3], inc=inc,
                       omega=float(np.mod(theta[5], 2 * np.pi)),
                       big_omega=float(np.mod(theta[6], 2 * np.pi)),
                       s_ccd=theta[7], s_tr=theta[8], rv_jitter=theta[9:], k=float(k),
                       rv_offsets=rvfit.params[:n], inner_k=rvfit.params[n:],
                       lnlike_astro=astro.lnlike, lnlike_rv=rvfit.lnlike, astro_fit=astro,
                       theta=theta)
