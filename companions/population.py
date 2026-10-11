"""Population inference: occurrence of outer companions and their mutual-inclination width.

host_completeness: probability that a host's companion drawn from a population (mass and period
log-uniform, inclination aligned or isotropic) is detected, from the calibrated logistic
completeness in expected full-orbit Delta-chi2 (results/completeness_models.csv, written by
scripts/target_sensitivity.py from the injection-recovery runs):
    dchi2_exp = N_obs a0^2 (1 + cos^2 i) / (4 sigma_eff^2).
occurrence_posterior: fraction f of hosts with such a companion, binomial likelihood with
per-host completeness p_i: L(f) = prod_det f p_i * prod_nondet (1 - f p_i), flat prior on [0, 1].
mutual_inclination_posterior: Rayleigh width sigma of the true mutual inclination dI, with the
relative node dOmega uniform (astrometry + transit do not constrain it). The companion's sky
inclination follows cos i_c = cos i_b cos dI + sin i_b sin dI cos dOmega, with the transit's
i_b <-> 180 - i_b ambiguity marginalised and Gaussian measurement errors on i_c.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

MJUP_MSUN = 9.5458e-4


@dataclass
class OccurrencePosterior:
    grid: np.ndarray
    pdf: np.ndarray

    def _cdf(self):
        c = np.cumsum(self.pdf)
        return c / c[-1]

    def quantile(self, q):
        return float(np.interp(q, self._cdf(), self.grid))

    @property
    def median(self):
        return self.quantile(0.5)

    @property
    def sigma(self):
        w = self.pdf / self.pdf.sum()
        m = np.sum(w * self.grid)
        return float(np.sqrt(np.sum(w * (self.grid - m) ** 2)))

    @property
    def upper95(self):
        return self.quantile(0.95)


def occurrence_posterior(detected, completeness, grid=None) -> OccurrencePosterior:
    det = np.asarray(detected, bool)
    p = np.clip(np.asarray(completeness, float), 0, 1)
    f = np.linspace(0, 1, 4001) if grid is None else grid
    with np.errstate(divide="ignore"):
        lnl = (np.log(np.outer(f, p[det])).sum(axis=1)
               + np.log1p(-np.outer(f, p[~det]).clip(max=1 - 1e-12)).sum(axis=1))
    lnl = np.where(np.isfinite(lnl), lnl, -np.inf)
    pdf = np.exp(lnl - lnl.max())
    return OccurrencePosterior(f, pdf)


class CompletenessModel:
    """Logistic detection probability per period bin and geometry."""

    def __init__(self, table: pd.DataFrame):
        self.t = table.sort_values("p_lo").reset_index(drop=True)

    @classmethod
    def from_csv(cls, path):
        return cls(pd.read_csv(path))

    def prob(self, dchi2, period, aligned=True):
        sub = self.t[self.t.aligned == aligned] if (self.t.aligned == aligned).any() else self.t
        edges = np.r_[sub.p_lo.to_numpy(), sub.p_hi.to_numpy()[-1]]
        idx = np.clip(np.searchsorted(edges, period, side="right") - 1, 0, len(sub) - 1)
        k, ln_d50 = sub.k.to_numpy()[idx], sub.ln_d50.to_numpy()[idx]
        z = np.clip(k * (np.log(np.maximum(dchi2, 1e-12)) - ln_d50), -50, 50)
        p = 1 / (1 + np.exp(-z))
        return np.where(np.asarray(period) > edges[-1], 0.0, p)


def host_completeness(host: dict, model: CompletenessModel, mass_range=(1, 13), period_range=(1, 4),
                      geometry="aligned", n_mc=4000, seed=0) -> float:
    """host: n_obs, sigma_eff [mas per CCD], plx [mas], mstar [M_sun]."""
    rng = np.random.default_rng(seed)
    m = np.exp(rng.uniform(*np.log(mass_range), n_mc)) * MJUP_MSUN
    per = np.exp(rng.uniform(*np.log(period_range), n_mc))
    aligned = geometry == "aligned"
    cos_i = np.zeros(n_mc) if aligned else rng.uniform(-1, 1, n_mc)
    a0 = ((host["mstar"] + m) * per**2) ** (1 / 3) * m / (host["mstar"] + m) * host["plx"]
    dchi2 = host["n_obs"] * a0**2 * (1 + cos_i**2) / (4 * host["sigma_eff"] ** 2)
    return float(np.mean(model.prob(dchi2, per, aligned)))


@dataclass
class InclinationPosterior:
    grid: np.ndarray            # sigma [deg]
    pdf: np.ndarray
    samples: np.ndarray


def mutual_inclination_posterior(i_c_deg, i_c_err_deg, i_b_deg, sigma_grid=None, n_mc=4000,
                                 seed=0, det_weight=None) -> InclinationPosterior:
    """det_weight(i_c_deg) -> relative detection probability (optional selection correction)."""
    rng = np.random.default_rng(seed)
    grid = np.linspace(0.5, 90, 360) if sigma_grid is None else sigma_grid
    i_c, err, i_b = (np.deg2rad(np.atleast_1d(np.asarray(x, float))) for x in (i_c_deg, i_c_err_deg, i_b_deg))
    u = rng.uniform(size=n_mc)
    d_om = rng.uniform(0, 2 * np.pi, n_mc)
    flip = rng.uniform(size=(len(i_c), n_mc)) < 0.5
    lnl = np.zeros(len(grid))
    for g, sigma in enumerate(np.deg2rad(grid)):
        d_i = sigma * np.sqrt(-2 * np.log(u))
        d_i = np.where(d_i > np.pi, 2 * np.pi - np.mod(d_i, 2 * np.pi), d_i)
        ib = np.where(flip, np.pi - i_b[:, None], i_b[:, None])
        cos_ic = np.cos(ib) * np.cos(d_i) + np.sin(ib) * np.sin(d_i) * np.cos(d_om)
        ic_pred = np.arccos(np.clip(cos_ic, -1, 1))
        w = np.ones_like(ic_pred) if det_weight is None else det_weight(np.rad2deg(ic_pred))
        kern = np.exp(-0.5 * ((i_c[:, None] - ic_pred) / err[:, None]) ** 2) / err[:, None]
        like = np.sum(kern * w, axis=1) / np.sum(w, axis=1)
        lnl[g] = np.sum(np.log(np.maximum(like, 1e-300)))
    pdf = np.exp(lnl - lnl.max())
    samples = rng.choice(grid, size=20000, p=pdf / pdf.sum())
    return InclinationPosterior(grid, pdf, samples)
