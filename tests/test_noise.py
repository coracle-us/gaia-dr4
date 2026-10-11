import numpy as np
import pytest

from companions.linear import wls
from companions.noise import gls


def _problem(seed=0, n_transit=12, p=4):
    rng = np.random.default_rng(seed)
    transit = np.repeat(np.arange(n_transit), rng.integers(5, 10, n_transit))
    n = len(transit)
    design = rng.standard_normal((n, p))
    var_ccd = rng.uniform(0.01, 0.04, n)
    y = design @ rng.standard_normal(p) + 0.1 * rng.standard_normal(n)
    return design, y, var_ccd, transit


def _dense_cov(var_ccd, transit, s_tr2):
    same = transit[:, None] == transit[None, :]
    return np.diag(var_ccd) + s_tr2 * same


def test_gls_matches_dense_matrix_algebra():
    design, y, var_ccd, transit = _problem()
    s_tr2 = 0.03
    cov = _dense_cov(var_ccd, transit, s_tr2)
    cinv = np.linalg.inv(cov)
    normal = design.T @ cinv @ design
    params = np.linalg.solve(normal, design.T @ cinv @ y)
    r = y - design @ params
    fit = gls(design, y, var_ccd, transit, s_tr2)
    np.testing.assert_allclose(fit.params, params, rtol=1e-10)
    np.testing.assert_allclose(fit.cov, np.linalg.inv(normal), rtol=1e-10)
    assert fit.chi2 == pytest.approx(r @ cinv @ r, rel=1e-10)
    assert fit.logdet == pytest.approx(np.linalg.slogdet(cov)[1], rel=1e-10)


def test_gls_without_transit_noise_equals_wls():
    design, y, var_ccd, transit = _problem(seed=1)
    prior = np.array([np.nan, np.nan, np.nan, 0.5])
    a = gls(design, y, var_ccd, transit, 0.0, prior)
    b = wls(design, y, var_ccd, prior)
    np.testing.assert_allclose(a.params, b.params, rtol=1e-12)
    assert a.chi2 == pytest.approx(b.chi2, rel=1e-12)
    assert a.logdet == pytest.approx(np.sum(np.log(var_ccd)), rel=1e-12)


def test_log_likelihood_is_gaussian_log_density():
    from scipy.stats import multivariate_normal
    design, y, var_ccd, transit = _problem(seed=2)
    fit = gls(design, y, var_ccd, transit, 0.02)
    cov = _dense_cov(var_ccd, transit, 0.02)
    expected = multivariate_normal(mean=design @ fit.params, cov=cov).logpdf(y)
    assert fit.lnlike == pytest.approx(expected, rel=1e-10)


def _simulate_single_star(prerelease_table, s_ccd, s_tr, sigma=0.08, seed=0):
    from companions.epoch import EpochData
    from companions.linear import single_star_design
    from conftest import GAIA4
    ep = EpochData.from_table(prerelease_table, GAIA4)
    rng = np.random.default_rng(seed)
    ep.sigma = np.full_like(ep.w, sigma)
    n_tr = ep.transit.max() + 1
    ep.w = (single_star_design(ep, "5p") @ np.array([0.1, -0.2, 13.6, -75.5, 17.9])
            + np.hypot(sigma, s_ccd) * rng.standard_normal(len(ep.w))
            + s_tr * rng.standard_normal(n_tr)[ep.transit])
    return ep


def test_ml_noise_fit_recovers_injected_jitters(prerelease_table):
    from companions.linear import single_star_design
    from companions.noise import fit_noise
    ep = _simulate_single_star(prerelease_table, s_ccd=0.05, s_tr=0.06, seed=11)
    res = fit_noise(single_star_design(ep, "5p"), ep.w, ep.sigma, ep.transit)
    assert res.s_ccd == pytest.approx(0.05, abs=0.015)
    assert res.s_tr == pytest.approx(0.06, abs=0.02)


def test_ml_noise_fit_finds_no_excess_when_none_injected(prerelease_table):
    from companions.linear import single_star_design
    from companions.noise import fit_noise
    ep = _simulate_single_star(prerelease_table, s_ccd=0.0, s_tr=0.0, seed=12)
    res = fit_noise(single_star_design(ep, "5p"), ep.w, ep.sigma, ep.transit)
    assert res.s_ccd < 0.02 and res.s_tr < 0.02
    # the ML fit can never be worse than the zero-jitter fit
    from companions.noise import gls
    assert res.fit.lnlike >= gls(single_star_design(ep, "5p"), ep.w, ep.sigma**2,
                                 ep.transit, 0.0).lnlike - 1e-9
