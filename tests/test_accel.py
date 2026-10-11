import numpy as np
import pytest

from companions.accel import acceleration_columns, fit_acceleration_ml
from companions.epoch import EpochData
from companions.linear import single_star_design
from companions.noise import fit_noise, gls
from companions.orbit import orbit_design, thiele_innes
from conftest import GAIA4

STAR = np.array([0.1, -0.2, 13.6, -75.5, 17.9])


def _epochs(prerelease_table, sigma=0.08):
    ep = EpochData.from_table(prerelease_table, GAIA4)
    ep.sigma = np.full_like(ep.w, sigma)
    return ep


def test_acceleration_columns_follow_dr3_definition(prerelease_table):
    # DR3 nss_acceleration_astro: da* = ... + g_a* t^2 / 2 (+ gdot_a* t^3 / 6), same for dec
    ep = _epochs(prerelease_table)
    s, c, t = np.sin(ep.psi), np.cos(ep.psi), ep.t
    np.testing.assert_allclose(acceleration_columns(ep, 7), np.column_stack([s * t**2 / 2, c * t**2 / 2]))
    np.testing.assert_allclose(acceleration_columns(ep, 9)[:, 2:],
                               np.column_stack([s * t**3 / 6, c * t**3 / 6]))


def test_noise_free_acceleration_and_jerk_recovered_exactly(prerelease_table):
    ep = _epochs(prerelease_table)
    accel = np.array([0.3, -0.2, 0.05, 0.08])
    ep.w = single_star_design(ep, "5p") @ STAR + acceleration_columns(ep, 9) @ accel
    res = fit_acceleration_ml(ep, order=9, model="5p")
    np.testing.assert_allclose(res.fit.params[5:], accel, atol=1e-8)


def test_long_period_orbit_detected_as_acceleration(prerelease_table):
    ep = _epochs(prerelease_table)
    rng = np.random.default_rng(51)
    ti = np.array(thiele_innes(0.6, np.deg2rad(70), 0.4, 1.0))
    ep.w = (single_star_design(ep, "5p") @ STAR + orbit_design(ep, 14.0, 0.2, 3.0) @ ti
            + 0.08 * rng.standard_normal(len(ep.w)))
    assert fit_acceleration_ml(ep, order=7, model="5p").dlnl2 > 50


def test_pure_noise_gives_chi2_like_acceleration_statistic(prerelease_table):
    vals = []
    for seed in range(5):
        ep = _epochs(prerelease_table)
        rng = np.random.default_rng(100 + seed)
        ep.w = single_star_design(ep, "5p") @ STAR + 0.08 * rng.standard_normal(len(ep.w))
        vals.append(fit_acceleration_ml(ep, order=9, model="5p").dlnl2)
    # 4 extra linear parameters: dlnl2 ~ chi2_4 (99.9% point 18.5)
    assert max(vals) < 18.5


def test_dlnl2_is_relative_to_ml_single_star(prerelease_table):
    ep = _epochs(prerelease_table)
    rng = np.random.default_rng(7)
    ep.w = (single_star_design(ep, "5p") @ STAR + acceleration_columns(ep, 7) @ [0.05, 0.02]
            + 0.08 * rng.standard_normal(len(ep.w)))
    res = fit_acceleration_ml(ep, order=7, model="5p")
    star = fit_noise(single_star_design(ep, "5p"), ep.w, ep.sigma, ep.transit)
    assert res.dlnl2 == pytest.approx(2 * (res.fit.lnlike - star.fit.lnlike), abs=1e-6)
