import numpy as np
import pytest

from companions.epoch import EpochData
from companions.linear import fit_single_star, single_star_design
from companions.orbit import campbell_from_thiele_innes, orbit_design, thiele_innes
from companions.search import period_search, refine_orbit
from conftest import GAIA4

MJUP_MSUN = 9.5458e-4


def _inject(ep, period, e, tp, a_mas, inc_deg, noise_mas, seed=0):
    ti = np.array(thiele_innes(a_mas, np.deg2rad(inc_deg), 1.0, 2.0))
    star = np.array([0.1, -0.2, 13.6, -75.5, 17.9])
    rng = np.random.default_rng(seed)
    ep.sigma = np.full_like(ep.w, noise_mas)
    ep.w = (single_star_design(ep, "5p") @ star + orbit_design(ep, period, e, tp) @ ti
            + noise_mas * rng.standard_normal(len(ep.w)))
    return ti


def test_noise_free_injection_found_by_grid_and_refined_exactly(prerelease_table):
    ep = EpochData.from_table(prerelease_table, GAIA4)
    ti = _inject(ep, period=1.56, e=0.3, tp=0.2, a_mas=0.4, inc_deg=60, noise_mas=1e-6)
    grid = period_search(ep, model="5p", noise=0.0)
    assert grid.best_period == pytest.approx(1.56, rel=0.03)
    best = refine_orbit(ep, grid, model="5p", noise=0.0)
    assert best.period == pytest.approx(1.56, rel=1e-6)
    assert best.e == pytest.approx(0.3, abs=1e-5)
    np.testing.assert_allclose(best.fit.params[5:], ti, atol=1e-5)


def test_noisy_edge_on_injection_recovered(prerelease_table):
    # per-CCD noise 0.1 mas, alpha 0.15 mas edge-on: expected dchi2 ~ N a^2/(4 s^2) ~ 460
    ep = EpochData.from_table(prerelease_table, GAIA4)
    _inject(ep, period=2.3, e=0.1, tp=-0.4, a_mas=0.15, inc_deg=90, noise_mas=0.1, seed=3)
    grid = period_search(ep, model="5p", noise=0.0)
    best = refine_orbit(ep, grid, model="5p", noise=0.0)
    assert best.period == pytest.approx(2.3, rel=0.05)
    assert best.dchi2 > 200
    a, inc, _, _ = campbell_from_thiele_innes(*best.fit.params[5:])
    assert np.rad2deg(inc) == pytest.approx(90, abs=15)


def test_pure_noise_gives_small_dchi2(prerelease_table):
    ep = EpochData.from_table(prerelease_table, GAIA4)
    _inject(ep, period=1.0, e=0.0, tp=0.0, a_mas=0.0, inc_deg=0, noise_mas=0.1, seed=5)
    best = refine_orbit(ep, period_search(ep, model="5p", noise=0.0), model="5p", noise=0.0)
    assert best.dchi2 < 50


def test_dchi2_is_relative_to_single_star_fit(prerelease_table):
    ep = EpochData.from_table(prerelease_table, GAIA4)
    _inject(ep, period=1.56, e=0.3, tp=0.2, a_mas=0.2, inc_deg=40, noise_mas=0.1, seed=7)
    best = refine_orbit(ep, period_search(ep, model="5p", noise=0.0), model="5p", noise=0.0)
    chi2_star = fit_single_star(ep, model="5p", noise=0.0).chi2
    assert best.dchi2 == pytest.approx(chi2_star - best.fit.chi2)


@pytest.mark.slow
def test_gaia4_prerelease_recovers_known_planet(prerelease_table):
    """Gaia-4b (Stefansson+2025): P = 571 d, 11.8 M_J around a 0.644 M_sun star.
    Expected photocentre semi-major axis ~0.28 mas at parallax 13.62 mas."""
    ep = EpochData.from_table(prerelease_table, GAIA4)
    best = refine_orbit(ep, period_search(ep, model="6p", noise=0.0), model="6p", noise=0.0)
    assert best.period * 365.25 == pytest.approx(571, rel=0.02)
    a, _, _, _ = campbell_from_thiele_innes(*best.fit.params[6:])
    assert a == pytest.approx(0.28, rel=0.2)
    assert best.dchi2 > 50
