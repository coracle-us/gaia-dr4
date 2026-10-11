import numpy as np
import pytest

from companions.epoch import EpochData
from companions.linear import single_star_design
from companions.orbit import campbell_from_thiele_innes, orbit_design, thiele_innes
from companions.search import period_search, refine_orbit, search_orbit_ml
from conftest import GAIA4


def _simulate(prerelease_table, a_mas, s_ccd, s_tr, sigma=0.08, period=2.0, e=0.2, tp=0.3,
              inc_deg=90, seed=0):
    ep = EpochData.from_table(prerelease_table, GAIA4)
    rng = np.random.default_rng(seed)
    ep.sigma = np.full_like(ep.w, sigma)
    ti = np.array(thiele_innes(a_mas, np.deg2rad(inc_deg), 0.7, 1.3))
    ep.w = (single_star_design(ep, "5p") @ np.array([0.1, -0.2, 13.6, -75.5, 17.9])
            + orbit_design(ep, period, e, tp) @ ti
            + np.hypot(sigma, s_ccd) * rng.standard_normal(len(ep.w))
            + s_tr * rng.standard_normal(ep.transit.max() + 1)[ep.transit])
    return ep


def test_ml_search_recovers_orbit_and_noise(prerelease_table):
    ep = _simulate(prerelease_table, a_mas=0.12, s_ccd=0.03, s_tr=0.04, seed=21)
    res = search_orbit_ml(ep, model="5p")
    assert res.period == pytest.approx(2.0, rel=0.05)
    assert res.s_tr == pytest.approx(0.04, abs=0.02)
    assert res.dlnl2 > 50
    a, inc, _, _ = campbell_from_thiele_innes(*res.fit.params[5:])
    assert a == pytest.approx(0.12, rel=0.3)


def test_transit_correlated_noise_alone_is_not_a_detection(prerelease_table):
    """Per-transit noise fools a fixed-noise per-CCD chi2 test but not the ML test."""
    ep = _simulate(prerelease_table, a_mas=0.0, s_ccd=0.0, s_tr=0.06, seed=22)
    fixed = refine_orbit(ep, period_search(ep, model="5p", noise=0.0), model="5p", noise=0.0)
    ml = search_orbit_ml(ep, model="5p")
    assert fixed.dchi2 > 50          # naive test would claim a detection
    assert ml.dlnl2 < 50


def test_dlnl2_compares_ml_fits_of_both_models(prerelease_table):
    ep = _simulate(prerelease_table, a_mas=0.1, s_ccd=0.02, s_tr=0.02, seed=23)
    res = search_orbit_ml(ep, model="5p")
    assert res.dlnl2 == pytest.approx(2 * (res.fit.lnlike - res.star.fit.lnlike))
    assert res.fit.lnlike >= res.star.fit.lnlike
