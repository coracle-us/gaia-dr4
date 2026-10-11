import numpy as np
import pytest

from companions.epoch import EpochData
from companions.linear import single_star_design
from companions.orbit import (campbell_from_thiele_innes, fit_orbit, orbit_design,
                              solve_kepler, thiele_innes)
from conftest import GAIA4


@pytest.mark.parametrize("e", [0.0, 0.1, 0.5, 0.9, 0.97])
def test_solve_kepler_satisfies_keplers_equation(e):
    m = np.linspace(-20, 20, 2001)
    big_e = solve_kepler(m, e)
    np.testing.assert_allclose(big_e - e * np.sin(big_e), m, atol=1e-12)


def test_thiele_innes_roundtrip():
    rng = np.random.default_rng(1)
    for _ in range(200):
        a, inc = rng.uniform(0.01, 5), rng.uniform(0.01, np.pi - 0.01)
        om, big_om = rng.uniform(0, 2 * np.pi), rng.uniform(0, np.pi)
        ti = thiele_innes(a, inc, om, big_om)
        a2, inc2, om2, big_om2 = campbell_from_thiele_innes(*ti)
        assert a2 == pytest.approx(a, rel=1e-10)
        assert inc2 == pytest.approx(inc, abs=1e-8)
        np.testing.assert_allclose(thiele_innes(a2, inc2, om2, big_om2), ti, atol=1e-10)


def test_face_on_circular_orbit_traces_circle_of_radius_a():
    # i = 0, e = 0: photocentre offset (B X + G Y, A X + F Y) has constant length a
    a = 0.3
    ti = thiele_innes(a, 0.0, 0.4, 1.1)
    big_e = np.linspace(0, 2 * np.pi, 50)
    x, y = np.cos(big_e), np.sin(big_e)
    A, B, F, G = ti
    np.testing.assert_allclose(np.hypot(B * x + G * y, A * x + F * y), a, atol=1e-12)


def test_fit_at_true_period_recovers_injected_orbit_exactly(prerelease_table):
    ep = EpochData.from_table(prerelease_table, GAIA4)
    P, e, tp = 1.56, 0.3, 0.2                      # yr, -, yr from J2017.5
    ti_true = np.array(thiele_innes(0.4, np.deg2rad(60), 1.0, 2.0))
    star_true = np.array([0.1, -0.2, 13.6, -75.5, 17.9])
    ep.w = (single_star_design(ep, "5p") @ star_true
            + orbit_design(ep, P, e, tp) @ ti_true)
    fit = fit_orbit(ep, P, e, tp, noise=0.0, model="5p")
    np.testing.assert_allclose(fit.params[:5], star_true, atol=1e-9)
    np.testing.assert_allclose(fit.params[5:], ti_true, atol=1e-9)
    assert fit.chi2 == pytest.approx(0.0, abs=1e-12)
