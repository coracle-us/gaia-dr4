import numpy as np
import pytest

from companions.epoch import EpochData
from companions.linear import single_star_design
from companions.orbit import campbell_from_thiele_innes, orbit_design, thiele_innes
from companions.posterior import companion_mass_mj, sample_orbit
from companions.search import search_orbit_ml
from conftest import GAIA4

MJUP_MSUN = 9.5458e-4


def test_companion_mass_recovers_jupiter_around_sun_at_10pc():
    mc = MJUP_MSUN
    a_rel = (1.0 + mc) ** (1 / 3) * 11.862 ** (2 / 3)        # AU, Kepler with P in yr
    a0_mas = a_rel * mc / (1 + mc) / 10.0 * 1e3                # photocentre = star, d = 10 pc
    assert companion_mass_mj(a0_mas, 100.0, 11.862, 1.0) == pytest.approx(1.0, rel=1e-6)


def test_companion_mass_is_vectorised():
    out = companion_mass_mj(np.array([0.1, 0.3]), 13.6, 1.56, 0.644)
    assert out.shape == (2,) and out[1] > out[0] > 0


def test_campbell_conversion_is_vectorised():
    ti = np.array([thiele_innes(0.3, 1.0, 0.5, 2.0), thiele_innes(0.2, 2.0, 4.0, 0.3)]).T
    a, inc, om, big_om = campbell_from_thiele_innes(*ti)
    np.testing.assert_allclose(a, [0.3, 0.2], rtol=1e-10)
    np.testing.assert_allclose(inc, [1.0, 2.0], atol=1e-8)


@pytest.mark.slow
def test_posterior_covers_injected_orbit(prerelease_table):
    ep = EpochData.from_table(prerelease_table, GAIA4)
    rng = np.random.default_rng(41)
    truth = dict(period=1.5, e=0.3, tp=0.25, a=0.25, inc=np.deg2rad(60))
    ti = np.array(thiele_innes(truth["a"], truth["inc"], 0.9, 1.7))
    ep.sigma = np.full_like(ep.w, 0.08)
    ep.w = (single_star_design(ep, "5p") @ np.array([0.1, -0.2, 13.6, -75.5, 17.9])
            + orbit_design(ep, truth["period"], truth["e"], truth["tp"]) @ ti
            + 0.08 * rng.standard_normal(len(ep.w))
            + 0.02 * rng.standard_normal(ep.transit.max() + 1)[ep.transit])
    ml = search_orbit_ml(ep, model="5p")
    post = sample_orbit(ep, ml, model="5p", n_walkers=24, n_steps=900, burn=300, seed=1)
    for name, true in [("period", truth["period"]), ("e", truth["e"]), ("a0", truth["a"]),
                       ("inc", truth["inc"])]:
        x = getattr(post, name)
        z = abs(np.median(x) - true) / np.std(x)
        assert z < 3, f"{name}: median {np.median(x):.4f} sd {np.std(x):.4f} truth {true:.4f}"
    assert post.linear.shape == (len(post.period), 9)
