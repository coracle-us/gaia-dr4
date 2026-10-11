import numpy as np
import pytest

from companions.linear import fit_single_star
from companions.orbit import fit_orbit
from companions.precision import sigma_ccd_holl2023
from companions.scanlaw import T_DR4_END, ScanLaw
from companions.simulate import simulate_epochs

STAR = np.array([0.1, -0.2, 10.0, -20.0, 5.0])


@pytest.fixture(scope="module")
def scanlaw():
    try:
        return ScanLaw.load()
    except FileNotFoundError:
        pytest.skip("gaiascanlaw data not downloaded")


def test_precision_curve_matches_holl2023_table():
    # gaiamock digitisation of Holl+2023 Fig. 3 ("EDR3 adjusted"), per CCD
    assert sigma_ccd_holl2023(12.0) == pytest.approx(0.125)
    assert sigma_ccd_holl2023(16.0) == pytest.approx(0.36)
    np.testing.assert_allclose(sigma_ccd_holl2023(np.array([10.0, 14.0])), [0.13, 0.15])


def test_simulated_sampling_structure(scanlaw):
    rng = np.random.default_rng(0)
    ep = simulate_epochs(scanlaw, 210.0, 30.0, 12.0, STAR, rng, reject=0.1)
    n_pred = len(scanlaw.transits(210.0, 30.0)[0])
    n_tr = ep.transit.max() + 1
    assert np.all(np.bincount(ep.transit) == 9)
    assert 0.75 * n_pred < n_tr <= n_pred
    assert ep.t.max() < T_DR4_END
    assert np.allclose(ep.sigma, sigma_ccd_holl2023(12.0))


def test_noise_free_single_star_recovered_exactly(scanlaw):
    ep = simulate_epochs(scanlaw, 80.0, -60.0, 11.0, STAR, np.random.default_rng(1), noise=False)
    fit = fit_single_star(ep, model="5p", noise=0.0)
    np.testing.assert_allclose(fit.params, STAR, atol=1e-9)


def test_noise_has_expected_amplitude(scanlaw):
    ep = simulate_epochs(scanlaw, 300.0, 10.0, 13.0, STAR, np.random.default_rng(2))
    fit = fit_single_star(ep, model="5p", noise=0.0)
    assert fit.chi2 / (len(ep.w) - 5) == pytest.approx(1.0, abs=0.12)


def test_injected_orbit_is_exact_at_truth(scanlaw):
    orbit = dict(period=2.0, e=0.3, tp=0.4, a0=0.2, inc=np.deg2rad(90), omega=1.0, big_omega=2.0)
    ep = simulate_epochs(scanlaw, 45.0, 20.0, 12.0, STAR, np.random.default_rng(3),
                         orbit=orbit, noise=False)
    fit = fit_orbit(ep, 2.0, 0.3, 0.4, noise=0.0, model="5p")
    assert fit.chi2 == pytest.approx(0.0, abs=1e-12)


def test_transit_correlated_noise_is_shared_within_transits(scanlaw):
    ep = simulate_epochs(scanlaw, 150.0, -5.0, 12.0, np.zeros(5), np.random.default_rng(4),
                         s_tr=0.5, sigma_override=1e-6)
    for k in range(5):
        assert np.ptp(ep.w[ep.transit == k]) < 1e-4
