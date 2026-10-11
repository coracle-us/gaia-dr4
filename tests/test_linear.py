import numpy as np
import pytest
from gaiasupdate.epoch_astrometry import GaiaEpochAstrometryArchive

from companions.epoch import EpochData
from companions.linear import fit_single_star
from conftest import GAIA4

PRERELEASE_IDS = [2237987199365376, 2309425390592896, 10973744521070720, 20694084440761600,
                  60730287810150016, 435469040545191680, 1457486023639239296,
                  1663617687609809280, 3926186255616949504, 3937211745905473024,
                  4181040337841125632, 4318465066420528000]


@pytest.mark.parametrize("source_id", PRERELEASE_IDS)
def test_6p_fit_reproduces_gaiasupdate_dr4_like(prerelease_table, gaiasupdate_fit, source_id):
    ref = gaiasupdate_fit(source_id)
    fit = fit_single_star(EpochData.from_table(prerelease_table, source_id), model="6p",
                          noise="agis")
    np.testing.assert_allclose(fit.params, ref["parameters"], rtol=1e-9, atol=1e-9)
    np.testing.assert_allclose(fit.cov, ref["parameter_covariance_matrix_formal"],
                               rtol=1e-7, atol=1e-15)


def test_5p_fit_reproduces_gaiasupdate(prerelease_table):
    df = prerelease_table[prerelease_table["source_id"] == GAIA4].to_pandas()
    ref = GaiaEpochAstrometryArchive.supdate(df, GAIA4, model="5p_single_source")
    fit = fit_single_star(EpochData.from_table(prerelease_table, GAIA4), model="5p",
                          noise="agis")
    np.testing.assert_allclose(fit.params, ref["parameters"], rtol=1e-9, atol=1e-9)


def test_chi2_and_residuals_are_consistent(prerelease_table):
    ep = EpochData.from_table(prerelease_table, GAIA4)
    fit = fit_single_star(ep, model="5p", noise=0.0)
    assert fit.residuals.shape == ep.w.shape
    assert fit.chi2 == pytest.approx(np.sum(fit.residuals**2 / ep.sigma**2))


def test_noise_option_sets_variance(prerelease_table):
    ep = EpochData.from_table(prerelease_table, GAIA4)
    np.testing.assert_allclose(fit_single_star(ep, noise=0.05).variance, ep.sigma**2 + 0.05**2)
    np.testing.assert_allclose(fit_single_star(ep, noise="agis").variance,
                               ep.sigma**2 + ep.agis_source_excess_noise**2)
