import numpy as np
import pytest
from astropy.table import Table

from companions.gaia_rv import rv_constancy, rvdata_from_gaia_epoch

SID = 42


def _table(t_days, rv_kms, err_kms, valid=None):
    n = len(t_days)
    return Table({"source_id": np.full(n, SID, np.int64), "obs_time_rv": np.asarray(t_days, float),
                  "radial_velocity": np.asarray(rv_kms, float),
                  "radial_velocity_error": np.asarray(err_kms, float),
                  "rv_validity_flag": np.ones(n, bool) if valid is None else np.asarray(valid)})


def test_time_and_units_conversion():
    # obs_time_rv = BJD(TCB) - 2455197.5; J2017.5 = 2457936.875 -> 2739.375 d
    rv = rvdata_from_gaia_epoch(_table([2739.375, 2739.375 + 365.25], [1.0, 2.0], [0.5, 0.5]), SID)
    np.testing.assert_allclose(rv.t, [0.0, 1.0])
    np.testing.assert_allclose(rv.rv, [1000.0, 2000.0])          # m/s
    np.testing.assert_allclose(rv.err, [500.0, 500.0])


def test_invalid_epochs_dropped():
    rv = rvdata_from_gaia_epoch(_table([2000, 2100, 2200], [1, np.nan, 3], [0.5, 0.5, 0.5],
                                       valid=[True, True, False]), SID)
    assert len(rv.t) == 1


def test_constant_star_passes_and_companion_fails():
    rng = np.random.default_rng(0)
    t = np.sort(rng.uniform(1800, 3800, 40))
    flat = rvdata_from_gaia_epoch(_table(t, 5 + 0.8 * rng.standard_normal(40), np.full(40, 0.8)), SID)
    assert rv_constancy(flat).pvalue > 1e-3
    wobble = 5 + 4.0 * np.sin(2 * np.pi * t / 700) + 0.8 * rng.standard_normal(40)
    bin_ = rvdata_from_gaia_epoch(_table(t, wobble, np.full(40, 0.8)), SID)
    assert rv_constancy(bin_).pvalue < 1e-6


def test_known_inner_planet_is_subtracted():
    rng = np.random.default_rng(1)
    t = np.sort(rng.uniform(1800, 3800, 40))
    p_b, tc_b, k_b = 3.0, 2500.0, 2.5                             # days, days, km/s
    rv_kms = 5 - k_b * np.sin(2 * np.pi * (t - tc_b) / p_b) + 0.5 * rng.standard_normal(40)
    rv = rvdata_from_gaia_epoch(_table(t, rv_kms, np.full(40, 0.5)), SID)
    assert rv_constancy(rv).pvalue < 1e-6
    inner = [(p_b / 365.25, (tc_b + 2455197.5 - 2457936.875) / 365.25, k_b * 1000)]
    assert rv_constancy(rv, inner=inner).pvalue > 1e-3
