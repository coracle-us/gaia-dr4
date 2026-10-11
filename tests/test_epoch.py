import numpy as np

from companions.epoch import EpochData
from conftest import GAIA4


def test_uses_only_ccds_used_by_agis(prerelease_table, gaiasupdate_fit):
    ep = EpochData.from_table(prerelease_table, GAIA4)
    assert len(ep.w) == gaiasupdate_fit(GAIA4)["n_measurements"] == 824


def test_relative_time_matches_gaiasupdate(prerelease_table, gaiasupdate_fit):
    ep = EpochData.from_table(prerelease_table, GAIA4)
    ref = gaiasupdate_fit(GAIA4)["timestamps_of_used_observations"]["relative_time_year"]
    np.testing.assert_allclose(np.sort(ep.t), np.sort(ref.to_numpy()), atol=1e-12)


def test_epochs_span_dr4_window_around_j2017p5(prerelease_table):
    ep = EpochData.from_table(prerelease_table, GAIA4)
    # DR4: 2014-07-25 .. 2020-01-20, i.e. about -2.93 .. +2.55 yr from J2017.5
    assert -2.95 < ep.t.min() < ep.t.max() < 2.56


def test_transit_index_groups_ccds_of_one_fov_transit(prerelease_table):
    ep = EpochData.from_table(prerelease_table, GAIA4)
    assert ep.transit.max() + 1 == len(np.unique(ep.transit))
    for k in np.unique(ep.transit)[:10]:
        assert np.ptp(ep.t[ep.transit == k]) * 365.25 * 86400 < 120  # CCDs within ~1 min
