import numpy as np
import pytest
from astropy.io.votable import parse_single_table

from companions.dr4io import epochs_to_dr4_table, write_votable
from companions.epoch import EpochData
from companions.scanlaw import ScanLaw
from companions.simulate import simulate_epochs
from conftest import GAIA4


@pytest.fixture(scope="module")
def scanlaw():
    try:
        return ScanLaw.load()
    except FileNotFoundError:
        pytest.skip("gaiascanlaw data not downloaded")


def _assert_same(a: EpochData, b: EpochData):
    np.testing.assert_allclose(a.t, b.t, atol=1e-12)
    for name in ("psi", "plx_factor", "colour_factor", "w", "sigma"):
        np.testing.assert_allclose(getattr(a, name), getattr(b, name), rtol=1e-6, atol=1e-9)
    np.testing.assert_array_equal(a.transit, b.transit)
    assert a.agis_source_excess_noise == pytest.approx(b.agis_source_excess_noise, abs=1e-6)


def test_simulated_epochs_roundtrip_through_dr4_votable(scanlaw, tmp_path):
    ep = simulate_epochs(scanlaw, 150.0, 20.0, 12.0, [0, 0, 8.0, -10.0, 4.0],
                         np.random.default_rng(0), source_id=123456789012345678)
    path = tmp_path / "EPOCH_ASTROMETRY-123456789012345678.xml"
    write_votable(epochs_to_dr4_table(ep, ra0=150.0, dec0=20.0), path)
    back = EpochData.from_table(parse_single_table(path).to_table(), 123456789012345678)
    _assert_same(ep, back)


def test_real_prerelease_epochs_roundtrip(prerelease_table, tmp_path):
    ep = EpochData.from_table(prerelease_table, GAIA4)
    path = tmp_path / "gaia4.xml"
    write_votable(epochs_to_dr4_table(ep, ra0=209.5, dec0=31.7), path)
    back = EpochData.from_table(parse_single_table(path).to_table(), GAIA4)
    _assert_same(ep, back)
