import numpy as np
import pytest

from companions.epoch import EpochData
from companions.scanlaw import ScanLaw, parallax_factor_al
from conftest import PRERELEASE_XML

PRERELEASE_IDS = [2309425390592896, 20694084440761600, 435469040545191680,
                  1457486023639239296, 1663617687609809280, 3937211745905473024,
                  4181040337841125632, 4318465066420528000, 3926186255616949504]


@pytest.fixture(scope="module")
def scanlaw():
    try:
        return ScanLaw.load()
    except FileNotFoundError:
        pytest.skip("gaiascanlaw data not downloaded (see companions/scanlaw.py)")


def _real_transits(prerelease_table, sid):
    """Per-FoV-transit mean time [yr from J2017.5], scan angle, parallax factor, ra0, dec0."""
    ep = EpochData.from_table(prerelease_table, sid)
    rows = prerelease_table[prerelease_table["source_id"] == sid]
    out = []
    for k in np.unique(ep.transit):
        m = ep.transit == k
        out.append((ep.t[m].mean(), np.mean(ep.psi[m]), ep.plx_factor[m][0]))
    return np.array(out), float(rows["ra0"][0]), float(rows["dec0"][0])


@pytest.mark.parametrize("sid", PRERELEASE_IDS)
def test_predicted_transits_are_statistically_representative(prerelease_table, scanlaw, sid):
    """Pixel-centre predictions (0.9 deg pixels vs 0.7 deg FoV width) reproduce only part of the
    real transits epoch-by-epoch, but timing, scan angles and counts must be right."""
    real, ra, dec = _real_transits(prerelease_table, sid)
    t_pred, psi_pred = scanlaw.transits(ra, dec)
    dt = (real[:, 0][:, None] - t_pred[None, :]) * 365.25
    nearest = np.argmin(np.abs(dt), axis=1)
    dt_near = dt[np.arange(len(real)), nearest]
    matched = np.abs(dt_near) < 0.02                              # days
    assert matched.mean() > 0.3
    assert np.median(np.abs(dt_near[matched])) * 86400 < 120       # seconds
    dpsi = np.angle(np.exp(1j * (real[matched, 1] - psi_pred[nearest[matched]])))
    assert np.median(np.abs(np.rad2deg(dpsi))) < 1.0
    assert 0.8 < len(t_pred) / len(real) < 1.5


@pytest.mark.parametrize("sid", PRERELEASE_IDS)
def test_parallax_factor_matches_real_epochs(prerelease_table, sid):
    real, ra, dec = _real_transits(prerelease_table, sid)
    f = parallax_factor_al(real[:, 0], real[:, 1], ra, dec)
    np.testing.assert_allclose(f, real[:, 2], atol=0.02)
