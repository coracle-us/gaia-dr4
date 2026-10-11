import numpy as np
import pytest

from companions.epoch import EpochData
from companions.joint import RVData, fit_joint_ml, rv_kepler, semi_amplitude_ms
from companions.linear import single_star_design
from companions.orbit import orbit_design, solve_kepler, thiele_innes
from conftest import GAIA4

MJUP_MSUN = 9.5458e-4
AU_PER_YR_MS = 4740.470463


def test_semi_amplitude_of_sun_due_to_jupiter():
    mc = MJUP_MSUN
    a_rel = (1 + mc) ** (1 / 3) * 11.862 ** (2 / 3)
    a_star_au = a_rel * mc / (1 + mc)
    a0_mas = a_star_au * 100.0             # at 10 pc, parallax 100 mas
    k = semi_amplitude_ms(a0_mas, 100.0, 11.862, 0.0, np.pi / 2)
    assert k == pytest.approx(12.47, abs=0.05)


def test_rv_is_time_derivative_of_thiele_innes_z_coordinate():
    """z = C X + H Y with C = a sin(w) sin(i), H = a cos(w) sin(i) (Gaia DR3 AstroSpectroSB1
    convention, z positive receding); RV must equal dz/dt."""
    period, e, tp, a_au, inc, om = 1.7, 0.4, 0.3, 0.02, np.deg2rad(65), np.deg2rad(130)
    t = np.linspace(-3, 3, 2001)

    def z(tt):
        big_e = solve_kepler(2 * np.pi * (tt - tp) / period, e)
        x, y = np.cos(big_e) - e, np.sqrt(1 - e**2) * np.sin(big_e)
        return a_au * np.sin(inc) * (np.sin(om) * x + np.cos(om) * y)

    h = 1e-6
    dzdt = (z(t + h) - z(t - h)) / (2 * h) * AU_PER_YR_MS
    k = semi_amplitude_ms(a_au * 50.0, 50.0, period, e, inc)
    np.testing.assert_allclose(rv_kepler(t, period, e, tp, k, om), dzdt, atol=1e-4)


def _synthetic(prerelease_table, seed=0):
    rng = np.random.default_rng(seed)
    truth = dict(period=1.55, e=0.35, tp=0.1, a0=0.30, inc=np.deg2rad(115),
                 omega=np.deg2rad(200), big_omega=np.deg2rad(160))
    ep = EpochData.from_table(prerelease_table, GAIA4)
    ep.sigma = np.full_like(ep.w, 0.08)
    ti = np.array(thiele_innes(truth["a0"], truth["inc"], truth["omega"], truth["big_omega"]))
    star = np.array([0.1, -0.2, 13.6, -75.5, 17.9])
    ep.w = (single_star_design(ep, "5p") @ star
            + orbit_design(ep, truth["period"], truth["e"], truth["tp"]) @ ti
            + 0.09 * rng.standard_normal(len(ep.w)))
    t_rv = np.sort(rng.uniform(6.0, 7.2, 25))                  # 2023-2024, after DR4
    k = semi_amplitude_ms(truth["a0"], 13.6, truth["period"], truth["e"], truth["inc"])
    err = np.full(25, 12.0)
    inst = rng.permutation(np.array([0] * 12 + [1] * 13))      # interleaved, like real data
    rv = (rv_kepler(t_rv, truth["period"], truth["e"], truth["tp"], k, truth["omega"])
          + np.where(inst == 0, 50.0, -30.0) + err * rng.standard_normal(25))
    return ep, RVData(t=t_rv, rv=rv, err=err, inst=inst), truth


def test_joint_fit_recovers_orbit_and_breaks_node_ambiguity(prerelease_table):
    ep, rv, truth = _synthetic(prerelease_table, seed=3)
    res = fit_joint_ml(ep, rv, model="5p")
    assert res.period == pytest.approx(truth["period"], rel=0.01)
    assert res.e == pytest.approx(truth["e"], abs=0.05)
    assert np.rad2deg(res.inc) == pytest.approx(115, abs=5)
    # astrometry alone cannot tell (omega, Omega) from (omega+pi, Omega+pi); RVs can
    assert np.rad2deg(res.omega) % 360 == pytest.approx(200, abs=15)
    assert np.rad2deg(res.big_omega) % 360 == pytest.approx(160, abs=10)
    np.testing.assert_allclose(res.rv_offsets, [50.0, -30.0], atol=15)


def test_joint_lnlike_is_sum_of_astrometric_and_rv_parts(prerelease_table):
    ep, rv, _ = _synthetic(prerelease_table, seed=4)
    res = fit_joint_ml(ep, rv, model="5p")
    assert res.lnlike == pytest.approx(res.lnlike_astro + res.lnlike_rv)


@pytest.mark.slow
def test_gaia4_joint_fit_selects_published_node_branch(prerelease_table):
    """Real data: DR4 pre-release astrometry + Stefansson+2025 RVs. The RVs must pick the
    Omega ~ 160-175 deg branch (published joint: 158.6 +- 5.1), not Omega + 180."""
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    from validate_gaia4_joint import RV_FILE, read_mrt_rv
    if not RV_FILE.exists():
        pytest.skip("Gaia-4 RV table not downloaded")
    df = read_mrt_rv(RV_FILE)
    rv = RVData.from_bjd_tdb(df["bjd_tdb"], df["rv"], df["err"], list(df["inst"]))
    res = fit_joint_ml(EpochData.from_table(prerelease_table, GAIA4), rv, model="6p")
    assert 145 < np.rad2deg(res.big_omega) < 195
    assert res.k == pytest.approx(368, rel=0.25)
