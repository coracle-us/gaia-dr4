import numpy as np
import pytest

from companions.cascade import HostInfo, Thresholds, run_cascade
from companions.posterior import companion_mass_mj
from companions.scanlaw import ScanLaw
from companions.simulate import simulate_epochs

STAR = [0.0, 0.0, 10.0, -20.0, 5.0]           # parallax 10 mas (100 pc)
THR = Thresholds(orbit=42.6, acc7=21.7, acc9=27.1)


@pytest.fixture(scope="module")
def scanlaw():
    try:
        return ScanLaw.load()
    except FileNotFoundError:
        pytest.skip("gaiascanlaw data not downloaded")


def _host(**kw):
    base = dict(hostname="TEST-1", mstar=1.0, inner_inc_deg=88.0, lit_rv_followup=True, ruwe=1.0)
    base.update(kw)
    return HostInfo(**base)


def test_thresholds_load_from_calibration_file():
    thr = Thresholds.from_csv(fap=1e-4)
    assert 35 < thr.orbit < 50 and 15 < thr.acc7 < 30 and 20 < thr.acc9 < 35


def test_single_star_is_not_a_detection(scanlaw):
    ep = simulate_epochs(scanlaw, 120.0, 25.0, 11.5, STAR, np.random.default_rng(10),
                         s_ccd=0.02, s_tr=0.02)
    res = run_cascade(ep, _host(), THR, model="5p")
    assert res.classification == "none"


def test_aligned_orbit_detected_and_characterised(scanlaw):
    orbit = dict(period=2.0, e=0.2, tp=0.3, a0=0.35, inc=np.deg2rad(86.0),
                 omega=1.0, big_omega=2.2)
    ep = simulate_epochs(scanlaw, 120.0, 25.0, 11.5, STAR, np.random.default_rng(11),
                         orbit=orbit, s_ccd=0.02, s_tr=0.02)
    res = run_cascade(ep, _host(), THR, model="5p")
    assert res.classification == "orbit"
    assert res.period == pytest.approx(2.0, rel=0.05)
    assert np.rad2deg(res.inc) == pytest.approx(86.0, abs=8)
    assert res.min_mutual_inc_deg < 10
    true_mass = companion_mass_mj(0.35, 10.0, 2.0, 1.0)
    assert res.mass_mj == pytest.approx(true_mass, rel=0.25)
    assert res.stable


def test_long_period_companion_classified_as_acceleration(scanlaw):
    orbit = dict(period=12.0, e=0.1, tp=2.0, a0=1.5, inc=np.deg2rad(70.0),
                 omega=0.5, big_omega=1.0)
    ep = simulate_epochs(scanlaw, 200.0, -30.0, 11.0, STAR, np.random.default_rng(12),
                         orbit=orbit, s_ccd=0.02, s_tr=0.02)
    res = run_cascade(ep, _host(), THR, model="5p")
    assert res.classification == "acceleration"


def test_vetting_flags(scanlaw):
    orbit = dict(period=1.5, e=0.1, tp=0.0, a0=0.4, inc=np.deg2rad(90.0), omega=0.0, big_omega=1.0)
    ep = simulate_epochs(scanlaw, 60.0, 10.0, 11.0, STAR, np.random.default_rng(13), orbit=orbit)
    res = run_cascade(ep, _host(lit_rv_followup=False, ruwe=2.0, prior_flags="hgca_3sigma"),
                      THR, model="5p", do_stability=False)
    assert "impostor_not_excluded" in res.flags
    assert "stellar_companion_check" in res.flags
    res2 = run_cascade(ep, _host(), THR, model="5p", do_stability=False)
    assert "twin_excluded_lit_rv" in res2.flags and "stellar_companion_check" not in res2.flags


def test_result_row_has_table_columns(scanlaw):
    ep = simulate_epochs(scanlaw, 120.0, 25.0, 11.5, STAR, np.random.default_rng(14))
    row = run_cascade(ep, _host(), THR, model="5p").row()
    for key in ("hostname", "classification", "dlnl2_orbit", "dlnl2_acc7", "dlnl2_acc9",
                "period_yr", "e", "a0_mas", "inc_deg", "mass_mj", "min_mutual_inc_deg",
                "flags", "s_ccd", "s_tr"):
        assert key in row


def test_mass_regime_flags(scanlaw):
    # a0 = 3 mas at 100 pc, P = 1 yr around 1 M_sun: ~0.3 M_sun companion -> stellar
    orbit = dict(period=1.0, e=0.1, tp=0.0, a0=3.0, inc=np.deg2rad(60.0), omega=0.0, big_omega=1.0)
    ep = simulate_epochs(scanlaw, 60.0, 10.0, 11.0, STAR, np.random.default_rng(15), orbit=orbit)
    res = run_cascade(ep, _host(), THR, model="5p", do_stability=False)
    assert res.classification == "orbit" and res.mass_mj > 80
    assert "stellar_mass" in res.flags
