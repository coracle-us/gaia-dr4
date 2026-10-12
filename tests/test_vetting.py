import numpy as np
import pytest
import pandas as pd
from astropy.table import Table

from companions.vetting import epoch_rv_vetting, summarize_dpac


def test_summarize_dpac_collects_solution_types_per_source():
    tables = {
        "gaiadr4.nss_two_body_orbit": pd.DataFrame(dict(source_id=[1, 1, 3],
                                                        nss_solution_type=["Orbital", "SB1", "Orbital"],
                                                        period=[500.0, 500.0, 1200.0])),
        "gaiadr4.nss_acceleration_astro": pd.DataFrame(dict(source_id=[2], nss_solution_type=["Acceleration7"])),
    }
    out = summarize_dpac(tables, [1, 2, 3, 4]).set_index("source_id")
    assert out.loc[1, "dpac_solutions"] == "Orbital;SB1"
    assert out.loc[1, "dpac_period_d"] == 500.0
    assert out.loc[2, "dpac_solutions"] == "Acceleration7"
    assert out.loc[4, "dpac_solutions"] == ""


def test_epoch_rv_vetting_flags_variable_host_after_inner_planet_removal():
    rng = np.random.default_rng(0)
    t = np.sort(rng.uniform(1800, 3800, 35))
    p_b, tc_b, k_b = 3.0, 2500.0, 1.5                              # d, d (obs_time_rv scale), km/s
    quiet = 7 - k_b * np.sin(2 * np.pi * (t - tc_b) / p_b) + 0.6 * rng.standard_normal(35)
    binary = quiet + 6 * np.sin(2 * np.pi * t / 600)
    tab = Table({"source_id": np.r_[np.full(35, 10), np.full(35, 20)],
                 "obs_time_rv": np.r_[t, t], "radial_velocity": np.r_[quiet, binary],
                 "radial_velocity_error": np.full(70, 0.6), "rv_validity_flag": np.ones(70, bool)})
    tables = {10: tab[tab["source_id"] == 10], 20: tab[tab["source_id"] == 20]}
    hosts = pd.DataFrame(dict(source_id=[10, 20, 30], pl_orbper=[p_b] * 3,
                              pl_tranmid=[tc_b + 2455197.5] * 3, pl_rvamp=[k_b * 1000] * 3))
    out = epoch_rv_vetting(tables, hosts).set_index("source_id")
    assert out.loc[10, "gaia_rv_pvalue"] > 1e-3 and not out.loc[10, "gaia_rv_variable"]
    assert out.loc[20, "gaia_rv_pvalue"] < 1e-6 and out.loc[20, "gaia_rv_variable"]
    assert np.isnan(out.loc[30, "gaia_rv_pvalue"])


def test_known_period_rv_test_detects_orbit_at_that_period():
    from companions.gaia_rv import rvdata_from_gaia_epoch
    from companions.vetting import rv_known_period_test
    rng = np.random.default_rng(3)
    t = np.sort(rng.uniform(1800, 3800, 30))
    p_d = 37.315
    rv_kms = -42 + 25 * np.sin(2 * np.pi * t / p_d + 0.4) + 1.0 * rng.standard_normal(30)
    tab = Table({"source_id": np.full(30, 5), "obs_time_rv": t, "radial_velocity": rv_kms,
                 "radial_velocity_error": np.full(30, 1.0), "rv_validity_flag": np.ones(30, bool)})
    res = rv_known_period_test(rvdata_from_gaia_epoch(tab, 5), p_d)
    assert res["k_ms"] / 1000 == pytest.approx(25, rel=0.1)
    assert res["dchi2"] > 100
    flat = Table({"source_id": np.full(30, 6), "obs_time_rv": t, "radial_velocity": -42 + rng.standard_normal(30),
                  "radial_velocity_error": np.full(30, 1.0), "rv_validity_flag": np.ones(30, bool)})
    assert rv_known_period_test(rvdata_from_gaia_epoch(flat, 6), p_d)["dchi2"] < 15


def _rv_table(sid, t, rv_kms, err=1.0):
    return Table({"source_id": np.full(len(t), sid), "obs_time_rv": t, "radial_velocity": rv_kms,
                  "radial_velocity_error": np.full(len(t), err), "rv_validity_flag": np.ones(len(t), bool)})


def test_known_period_rv_test_reports_k_uncertainty_and_eclipse_phase():
    from companions.gaia_rv import rvdata_from_gaia_epoch
    from companions.vetting import rv_known_period_test
    rng = np.random.default_rng(4)
    t = np.sort(rng.uniform(1800, 3800, 40))                     # obs_time_rv [d]
    p_d, t0_bjd = 37.315, 2458515.0661
    t0 = t0_bjd - 2455197.5
    # primary behind the companion at T0: RV through systemic, decreasing
    rv = -42 - 25 * np.sin(2 * np.pi * (t - t0) / p_d) + rng.standard_normal(40)
    res = rv_known_period_test(rvdata_from_gaia_epoch(_rv_table(5, t, rv), 5), p_d, t0_bjd=t0_bjd)
    assert 0 < res["k_err_ms"] < 1000 and res["k_upper_ms"] > res["k_ms"]
    assert abs(res["t0_offset_d"]) < 0.5
    shifted = rv_known_period_test(rvdata_from_gaia_epoch(_rv_table(5, t, rv), 5), p_d,
                                   t0_bjd=t0_bjd + p_d / 4)
    assert abs(abs(shifted["t0_offset_d"]) - p_d / 4) < 0.5
    assert np.isnan(rv_known_period_test(rvdata_from_gaia_epoch(_rv_table(5, t, rv), 5), p_d)["t0_offset_d"])


def test_known_period_rv_table_has_one_row_per_candidate_period():
    from companions.vetting import known_period_rv_table
    rng = np.random.default_rng(5)
    t = np.sort(rng.uniform(1800, 3800, 30))
    tables = {7: _rv_table(7, t, 3 + rng.standard_normal(30))}
    hosts = pd.DataFrame(dict(source_id=[7, 8, 9], hostname=["A", "B", "C"],
                              rv_test_periods_d=["89.9651;179.9302", "10", ""],
                              rv_t0_bjd=[2458703.8277, np.nan, np.nan]))
    out = known_period_rv_table(tables, hosts)
    a = out[out.hostname == "A"]
    assert list(a.period_d) == [89.9651, 179.9302]
    assert (a.n == 30).all() and (a.k_upper_ms > 0).all() and a.t0_offset_d.notna().all()
    b = out[out.hostname == "B"]
    assert len(b) == 1 and b.n.iloc[0] == 0                    # no epoch RVs downloaded
    assert "C" not in set(out.hostname)


def _conjunction_bjd(p_d, e, omega, tp_bjd, secondary=False):
    nu = (1.5 if secondary else 0.5) * np.pi - omega
    big_e = 2 * np.arctan(np.sqrt((1 - e) / (1 + e)) * np.tan(nu / 2))
    return tp_bjd + ((big_e - e * np.sin(big_e)) % (2 * np.pi)) / (2 * np.pi) * p_d


def test_keplerian_rv_fit_at_fixed_period_recovers_eccentric_orbit_and_eclipse_phase():
    from companions.gaia_rv import rvdata_from_gaia_epoch
    from companions.orbit import solve_kepler
    from companions.vetting import rv_keplerian_fixed_period
    rng = np.random.default_rng(6)
    p_d, e, omega, k, tp_bjd = 72.5336, 0.48, np.deg2rad(201.5), 22.7, 2457418.905   # DR3 SB1 of EPIC 211892898
    t = np.sort(rng.uniform(1800, 3800, 40))                    # obs_time_rv [d]
    big_e = solve_kepler(2 * np.pi * (t + 2455197.5 - tp_bjd) / p_d, e)
    nu = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(big_e / 2), np.sqrt(1 - e) * np.cos(big_e / 2))
    rv = 83 + k * (np.cos(nu + omega) + e * np.cos(omega)) + rng.standard_normal(40)
    data = rvdata_from_gaia_epoch(_rv_table(9, t, rv), 9)
    t_prim = _conjunction_bjd(p_d, e, omega, tp_bjd)
    res = rv_keplerian_fixed_period(data, p_d, t0_bjd=t_prim)
    assert res["k_ms"] / 1000 == pytest.approx(k, rel=0.05)
    assert res["e"] == pytest.approx(e, abs=0.05)
    assert np.rad2deg(np.angle(np.exp(1j * (np.deg2rad(res["omega_deg"]) - omega)))) == pytest.approx(0, abs=8)
    assert abs(res["t0_offset_d"]) < 1.0 and res["dchi2"] > 1000
    t_sec = _conjunction_bjd(p_d, e, omega, tp_bjd, secondary=True)
    assert abs(rv_keplerian_fixed_period(data, p_d, t0_bjd=t_sec)["t0_offset_d"]) > 10


def test_known_period_rv_table_includes_keplerian_columns():
    from companions.vetting import known_period_rv_table
    rng = np.random.default_rng(7)
    t = np.sort(rng.uniform(1800, 3800, 30))
    out = known_period_rv_table({7: _rv_table(7, t, 3 + rng.standard_normal(30))},
                                pd.DataFrame(dict(source_id=[7], hostname=["A"], rv_test_periods_d=["73.245"],
                                                  rv_t0_bjd=[2457747.385])))
    for col in ("kep_k_ms", "kep_e", "kep_omega_deg", "kep_dchi2", "kep_t0_offset_d"):
        assert col in out.columns and np.isfinite(out[col].iloc[0])
