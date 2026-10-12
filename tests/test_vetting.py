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
