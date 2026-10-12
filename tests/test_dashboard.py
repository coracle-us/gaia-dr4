import json

import numpy as np
import pandas as pd
import pytest

from companions.dashboard import build_dashboard, overall_status, vetting_checklist

THR = dict(orbit=42.6, acc7=21.7, acc9=27.1)


def _orbit_row(**kw):
    row = dict(hostname="HAT-P-99", source_id=1, classification="orbit", dlnl2_orbit=120.0,
               dlnl2_acc7=3.0, dlnl2_acc9=4.0, period_yr=2.1, e=0.1, mass_mj=4.0,
               flags="twin_excluded_lit_rv;stability_confirmed", gaia_rv_n=20,
               gaia_rv_variable=False, gaia_rv_pvalue=0.4, dpac_solutions="")
    row.update(kw)
    return row


def _checks(row, manual=None):
    return {c["check"]: c for c in vetting_checklist(row, THR, manual or {})}


def test_clean_strong_orbit_waits_only_on_independent_rederivation():
    c = _checks(_orbit_row())
    assert c["significance"]["status"] == "pass"
    assert c["stability"]["status"] == "pass"
    assert c["impostor_twin"]["status"] == "pass"
    assert c["stellar_companion"]["status"] == "pass"
    assert c["gaia_epoch_rv"]["status"] == "pass"
    assert c["rederivation"]["status"] == "pending"
    assert overall_status(c.values()) == "pending"


def test_problem_flags_map_to_warnings_and_failures():
    c = _checks(_orbit_row(dlnl2_orbit=45.0, gaia_rv_variable=True, dpac_solutions="Orbital",
                           flags="impostor_not_excluded;stellar_companion_check;"
                                 "stability_not_confirmed;marginal;high_e_check;stellar_mass"))
    assert c["significance"]["status"] == "warn"
    assert c["stability"]["status"] == "warn"
    assert c["impostor_twin"]["status"] == "pending"
    assert c["stellar_companion"]["status"] == "warn"
    assert c["gaia_epoch_rv"]["status"] == "warn"
    assert c["eccentricity"]["status"] == "warn"
    assert c["mass"]["status"] == "fail"
    assert c["dpac"]["status"] == "info"
    assert overall_status(c.values()) == "fail"


def test_manual_decisions_override_automatic_status():
    manual = {("HAT-P-99", "rederivation"): ("pass", "fresh script, P = 2.08 yr"),
              ("HAT-P-99", "impostor_twin"): ("pass", "CCF single-lined")}
    c = _checks(_orbit_row(flags="impostor_not_excluded;stability_confirmed"), manual)
    assert c["rederivation"]["status"] == "pass"
    assert "fresh script" in c["rederivation"]["detail"]
    assert c["impostor_twin"]["status"] == "pass"
    assert overall_status(c.values()) == "pass"


def test_missing_epoch_rvs_are_pending_not_pass():
    c = _checks(_orbit_row(gaia_rv_n=np.nan, gaia_rv_variable=np.nan))
    assert c["gaia_epoch_rv"]["status"] == "pending"


def _diag(n=30, period=2.0):
    t = np.linspace(-2.5, 2.5, n)
    model = 0.3 * np.sin(2 * np.pi * t / period)
    tr = dict(t=t, phase=np.mod(t / period, 1), err=np.full(n, 0.05), orbit_resid=model + 0.01,
              orbit_model=model, accel_resid=np.zeros(n), accel_model=np.zeros(n))
    ang = np.linspace(0, 2 * np.pi, 50)
    return dict(transits={k: list(map(float, v)) for k, v in tr.items()},
                orbit=dict(period_yr=period, e=0.1, tp=0.0, dlnl2=120.0, s_ccd=0.01, s_tr=0.02,
                           a0_mas=0.3, inc_deg=88.0),
                accel=dict(order=7, dlnl2=3.0),
                period_profile=dict(period_yr=list(np.linspace(0.1, 4, 80)),
                                    dchi2=list(np.linspace(0, 100, 80)), best_period_yr=period),
                sky_orbit=dict(ra_mas=list(0.3 * np.cos(ang)), dec_mas=list(0.2 * np.sin(ang))),
                thresholds=THR)


def _run_dir(tmp_path, stage, n_done=2):
    out = tmp_path / "run"
    (out / "diagnostics").mkdir(parents=True)
    hosts = pd.DataFrame(dict(hostname=["HAT-P-99", "WASP-<x>", "TIC 147035057"],
                              source_id=[1, 2, 3], phot_g_mean_mag=[10.0, 11.0, 8.8],
                              dist_pc=[100.0, 200.0, 166.0], named_target=[False, False, True],
                              m50_mj_P3=[3.0, 8.0, np.nan]))
    hosts.to_csv(out / "hosts.csv", index=False)
    rows = [dict(_orbit_row(), status="fitted"),
            dict(hostname="WASP-<x>", source_id=2, classification="none", dlnl2_orbit=10.0,
                 dlnl2_acc7=1.0, dlnl2_acc9=2.0, flags="", status="fitted"),
            dict(hostname="TIC 147035057", source_id=3, classification="none", dlnl2_orbit=5.0,
                 dlnl2_acc7=1.0, dlnl2_acc9=1.0, period_yr=np.nan, flags="", status="fitted")]
    with open(out / "progress.jsonl", "w") as f:
        for r in rows[:n_done]:
            f.write(json.dumps({k: (None if isinstance(v, float) and np.isnan(v) else v)
                                for k, v in r.items()}) + "\n")
    (out / "diagnostics" / "1.json").write_text(json.dumps(_diag()))
    named = _diag(period=37.31 / 365.25)
    named["orbit"]["a0_mas"] = 0.0123
    (out / "diagnostics" / "3.json").write_text(json.dumps(named))
    status = dict(stage=stage, source="simulate", n_hosts=3, n_jobs=3, n_done=n_done,
                  started="2026-12-02T12:00:00", updated="2026-12-02T12:10:00",
                  git_commit="abc123", thresholds=THR)
    (out / "status.json").write_text(json.dumps(status))
    if stage == "done":
        pd.DataFrame(rows).merge(hosts[["source_id", "phot_g_mean_mag", "dist_pc", "named_target",
                                        "m50_mj_P3"]], on="source_id").to_csv(
            out / "results.csv", index=False)
    return out


def test_in_progress_page_shows_progress_and_refreshes(tmp_path):
    out = _run_dir(tmp_path, stage="fit", n_done=2)
    page = build_dashboard(out).read_text(encoding="utf-8")
    assert 'http-equiv="refresh"' in page
    assert "2 / 3" in page                                      # hosts fitted so far
    assert "HAT-P-99" in page                                   # candidate already listed
    assert "WASP-&lt;x&gt;" in page and "WASP-<x>" not in page  # names are escaped


def test_finished_page_has_candidate_plots_named_targets_and_no_refresh(tmp_path):
    out = _run_dir(tmp_path, stage="done", n_done=3)
    preds = tmp_path / "preds.csv"
    pd.DataFrame(dict(target=["TIC 147035057"], hypothesis=["H1 single-lined star"],
                      period_d=[37.315], a0_mas=[0.49], dchi2_exp=[3959.0],
                      k1_kms=[24.9])).to_csv(preds, index=False)
    manual = tmp_path / "manual.csv"
    pd.DataFrame(dict(hostname=["HAT-P-99"], check=["rederivation"], status=["pass"],
                      note=["independent fit agrees"])).to_csv(manual, index=False)
    page = build_dashboard(out, manual_path=manual, predictions_path=preds).read_text(encoding="utf-8")
    assert 'http-equiv="refresh"' not in page
    for name in ("period_profile", "residuals_time", "residuals_phase", "sky_orbit"):
        assert (out / "dashboard" / "plots" / f"1_{name}.png").exists()
        assert f'alt="1_{name}.png" src="data:image/png;base64,' in page   # self-contained page
    assert "independent fit agrees" in page
    assert "H1 single-lined star" in page                     # named-target predictions
    assert "TIC 147035057" in page
    assert "37.31" in page and "0.012" in page               # best orbit in its period window


def test_epoch_rv_check_not_applicable_when_vetting_was_not_run():
    row = _orbit_row()
    for k in ("gaia_rv_n", "gaia_rv_variable", "gaia_rv_pvalue"):
        row.pop(k)
    assert _checks(row)["gaia_epoch_rv"]["status"] == "na"


def test_double_lined_host_fails_the_twin_check():
    c = _checks(_orbit_row(fraction_double_lined=0.6, has_epoch_rvs=True))
    assert c["impostor_twin"]["status"] == "fail"
    single = _checks(_orbit_row(fraction_double_lined=0.0, has_epoch_rvs=True,
                                flags="impostor_not_excluded;stability_confirmed"))
    assert single["impostor_twin"]["status"] == "pending"
    assert "single-lined" in single["impostor_twin"]["detail"]


def test_named_target_shows_rv_tests_at_known_periods(tmp_path):
    out = _run_dir(tmp_path, stage="done", n_done=3)
    pd.DataFrame(dict(source_id=[3, 3], hostname=["TIC 147035057"] * 2, period_d=[37.31497, 74.6],
                      n=[40, 40], k_ms=[24900.0, 300.0], k_err_ms=[200.0, 200.0],
                      k_upper_ms=[25500.0, 900.0], dchi2=[3000.0, 2.0],
                      t0_offset_d=[0.3, np.nan], kep_k_ms=[25100.0, 400.0], kep_e=[0.437, 0.1],
                      kep_omega_deg=[201.0, 10.0], kep_dchi2=[3100.0, 3.0],
                      kep_t0_offset_d=[-0.9, 5.0])).to_csv(out / "named_rv.csv", index=False)
    page = build_dashboard(out).read_text(encoding="utf-8")
    assert "37.315" in page and "24900" in page
    assert "0.437" in page and "25100" in page                # Keplerian fit at the fixed period
