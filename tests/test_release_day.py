import json
import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

N_NAMED = len(pd.read_csv(Path(__file__).resolve().parents[1] / "data/extra_targets.csv"))


@pytest.fixture(scope="module")
def dryrun(tmp_path_factory):
    try:
        from companions.scanlaw import ScanLaw
        ScanLaw.load()
    except FileNotFoundError:
        pytest.skip("gaiascanlaw data not downloaded")
    import release_day
    out = tmp_path_factory.mktemp("dryrun")
    release_day.main(["--source", "simulate", "--limit", "8", "--inject", "2",
                      "--inject-factor", "3", "--procs", "4", "--no-posterior",
                      "--out", str(out), "--seed", "3"])
    return out


def test_results_table_has_every_host(dryrun):
    df = pd.read_csv(dryrun / "results.csv")
    assert len(df) == 8 + N_NAMED                 # 8 archive hosts + the named extra targets
    for col in ("hostname", "source_id", "classification", "dlnl2_orbit", "flags",
                "m50_mj_P3", "status"):
        assert col in df.columns
    assert set(df.status) == {"fitted"}


def test_injected_planets_are_recovered(dryrun):
    df = pd.read_csv(dryrun / "results.csv")
    inj = pd.read_csv(dryrun / "injections.csv")
    assert len(inj) == 2
    m = inj.merge(df, on="source_id")
    assert (m.classification == "orbit").all()
    assert ((m.period_yr / m.period_true_yr - 1).abs() < 0.1).all()


def test_run_summary_written(dryrun):
    summary = json.loads((dryrun / "run_summary.json").read_text())
    assert summary["n_hosts"] == 8 + N_NAMED and summary["source"] == "simulate"
    assert "git_commit" in summary and "thresholds" in summary


def test_cache_mode_runs_on_raw_downloads_with_dpac_vetting(tmp_path):
    """Offline release-day run from raw DataLink zips (stub archive) + a DPAC table."""
    try:
        from companions.scanlaw import ScanLaw
        sl = ScanLaw.load()
    except FileNotFoundError:
        pytest.skip("gaiascanlaw data not downloaded")
    import release_day
    from companions.archive import DataLinkSource
    from test_archive import StubGaia
    hosts = release_day.build_hosts().head(3)
    ids = hosts["source_id"].astype(int).tolist()
    raw = tmp_path / "raw"
    DataLinkSource("Gaia DR4", raw / "epoch_astrometry", batch_size=10, retries=1, backoff_s=0,
                   client=StubGaia(sl)).fetch(ids)
    (raw / "dpac").mkdir(parents=True)
    pd.DataFrame(dict(source_id=[ids[0]], nss_solution_type=["Orbital"], period=[700.0])).to_csv(
        raw / "dpac" / "gaiadr4.nss_two_body_orbit.csv", index=False)
    pd.DataFrame(dict(source_id=[ids[1]], has_epoch_rvs=[True], fraction_double_lined=[0.4])).to_csv(
        raw / "dpac" / "gaiadr4.all_source_rvs.csv", index=False)
    out = tmp_path / "out"
    release_day.main(["--source", "cache", "--limit", "3", "--procs", "2", "--no-posterior",
                      "--raw-dir", str(raw), "--out", str(out)])
    df = pd.read_csv(out / "results.csv").set_index("source_id")
    assert (df.loc[ids, "status"] == "fitted").all()
    # named extra targets were never downloaded by the stub archive
    assert (df[~df.index.isin(ids)].status == "no epoch data").all()
    assert df.loc[ids[0], "dpac_solutions"] == "Orbital"
    assert df.loc[ids[1], "fraction_double_lined"] == pytest.approx(0.4)    # DR4 SB2 indicator
    assert (out / "named_rv.csv").exists()


def test_named_extra_targets_always_included_with_period_windows():
    import release_day
    hosts = release_day.build_hosts(limit=3)
    names = set(hosts.hostname)
    assert {"TIC 147035057", "TIC 159614322"} <= names
    assert len(hosts) == 3 + N_NAMED
    row = hosts.set_index("hostname").loc["TIC 147035057"]
    assert 0.09 < row.p_min_yr < 0.1022 < row.p_max_yr < 0.11
    assert row.rv_known_period_d == pytest.approx(37.31497)
    assert release_day.search_kw_for(row) == dict(p_min=row.p_min_yr, p_max=row.p_max_yr)
    assert release_day.search_kw_for(hosts.set_index("hostname").loc["TIC 159614322"]) == {}


def test_named_targets_are_marked():
    import release_day
    hosts = release_day.build_hosts(limit=3).set_index("hostname")
    assert hosts.loc["TIC 147035057", "named_target"] and hosts.loc["TIC 159614322", "named_target"]
    assert hosts["named_target"].sum() == N_NAMED


def test_run_writes_live_progress_and_dashboard(dryrun):
    status = json.loads((dryrun / "status.json").read_text())
    assert status["stage"] == "done" and status["n_jobs"] == 8 + N_NAMED
    lines = (dryrun / "progress.jsonl").read_text().splitlines()
    assert len(lines) == 8 + N_NAMED and all(json.loads(x)["status"] == "fitted" for x in lines)
    assert len(pd.read_csv(dryrun / "hosts.csv")) == 8 + N_NAMED
    inj = pd.read_csv(dryrun / "injections.csv")
    df = pd.read_csv(dryrun / "results.csv")
    named = df[df.named_target].source_id
    for sid in list(inj.source_id) + list(named):
        assert (dryrun / "diagnostics" / f"{sid}.json").exists()
    page = (dryrun / "dashboard" / "index.html").read_text(encoding="utf-8")
    assert "TIC 147035057" in page and "Injected planets" in page


def test_named_targets_carry_rv_test_periods_and_reference_epochs():
    import release_day
    hosts = release_day.build_hosts(limit=1).set_index("hostname")
    assert hosts.loc["TIC 159614322", "rv_test_periods_d"] == "89.9651;179.9302"
    assert hosts.loc["TIC 147035057", "rv_t0_bjd"] == pytest.approx(2458515.0661)
