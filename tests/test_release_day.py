import json
import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))


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
    assert len(df) == 8
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
    assert summary["n_hosts"] == 8 and summary["source"] == "simulate"
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
    out = tmp_path / "out"
    release_day.main(["--source", "cache", "--limit", "3", "--procs", "2", "--no-posterior",
                      "--raw-dir", str(raw), "--out", str(out)])
    df = pd.read_csv(out / "results.csv").set_index("source_id")
    assert (df.status == "fitted").all()
    assert df.loc[ids[0], "dpac_solutions"] == "Orbital"
