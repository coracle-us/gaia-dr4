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
