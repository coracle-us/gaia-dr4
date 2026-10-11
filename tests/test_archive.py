import os
import zipfile
from pathlib import Path

import numpy as np
import pytest

from companions.archive import DataLinkSource, SimulatedSource, ZipCacheSource
from companions.dr4io import epochs_to_dr4_table, write_votable
from companions.epoch import EpochData
from companions.scanlaw import ScanLaw
from companions.simulate import simulate_epochs


@pytest.fixture(scope="module")
def scanlaw():
    try:
        return ScanLaw.load()
    except FileNotFoundError:
        pytest.skip("gaiascanlaw data not downloaded")


class StubGaia:
    """Mimics astroquery Gaia.load_data(dump_to_file=True): writes a zip in the cwd holding one
    RAW VOTable for all requested ids and returns {filename: [Table]}. Fails the first
    `fail_first` calls to exercise retries."""

    def __init__(self, scanlaw, fail_first=0):
        self.scanlaw, self.fail_first, self.calls = scanlaw, fail_first, []

    def load_data(self, ids, data_release, retrieval_type, data_structure, format, dump_to_file):
        self.calls.append(list(ids))
        if len(self.calls) <= self.fail_first:
            raise ConnectionError("archive busy")
        from astropy.table import vstack
        tabs = []
        for sid in ids:
            ep = simulate_epochs(self.scanlaw, 100.0, 10.0, 12.0, [0, 0, 5, 1, 1],
                                 np.random.default_rng(int(sid) % 1000), source_id=int(sid))
            tabs.append(epochs_to_dr4_table(ep, 100.0, 10.0))
        tab = vstack(tabs)
        name = f"EPOCH_ASTROMETRY-{data_release} batch.xml"
        write_votable(tab, Path(name))
        zname = f"datalink_output_{len(self.calls)}.zip"
        with zipfile.ZipFile(zname, "w") as z:
            z.write(name)
        os.remove(name)
        return {name: [tab]}


def test_datalink_source_batches_caches_and_retries(scanlaw, tmp_path):
    stub = StubGaia(scanlaw, fail_first=1)
    src = DataLinkSource(data_release="Gaia DR4", raw_dir=tmp_path, batch_size=2,
                         retries=3, backoff_s=0, client=stub)
    ids = [11, 12, 13]
    tables = src.fetch(ids)
    assert sorted(tables) == ids
    assert [len(c) for c in stub.calls] == [2, 2, 1]               # failed call retried
    assert len(list(tmp_path.glob("*.zip"))) == 2                  # raw zips kept unmodified
    ep = EpochData.from_table(tables[12], 12)
    assert len(ep.w) > 100
    # second fetch is served from the cache, no new archive calls
    n_calls = len(stub.calls)
    again = src.fetch(ids)
    assert len(stub.calls) == n_calls and sorted(again) == ids


def test_zip_cache_source_reads_raw_downloads(scanlaw, tmp_path):
    src = DataLinkSource(data_release="Gaia DR4", raw_dir=tmp_path, batch_size=5,
                         retries=1, backoff_s=0, client=StubGaia(scanlaw))
    src.fetch([21, 22])
    cache = ZipCacheSource(tmp_path)
    tables = cache.fetch([21, 22, 99])
    assert sorted(tables) == [21, 22]                              # 99 never downloaded


def test_simulated_source_writes_dr4_files_and_injects(scanlaw, tmp_path):
    import pandas as pd
    hosts = pd.DataFrame(dict(source_id=[1, 2], ra=[50.0, 250.0], dec=[10.0, -40.0],
                              phot_g_mean_mag=[11.0, 12.0], dist_pc=[80.0, 120.0]))
    orbit = dict(period=2.0, e=0.1, tp=0.0, a0=0.5, inc=np.pi / 2, omega=0.0, big_omega=1.0)
    src = SimulatedSource(scanlaw, hosts, out_dir=tmp_path, injections={2: orbit}, seed=0)
    tables = src.fetch([1, 2])
    assert (tmp_path / "1.xml").exists() and (tmp_path / "2.xml").exists()
    assert src.truth[2]["a0"] == 0.5 and 1 not in src.truth
    assert len(EpochData.from_table(tables[2], 2).w) > 100


def test_individual_files_get_source_id_from_file_name():
    from astropy.table import Table
    from companions.archive import _split_by_source
    t = Table({"transit_id": [1, 2], "g_transit_mag": [12.0, 12.1]})
    out = _split_by_source([("EPOCH_PHOTOMETRY-Gaia DR3 141810650360843008.xml", t)])
    assert list(out) == [141810650360843008] and len(out[141810650360843008]) == 2


@pytest.mark.slow
def test_datalink_plumbing_against_live_archive_dr3(tmp_path):
    """Live check of the real astroquery path with DR3 epoch photometry (DR4 not yet public)."""
    src = DataLinkSource(data_release="Gaia DR3", retrieval_type="EPOCH_PHOTOMETRY",
                         raw_dir=tmp_path, batch_size=2, retries=2, backoff_s=5)
    # (both RAW and INDIVIDUAL structures were checked by hand on 2026-10-11)
    ids = [141810650360843008, 142265642015674880]     # has_epoch_photometry = True
    tables = src.fetch(ids)
    assert sorted(tables) == sorted(ids)
    assert len(list(tmp_path.glob("*.zip"))) == 1
