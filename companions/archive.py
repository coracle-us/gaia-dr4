"""Sources of DR4 epoch tables, keyed by source_id: real DataLink, raw-zip cache, simulation.

DataLinkSource wraps astroquery `Gaia.load_data(..., dump_to_file=True)`, which writes the
archive's zip into the current directory and returns the parsed tables. Calls are made from
inside `raw_dir`, so the downloaded zips stay there unmodified (CLAUDE.md: never edit raw
downloads); `raw_dir/manifest.json` records which source_ids each zip holds. Batches are retried
with exponential back-off (expect archive load on release day); failed batches are reported and
their ids are simply absent from the result.
"""
from __future__ import annotations

import contextlib
import io
import json
import re
import time
import warnings
import zipfile
from pathlib import Path

import numpy as np
from astropy.io.votable import parse_single_table
from astropy.table import Table, vstack

from companions.dr4io import epochs_to_dr4_table, write_votable


def _split_by_source(named_tables) -> dict[int, Table]:
    """(file name, table) pairs -> {source_id: table}. RAW files carry a source_id column;
    INDIVIDUAL files (one source each) only have the id at the end of the file name."""
    out: dict[int, list] = {}
    for name, tab in named_tables:
        tab = tab if isinstance(tab, Table) else tab.to_table()
        if "source_id" in tab.colnames:
            sids = np.asarray(tab["source_id"])
            for sid in np.unique(sids):
                out.setdefault(int(sid), []).append(tab[sids == sid])
        else:
            m = re.search(r"(\d+)\D*$", Path(name).stem)
            if m:
                out.setdefault(int(m.group(1)), []).append(tab)
    return {k: (v[0] if len(v) == 1 else vstack(v)) for k, v in out.items()}


class ZipCacheSource:
    """Epoch tables from previously downloaded DataLink zips in raw_dir."""

    def __init__(self, raw_dir: Path):
        self.raw_dir = Path(raw_dir)

    def fetch(self, source_ids) -> dict[int, Table]:
        want = {int(s) for s in source_ids}
        tables = []
        for zpath in sorted(self.raw_dir.glob("*.zip")):
            with zipfile.ZipFile(zpath) as z:
                for name in z.namelist():
                    if name.lower().endswith((".xml", ".vot")):
                        tables.append((name, parse_single_table(io.BytesIO(z.read(name))).to_table()))
        return {k: v for k, v in _split_by_source(tables).items() if k in want}


class DataLinkSource:
    def __init__(self, data_release: str, raw_dir: Path, retrieval_type: str = "EPOCH_ASTROMETRY",
                 data_structure: str = "RAW", batch_size: int = 100, retries: int = 5,
                 backoff_s: float = 60.0, client=None):
        self.data_release, self.retrieval_type = data_release, retrieval_type
        self.data_structure, self.batch_size = data_structure, batch_size
        self.retries, self.backoff_s = retries, backoff_s
        self.raw_dir = Path(raw_dir)
        self.raw_dir.mkdir(parents=True, exist_ok=True)
        self.manifest_path = self.raw_dir / "manifest.json"
        self._client = client
        self.failed: list[int] = []

    @property
    def client(self):
        if self._client is None:
            from astroquery.gaia import Gaia
            self._client = Gaia
        return self._client

    def _manifest(self) -> dict:
        return json.loads(self.manifest_path.read_text()) if self.manifest_path.exists() else {}

    def fetch(self, source_ids) -> dict[int, Table]:
        ids = [int(s) for s in source_ids]
        manifest = self._manifest()
        cached = {s for v in manifest.values() for s in v}
        todo = [s for s in ids if s not in cached]
        for i in range(0, len(todo), self.batch_size):
            batch = todo[i:i + self.batch_size]
            for attempt in range(self.retries):
                try:
                    with contextlib.chdir(self.raw_dir):
                        before = set(Path(".").glob("*.zip"))
                        self.client.load_data(ids=batch, data_release=self.data_release,
                                              retrieval_type=self.retrieval_type,
                                              data_structure=self.data_structure,
                                              format="votable", dump_to_file=True)
                        new = sorted(set(Path(".").glob("*.zip")) - before)
                    for z in new:
                        manifest[z.name] = batch
                    self.manifest_path.write_text(json.dumps(manifest, indent=1))
                    break
                except Exception as err:                  # archive busy / timeouts
                    if attempt == self.retries - 1:
                        warnings.warn(f"batch of {len(batch)} failed after {self.retries} "
                                      f"attempts: {err}")
                        self.failed.extend(batch)
                    else:
                        time.sleep(self.backoff_s * 2**attempt)
        return ZipCacheSource(self.raw_dir).fetch(ids)


class SimulatedSource:
    """Simulated DR4 epoch tables for hosts (columns: source_id, ra, dec, phot_g_mean_mag,
    dist_pc), written as DR4-layout VOTables to out_dir and read back, optionally with injected
    orbits {source_id: orbit dict (companions.simulate)}. Excess noise s_ccd, s_tr ~ U(0, 0.05)."""

    def __init__(self, scanlaw, hosts, out_dir: Path, injections=None, seed: int = 0):
        self.scanlaw, self.hosts = scanlaw, hosts.set_index("source_id")
        self.out_dir, self.seed = Path(out_dir), seed
        self.truth = dict(injections or {})

    def fetch(self, source_ids) -> dict[int, Table]:
        from companions.simulate import simulate_epochs
        out = {}
        for sid in source_ids:
            sid = int(sid)
            if sid not in self.hosts.index:
                continue
            h = self.hosts.loc[sid]
            rng = np.random.default_rng([self.seed, sid % (2**32)])
            s_ccd, s_tr = rng.uniform(0, 0.05, 2)
            star = [0.0, 0.0, 1000.0 / h.dist_pc, rng.normal(0, 30), rng.normal(0, 30)]
            ep = simulate_epochs(self.scanlaw, h.ra, h.dec, h.phot_g_mean_mag, star, rng,
                                 orbit=self.truth.get(sid), s_ccd=s_ccd, s_tr=s_tr, source_id=sid)
            path = self.out_dir / f"{sid}.xml"
            write_votable(epochs_to_dr4_table(ep, h.ra, h.dec), path)
            out[sid] = parse_single_table(path).to_table()
        return out
