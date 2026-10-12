"""Release-day pipeline: epoch astrometry for every target host -> cascade -> results table.

    python scripts/release_day.py --source datalink --out results/dr4          # 2 Dec 2026
    python scripts/release_day.py --source cache    --out results/dr4          # re-run offline
    python scripts/release_day.py --source simulate --inject 20 --out results/dryrun

Steps
 1. Hosts: data/targets.csv (one row per host, innermost transiting giant) joined with
    results/target_sensitivity.csv (m50/m90 limits), results/impostor_vetting.csv (literature
    RVs), data/target_prior_astrometry.csv (RUWE, HGCA/PMa/DR3-NSS evidence) and the transiting
    planet's inclination (NASA Exoplanet Archive pl_orbincl, cached in data/raw/release_prep/).
 2. Epoch astrometry: DataLink (config/release_day.json; raw zips kept in data/raw/dr4/),
    the local zip cache, or simulated DR4 data (data/sim/, optional injected planets).
 3. companions.cascade.run_cascade per host, in parallel; posterior for orbit detections.
    DataLink/cache runs also cross-match DPAC NSS tables (TAP) and test Gaia epoch-RV
    constancy after removing the transiting planet (companions.vetting).
 4. results.csv (all hosts: detections, accelerations, upper limits, flags), injections.csv
    (simulate only), run_summary.json (counts, thresholds, config, git commit, timings).
 5. Live status while running: status.json (stage, counts), hosts.csv, progress.jsonl (one row
    per host as it finishes), diagnostics/<source_id>.json (plot data for detections and named
    targets) and the HTML dashboard dashboard/index.html (companions.dashboard; rebuilt as hosts
    finish; `python scripts/dashboard.py OUT --watch 60` rebuilds it from another terminal).
Source ids: DR4 ids are assumed equal to DR3 ids ("dr4_source_id": "same_as_dr3"); if DR4
publishes a DR3->DR4 crossmatch, add the mapping in `map_source_ids`.
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import subprocess
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config/release_day.json"
MJUP_MSUN = 9.5458e-4


def load_inner_inclinations() -> pd.DataFrame:
    path = ROOT / "data/raw/release_prep/nea_inner.csv"
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        q = ("select pl_name, hostname, pl_orbper, pl_orbincl, pl_tranmid, pl_rvamp "
             "from pscomppars where tran_flag = 1")
        r = requests.get("https://exoplanetarchive.ipac.caltech.edu/TAP/sync",
                         params={"query": q, "format": "csv"}, timeout=300)
        r.raise_for_status()
        path.write_text(r.text, encoding="utf-8")
    return pd.read_csv(path)


EXTRA_TARGETS = ROOT / "data/extra_targets.csv"


def build_hosts(limit: int | None = None) -> pd.DataFrame:
    """Archive hosts (first `limit` if given) plus the named targets in data/extra_targets.csv,
    which are always included (columns p_min_yr / p_max_yr / rv_known_period_d / note)."""
    tg = pd.read_csv(ROOT / "data/targets.csv", dtype={"source_id": "Int64"})
    tg = tg.dropna(subset=["source_id", "ra", "dec", "phot_g_mean_mag", "dist_pc"])
    tg = tg.sort_values("pl_orbper").drop_duplicates("hostname")          # innermost giant
    inner = load_inner_inclinations()[["pl_name", "pl_orbincl", "pl_tranmid", "pl_rvamp"]]
    tg = tg.merge(inner, on="pl_name", how="left")
    sens = pd.read_csv(ROOT / "results/target_sensitivity.csv", dtype={"source_id": "Int64"})
    lim_cols = [c for c in sens.columns if c.startswith(("m50_", "m90_"))]
    tg = tg.merge(sens[["hostname"] + lim_cols], on="hostname", how="left")
    vet = pd.read_csv(ROOT / "results/impostor_vetting.csv")[["hostname", "lit_rv_followup"]]
    tg = tg.merge(vet, on="hostname", how="left")
    prior = pd.read_csv(ROOT / "data/target_prior_astrometry.csv")[["hostname", "prior_evidence"]]
    tg = tg.merge(prior, on="hostname", how="left")
    tg["lit_rv_followup"] = tg["lit_rv_followup"].astype("boolean").fillna(False).astype(bool)
    tg["prior_evidence"] = tg["prior_evidence"].fillna("")
    if limit:
        tg = tg.head(limit)
    tg["named_target"] = False
    if EXTRA_TARGETS.exists():
        ex = pd.read_csv(EXTRA_TARGETS, dtype={"source_id": "Int64"})
        ex = ex[~ex.source_id.isin(tg.source_id)]
        ex["prior_evidence"] = ""
        ex["named_target"] = True
        tg = pd.concat([tg, ex], ignore_index=True)
    for col, fill in (("note", ""), ("p_min_yr", np.nan), ("p_max_yr", np.nan),
                      ("rv_known_period_d", np.nan)):
        tg[col] = tg[col].fillna(fill) if col in tg else fill
    tg["lit_rv_followup"] = tg["lit_rv_followup"].astype(bool)
    tg["epoch_rv_expected"] = tg["epoch_rv_expected"].astype(bool)
    tg["named_target"] = tg["named_target"].astype(bool)
    return tg.reset_index(drop=True)


def search_kw_for(host) -> dict:
    """Per-host period window for named targets with a photometric period."""
    if np.isfinite(host.get("p_min_yr", np.nan)) and np.isfinite(host.get("p_max_yr", np.nan)):
        return dict(p_min=float(host["p_min_yr"]), p_max=float(host["p_max_yr"]))
    return {}


def map_source_ids(hosts: pd.DataFrame, mode: str) -> pd.Series:
    if mode == "same_as_dr3":
        return hosts["source_id"].astype("int64")
    raise NotImplementedError(f"source-id mapping {mode!r}")


def choose_injections(hosts: pd.DataFrame, n: int, factor: float, rng) -> dict:
    """Aligned planets on random eligible hosts, mass = factor x m90 at the drawn period."""
    pool = hosts[np.isfinite(hosts["m90_mj_P3"]) & (hosts["m90_mj_P3"] < 30)]
    pick = pool.sample(n=min(n, len(pool)), random_state=int(rng.integers(1 << 31)))
    inj = {}
    for _, h in pick.iterrows():
        period = rng.uniform(1.0, 3.0)
        mass = factor * h["m90_mj_P3"]
        mstar = h["st_mass"] if np.isfinite(h["st_mass"]) else 1.0
        mc = mass * MJUP_MSUN
        a0 = ((mstar + mc) * period**2) ** (1 / 3) * mc / (mstar + mc) * 1000.0 / h["dist_pc"]
        inj[int(h["source_id"])] = dict(period=period, e=rng.uniform(0, 0.3),
                                         tp=rng.uniform(0, period), a0=a0,
                                         inc=np.deg2rad(rng.normal(90, 2)),
                                         omega=rng.uniform(0, 2 * np.pi),
                                         big_omega=rng.uniform(0, np.pi), mass_mj=mass)
    return inj


def dpac_and_rv_vetting(hosts: pd.DataFrame, cfg: dict, source: str, raw_root: Path) -> pd.DataFrame:
    """DPAC NSS matches (TAP) and Gaia epoch-RV constancy (DataLink) for every host.
    Failures are logged and leave the columns empty rather than stopping the run."""
    from companions.vetting import epoch_rv_vetting, summarize_dpac
    ids = hosts["dr4_source_id"].astype(int).tolist()
    out = pd.DataFrame(dict(source_id=ids))
    dpac_dir = raw_root / "dpac"
    try:
        tables = {}
        for table in cfg.get("dpac_tables", []):
            path = dpac_dir / f"{table}.csv"
            if not path.exists() and source == "datalink":
                from astroquery.gaia import Gaia
                parts = []
                for i in range(0, len(ids), 500):
                    id_list = ",".join(str(s) for s in ids[i:i + 500])
                    parts.append(Gaia.launch_job(f"select * from {table} where source_id in "
                                                 f"({id_list})").get_results().to_pandas())
                dpac_dir.mkdir(parents=True, exist_ok=True)
                pd.concat(parts).to_csv(path, index=False)
            if path.exists():
                tables[table] = pd.read_csv(path)
        out = out.merge(summarize_dpac(tables, ids), on="source_id", how="left")
    except Exception as err:
        print(f"WARNING: DPAC NSS cross-match failed: {err}")
    try:
        rv_hosts = hosts[hosts["epoch_rv_expected"].astype(bool)]
        rv_ids = rv_hosts["dr4_source_id"].astype(int).tolist()
        raw = raw_root / "epoch_rv"
        if source == "datalink":
            from companions.archive import DataLinkSource
            rv_tables = DataLinkSource(cfg["data_release"], raw, retrieval_type=cfg["rv_retrieval_type"],
                                       data_structure=cfg["data_structure"], batch_size=cfg["batch_size"],
                                       retries=cfg["retries"], backoff_s=cfg["backoff_s"]).fetch(rv_ids)
        else:
            from companions.archive import ZipCacheSource
            rv_tables = ZipCacheSource(raw).fetch(rv_ids)
        vet = epoch_rv_vetting(rv_tables, rv_hosts.assign(source_id=rv_hosts["dr4_source_id"]))
        out = out.merge(vet, on="source_id", how="left")
    except Exception as err:
        print(f"WARNING: epoch-RV vetting failed: {err}")
    return out


def _worker(args):
    from companions.cascade import HostInfo, Thresholds, run_cascade
    from companions.epoch import EpochData
    sid, table, host, thr, model, do_post, post_kw, search_kw, diag_dir, named = args
    try:
        ep = EpochData.from_table(table, sid)
        res = run_cascade(ep, HostInfo(**host), Thresholds(**thr), model=model,
                          do_posterior=do_post, posterior_kw=post_kw, search_kw=search_kw,
                          diagnostics=True if named else None)
        if res.diagnostics:
            (Path(diag_dir) / f"{sid}.json").write_text(json.dumps(res.diagnostics))
        row = res.row()
        row.update(status="fitted", n_obs=len(ep.w), n_transit=int(ep.transit.max() + 1))
    except Exception as err:                          # one bad host must not stop the run
        row = dict(hostname=host["hostname"], status=f"error: {type(err).__name__}: {err}")
    row["source_id"] = sid
    return row


def git_commit() -> str:
    try:
        return subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True,
                              text=True, check=True).stdout.strip()
    except Exception:
        return "unknown"


def _jsonable(row: dict) -> dict:
    out = {}
    for k, v in row.items():
        if isinstance(v, (np.floating, float)):
            v = float(v) if np.isfinite(v) else None
        elif isinstance(v, np.integer):
            v = int(v)
        elif isinstance(v, np.bool_):
            v = bool(v)
        out[k] = v
    return out


class RunStatus:
    """status.json + dashboard refresh (at most once a minute unless forced). A dashboard
    failure is reported, never fatal."""

    def __init__(self, out: Path, **fields):
        self.out, self._last_build = out, 0.0
        self.state = dict(started=dt.datetime.now().isoformat(timespec="seconds"), **fields)

    def update(self, build=False, **fields):
        self.state.update(fields, updated=dt.datetime.now().isoformat(timespec="seconds"))
        (self.out / "status.json").write_text(json.dumps(self.state, indent=1, default=str))
        if build or time.time() - self._last_build > 60:
            self._last_build = time.time()
            try:
                from companions.dashboard import build_dashboard
                build_dashboard(self.out)
            except Exception as err:
                print(f"WARNING: dashboard not rebuilt: {err}")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source", choices=["datalink", "cache", "simulate"], required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--config", type=Path, default=CONFIG)
    ap.add_argument("--limit", type=int, default=None, help="first N hosts only (testing)")
    ap.add_argument("--inject", type=int, default=0, help="simulate: number of injected planets")
    ap.add_argument("--inject-factor", type=float, default=1.5, help="injected mass / m90(P=3 yr)")
    ap.add_argument("--procs", type=int, default=None)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--no-posterior", action="store_true")
    ap.add_argument("--raw-dir", type=Path, default=ROOT / "data/raw/dr4",
                    help="root of raw DR4 downloads (epoch_astrometry/, epoch_rv/, dpac/)")
    args = ap.parse_args(argv)
    cfg = json.loads(args.config.read_text())
    t0 = time.time()
    rng = np.random.default_rng(args.seed)

    from companions.cascade import Thresholds
    thr = Thresholds.from_csv(ROOT / cfg.get("thresholds_file", "results/thresholds.csv"),
                              fap=cfg["fap_per_star"])
    hosts = build_hosts(limit=args.limit)
    hosts["dr4_source_id"] = map_source_ids(hosts, cfg["dr4_source_id"]).to_numpy()
    ids = hosts["dr4_source_id"].tolist()
    args.out.mkdir(parents=True, exist_ok=True)
    keep = ["hostname", "source_id", "phot_g_mean_mag", "dist_pc", "st_mass", "pl_name", "pl_orbper",
            "pl_orbincl", "lit_rv_followup", "ruwe", "prior_evidence", "note", "rv_known_period_d",
            "named_target"] + [c for c in hosts.columns if c.startswith(("m50_", "m90_"))]
    base = hosts.assign(source_id=hosts["dr4_source_id"])[keep]
    base.to_csv(args.out / "hosts.csv", index=False, float_format="%.6g")
    diag_dir = args.out / "diagnostics"
    diag_dir.mkdir(exist_ok=True)
    for old in diag_dir.glob("*.json"):              # derived plot data from an earlier run here
        old.unlink()
    progress = args.out / "progress.jsonl"
    progress.write_text("")
    status = RunStatus(args.out, source=args.source, n_hosts=len(hosts), git_commit=git_commit(),
                       thresholds=thr.__dict__)

    injections, model = {}, cfg["model"]
    if args.source == "datalink":
        from companions.archive import DataLinkSource
        src = DataLinkSource(cfg["data_release"], args.raw_dir / "epoch_astrometry",
                             retrieval_type=cfg["retrieval_type"], data_structure=cfg["data_structure"],
                             batch_size=cfg["batch_size"], retries=cfg["retries"],
                             backoff_s=cfg["backoff_s"])
    elif args.source == "cache":
        from companions.archive import ZipCacheSource
        src = ZipCacheSource(args.raw_dir / "epoch_astrometry")
    else:
        from companions.archive import SimulatedSource
        from companions.scanlaw import ScanLaw
        injections = choose_injections(hosts.assign(source_id=hosts["dr4_source_id"]),
                                       args.inject, args.inject_factor, rng)
        sim_hosts = hosts.assign(source_id=hosts["dr4_source_id"])
        src = SimulatedSource(ScanLaw.load(), sim_hosts, args.out / "sim_epochs",
                              injections=injections, seed=args.seed)
        model = "5p"                                  # simulated data carry no colour term
    status.update(build=True, stage="fetch",
                  fetch_manifest=str(src.manifest_path) if args.source == "datalink" else None)
    t_fetch = time.time()
    tables = src.fetch(ids)
    t_fetch = time.time() - t_fetch

    jobs = []
    for _, h in hosts.iterrows():
        sid = int(h["dr4_source_id"])
        if sid not in tables:
            continue
        inc_b = h["pl_orbincl"] if np.isfinite(h.get("pl_orbincl", np.nan)) else 90.0
        host = dict(hostname=h["hostname"], mstar=h["st_mass"] if np.isfinite(h["st_mass"]) else 1.0,
                    inner_inc_deg=float(inc_b), lit_rv_followup=bool(h["lit_rv_followup"]),
                    ruwe=float(h["ruwe"]) if np.isfinite(h["ruwe"]) else np.nan,
                    prior_flags=str(h["prior_evidence"]))
        jobs.append((sid, tables[sid], host, thr.__dict__, model, not args.no_posterior,
                     cfg["posterior"], search_kw_for(h), str(diag_dir), bool(h["named_target"])))
    status.update(build=True, stage="fit", n_jobs=len(jobs), n_done=0, n_fetched=len(tables),
                  fit_started=dt.datetime.now().isoformat(timespec="seconds"))
    t_fit = time.time()
    procs = args.procs or cfg["procs"]
    rows = []
    with Pool(procs) as pool, open(progress, "a", encoding="utf-8") as prog:
        for row in pool.imap_unordered(_worker, jobs, chunksize=1):
            rows.append(row)
            prog.write(json.dumps(_jsonable(row)) + "\n")
            prog.flush()
            status.update(n_done=len(rows))
    t_fit = time.time() - t_fit
    status.update(build=True, stage="vetting")

    res = pd.DataFrame(rows)
    if args.source in ("datalink", "cache"):
        res = res.merge(dpac_and_rv_vetting(hosts, cfg, args.source, args.raw_dir),
                        on="source_id", how="left")
    out = base.merge(res.drop(columns=["hostname"]), on="source_id", how="left")
    out["status"] = out["status"].fillna("no epoch data")
    out.to_csv(args.out / "results.csv", index=False, float_format="%.6g")
    if injections:
        inj = pd.DataFrame([dict(source_id=k, period_true_yr=v["period"], e_true=v["e"],
                                 a0_true_mas=v["a0"], inc_true_deg=np.rad2deg(v["inc"]),
                                 mass_true_mj=v["mass_mj"]) for k, v in injections.items()])
        inj.to_csv(args.out / "injections.csv", index=False, float_format="%.6g")

    counts = out["classification"].value_counts().to_dict()
    summary = dict(source=args.source, n_hosts=len(out), n_fitted=int((out.status == "fitted").sum()),
                   n_missing=int((out.status == "no epoch data").sum()),
                   n_errors=int(out.status.str.startswith("error").sum()),
                   classifications=counts, thresholds=thr.__dict__, config=cfg, model=model,
                   git_commit=git_commit(), seconds=dict(fetch=t_fetch, fit=t_fit, total=time.time() - t0),
                   n_injected=len(injections))
    (args.out / "run_summary.json").write_text(json.dumps(summary, indent=1, default=str))
    status.update(build=True, stage="done", seconds=summary["seconds"])
    print(json.dumps({k: summary[k] for k in ("source", "n_hosts", "n_fitted", "n_missing",
                                              "n_errors", "classifications", "seconds")}, default=str))
    det = out[out.classification.isin(["orbit", "acceleration"])]
    if len(det):
        cols = ["hostname", "classification", "dlnl2_orbit", "period_yr", "mass_mj",
                "min_mutual_inc_deg", "flags"]
        print(det[cols].round(2).to_string(index=False))


if __name__ == "__main__":
    main()
