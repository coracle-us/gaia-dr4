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
 4. results.csv (all hosts: detections, accelerations, upper limits, flags), injections.csv
    (simulate only), run_summary.json (counts, thresholds, config, git commit, timings).
Source ids: DR4 ids are assumed equal to DR3 ids ("dr4_source_id": "same_as_dr3"); if DR4
publishes a DR3->DR4 crossmatch, add the mapping in `map_source_ids`.
"""
from __future__ import annotations

import argparse
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


def build_hosts() -> pd.DataFrame:
    tg = pd.read_csv(ROOT / "data/targets.csv", dtype={"source_id": "Int64"})
    tg = tg.dropna(subset=["source_id", "ra", "dec", "phot_g_mean_mag", "dist_pc"])
    tg = tg.sort_values("pl_orbper").drop_duplicates("hostname")          # innermost giant
    inner = load_inner_inclinations()[["pl_name", "pl_orbincl"]]
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
    return tg.reset_index(drop=True)


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


def _worker(args):
    from companions.cascade import HostInfo, Thresholds, run_cascade
    from companions.epoch import EpochData
    sid, table, host, thr, model, do_post, post_kw = args
    try:
        ep = EpochData.from_table(table, sid)
        res = run_cascade(ep, HostInfo(**host), Thresholds(**thr), model=model,
                          do_posterior=do_post, posterior_kw=post_kw)
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
    args = ap.parse_args(argv)
    cfg = json.loads(args.config.read_text())
    t0 = time.time()
    rng = np.random.default_rng(args.seed)

    from companions.cascade import Thresholds
    thr = Thresholds.from_csv(ROOT / cfg.get("thresholds_file", "results/thresholds.csv"),
                              fap=cfg["fap_per_star"])
    hosts = build_hosts()
    if args.limit:
        hosts = hosts.head(args.limit).copy()
    hosts["dr4_source_id"] = map_source_ids(hosts, cfg["dr4_source_id"]).to_numpy()
    ids = hosts["dr4_source_id"].tolist()
    args.out.mkdir(parents=True, exist_ok=True)

    injections, model = {}, cfg["model"]
    if args.source == "datalink":
        from companions.archive import DataLinkSource
        src = DataLinkSource(cfg["data_release"], ROOT / "data/raw/dr4/epoch_astrometry",
                             retrieval_type=cfg["retrieval_type"], data_structure=cfg["data_structure"],
                             batch_size=cfg["batch_size"], retries=cfg["retries"],
                             backoff_s=cfg["backoff_s"])
    elif args.source == "cache":
        from companions.archive import ZipCacheSource
        src = ZipCacheSource(ROOT / "data/raw/dr4/epoch_astrometry")
    else:
        from companions.archive import SimulatedSource
        from companions.scanlaw import ScanLaw
        injections = choose_injections(hosts.assign(source_id=hosts["dr4_source_id"]),
                                       args.inject, args.inject_factor, rng)
        sim_hosts = hosts.assign(source_id=hosts["dr4_source_id"])
        src = SimulatedSource(ScanLaw.load(), sim_hosts, args.out / "sim_epochs",
                              injections=injections, seed=args.seed)
        model = "5p"                                  # simulated data carry no colour term
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
                     cfg["posterior"]))
    t_fit = time.time()
    procs = args.procs or cfg["procs"]
    with Pool(procs) as pool:
        rows = pool.map(_worker, jobs, chunksize=1)
    t_fit = time.time() - t_fit

    res = pd.DataFrame(rows)
    keep = ["hostname", "source_id", "phot_g_mean_mag", "dist_pc", "st_mass", "pl_name", "pl_orbper",
            "pl_orbincl", "lit_rv_followup", "ruwe", "prior_evidence"] + \
           [c for c in hosts.columns if c.startswith(("m50_", "m90_"))]
    base = hosts.assign(source_id=hosts["dr4_source_id"])[keep]
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
    print(json.dumps({k: summary[k] for k in ("source", "n_hosts", "n_fitted", "n_missing",
                                              "n_errors", "classifications", "seconds")}, default=str))
    det = out[out.classification.isin(["orbit", "acceleration"])]
    if len(det):
        cols = ["hostname", "classification", "dlnl2_orbit", "period_yr", "mass_mj",
                "min_mutual_inc_deg", "flags"]
        print(det[cols].round(2).to_string(index=False))


if __name__ == "__main__":
    main()
