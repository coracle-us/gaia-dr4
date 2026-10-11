"""Release-day null calibration on REAL DR4 single stars near the target hosts.

The simulated null (scripts/calibrate_null.py) assumed an excess-noise prior. On release day,
measure the detection statistics on control stars that share each host's magnitude and sky
region (so: noise, sampling, attitude systematics), are well fitted by a single star
(RUWE < 1.1) and have no DR4 NSS solution:
    python scripts/calibrate_null_real.py --k 5 --procs 30
Outputs: results/calibration_null_real.csv, results/thresholds_real.csv (compare with
results/thresholds.csv; set "thresholds_file" in config/release_day.json to adopt them).
Real "single" stars contain some genuine companions, so the real-data tail and its thresholds
are conservative (an upper bound on the false-alarm threshold).
ADQL templates and table names come from config/release_day.json ("control" section): verify
them against the DR4 archive on release day.
"""
from __future__ import annotations

import argparse
import json
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def angsep_deg(ra1, dec1, ra2, dec2):
    r1, d1, r2, d2 = map(np.deg2rad, (ra1, dec1, ra2, dec2))
    a = np.sin((d2 - d1) / 2) ** 2 + np.cos(d1) * np.cos(d2) * np.sin((r2 - r1) / 2) ** 2
    return np.rad2deg(2 * np.arcsin(np.sqrt(np.clip(a, 0, 1))))


def select_controls(cand: pd.DataFrame, hosts: pd.DataFrame, k=5, max_dg=0.5, max_sep_deg=3.0,
                    max_ruwe=1.1, seed=0) -> pd.DataFrame:
    """Up to k distinct control stars per host (columns: source_id, ra, dec, phot_g_mean_mag,
    ruwe), never a host itself, each used once."""
    rng = np.random.default_rng(seed)
    cand = cand[~cand.source_id.isin(hosts.source_id) & (cand.ruwe < max_ruwe)]
    used, rows = set(), []
    for _, h in hosts.iterrows():
        sep = angsep_deg(h.ra, h.dec, cand.ra.to_numpy(), cand.dec.to_numpy())
        ok = cand[(sep < max_sep_deg) & (np.abs(cand.phot_g_mean_mag - h.phot_g_mean_mag) <= max_dg)
                  & ~cand.source_id.isin(used)]
        pick = ok.iloc[rng.permutation(len(ok))[:k]]
        used.update(pick.source_id)
        rows.append(pick.assign(host_source_id=h.source_id))
    return pd.concat(rows, ignore_index=True) if rows else cand.iloc[:0]


def query_candidates(hosts: pd.DataFrame, ccfg: dict, cache_dir: Path) -> pd.DataFrame:
    from astroquery.gaia import Gaia
    cache_dir.mkdir(parents=True, exist_ok=True)
    parts = []
    for _, h in hosts.iterrows():
        path = cache_dir / f"{int(h.source_id)}.csv"
        if not path.exists():
            q = ccfg["candidate_adql"].format(table=ccfg["gaia_source_table"], ra=h.ra, dec=h.dec,
                                              radius=ccfg["max_sep_deg"],
                                              gmin=h.phot_g_mean_mag - ccfg["max_dg"],
                                              gmax=h.phot_g_mean_mag + ccfg["max_dg"],
                                              max_ruwe=ccfg["max_ruwe"])
            Gaia.launch_job_async(q).get_results().to_pandas().to_csv(path, index=False)
        parts.append(pd.read_csv(path))
    cand = pd.concat(parts, ignore_index=True).drop_duplicates("source_id")
    for table in ccfg.get("nss_tables", []):                         # drop known non-single stars
        ids = cand.source_id.tolist()
        nss = []
        for i in range(0, len(ids), 500):
            id_list = ",".join(str(int(s)) for s in ids[i:i + 500])
            nss.append(Gaia.launch_job(f"select source_id from {table} where source_id in ({id_list})")
                       .get_results().to_pandas())
        cand = cand[~cand.source_id.isin(pd.concat(nss).source_id)]
    return cand


def _stats(args):
    from companions.cascade import HostInfo, Thresholds, run_cascade
    from companions.epoch import EpochData
    sid, table, g, model = args
    try:
        r = run_cascade(EpochData.from_table(table, sid), HostInfo(str(sid), 1.0),
                        Thresholds(1e9, 1e9, 1e9), model=model, do_stability=False)
        return dict(source_id=sid, g_mag=g, dlnl2_orbit=r.dlnl2_orbit, dlnl2_acc7=r.dlnl2_acc7,
                    dlnl2_acc9=r.dlnl2_acc9, s_ccd_fit=r.s_ccd, s_tr_fit=r.s_tr)
    except Exception as err:
        return dict(source_id=sid, g_mag=g, error=str(err))


def main(argv=None) -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", type=Path, default=ROOT / "config/release_day.json")
    ap.add_argument("--k", type=int, default=5)
    ap.add_argument("--procs", type=int, default=30)
    args = ap.parse_args(argv)
    cfg = json.loads(args.config.read_text())
    ccfg = cfg["control"]
    import release_day
    from analyze_calibration import STATS, FAPS, threshold
    from companions.archive import DataLinkSource

    hosts = release_day.build_hosts()[["source_id", "ra", "dec", "phot_g_mean_mag"]]
    cand = query_candidates(hosts, ccfg, ROOT / "data/raw/dr4/control_candidates")
    ctl = select_controls(cand, hosts, k=args.k, max_dg=ccfg["max_dg"],
                          max_sep_deg=ccfg["max_sep_deg"], max_ruwe=ccfg["max_ruwe"])
    src = DataLinkSource(cfg["data_release"], ROOT / "data/raw/dr4/control_epoch_astrometry",
                         retrieval_type=cfg["retrieval_type"], data_structure=cfg["data_structure"],
                         batch_size=cfg["batch_size"], retries=cfg["retries"], backoff_s=cfg["backoff_s"])
    tables = src.fetch(ctl.source_id.tolist())
    jobs = [(int(r.source_id), tables[int(r.source_id)], r.phot_g_mean_mag, cfg["model"])
            for r in ctl.itertuples() if int(r.source_id) in tables]
    with Pool(args.procs) as pool:
        null = pd.DataFrame(pool.map(_stats, jobs, chunksize=2))
    null.to_csv(ROOT / "results/calibration_null_real.csv", index=False, float_format="%.6g")
    good = null.dropna(subset=["dlnl2_orbit"])
    rows = [dict(statistic=s, fap_per_star=f, threshold=threshold(good[s], f)[0],
                 method=threshold(good[s], f)[1], n_null=len(good)) for s in STATS for f in FAPS]
    thr = pd.DataFrame(rows)
    thr.to_csv(ROOT / "results/thresholds_real.csv", index=False, float_format="%.4g")
    sim = pd.read_csv(ROOT / "results/thresholds.csv")
    print(thr.merge(sim[["statistic", "fap_per_star", "threshold"]], on=["statistic", "fap_per_star"],
                    suffixes=("_real", "_sim")).to_string(index=False))


if __name__ == "__main__":
    main()
