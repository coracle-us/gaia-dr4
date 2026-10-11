"""Headline population numbers from a release-day results table.

    python scripts/population_analysis.py results/dr4/results.csv --out results/dr4

Population: outer companions of 1-13 M_J (log-uniform) with P = 1-4 yr (log-uniform), the range
where DR4 orbits are reliable (Lammers & Winn 2025; our injections). Detections: classification
"orbit" with posterior mass median in range and no stellar/brown-dwarf mass flag, no
"marginal" flag. Hosts with the "stellar_companion_check" flag but no detection stay in the
sample (conservative). Completeness per host: companions.population.host_completeness with the
calibrated models (results/completeness_models.csv), N_obs from the fit, sigma_eff^2 =
sigma_ccd(G)^2 + s_ccd^2 (fitted per-CCD jitter).
Outputs: population_summary.json (occurrence for all / hot (P_b < 10 d) / warm hosts, aligned and
isotropic geometry; mutual-inclination Rayleigh width), host_completeness.csv.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from companions.population import (CompletenessModel, host_completeness, mutual_inclination_posterior,
                                   occurrence_posterior)
from companions.precision import sigma_ccd_holl2023

ROOT = Path(__file__).resolve().parents[1]
MASS_RANGE, PERIOD_RANGE = (1.0, 13.0), (1.0, 4.0)


def is_detection(r) -> bool:
    flags = r["flags"] if isinstance(r["flags"], str) else ""
    m = r.get("post_mass_mj_p50", r.get("mass_mj", np.nan))
    if not np.isfinite(m):
        m = r.get("mass_mj", np.nan)
    return (r["classification"] == "orbit" and MASS_RANGE[0] <= m <= MASS_RANGE[1]
            and PERIOD_RANGE[0] <= r["period_yr"] <= PERIOD_RANGE[1]
            and not any(f in flags for f in ("stellar_mass", "brown_dwarf_mass", "marginal")))


def main(argv=None) -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("results", type=Path)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)
    res = pd.read_csv(args.results)
    res = res[res.status == "fitted"].copy()
    model = CompletenessModel.from_csv(ROOT / "results/completeness_models.csv")
    res["detected"] = res.apply(is_detection, axis=1)
    for geom in ("aligned", "isotropic"):
        res[f"completeness_{geom}"] = [
            host_completeness(dict(n_obs=r.n_obs,
                                   sigma_eff=float(np.hypot(sigma_ccd_holl2023(r.phot_g_mean_mag), r.s_ccd)),
                                   plx=1000.0 / r.dist_pc,
                                   mstar=r.st_mass if np.isfinite(r.st_mass) else 1.0),
                              model, MASS_RANGE, PERIOD_RANGE, geometry=geom, seed=i)
            for i, r in enumerate(res.itertuples())]
    args.out.mkdir(parents=True, exist_ok=True)
    res[["hostname", "source_id", "pl_orbper", "detected", "completeness_aligned",
         "completeness_isotropic"]].to_csv(args.out / "host_completeness.csv", index=False,
                                           float_format="%.4g")
    summary = dict(population=dict(mass_mj=MASS_RANGE, period_yr=PERIOD_RANGE), occurrence={})
    groups = {"all": res, "hot (P_b < 10 d)": res[res.pl_orbper < 10],
              "warm (P_b >= 10 d)": res[res.pl_orbper >= 10]}
    for name, g in groups.items():
        for geom in ("aligned", "isotropic"):
            post = occurrence_posterior(g.detected, g[f"completeness_{geom}"])
            summary["occurrence"][f"{name}, {geom}"] = dict(
                n_hosts=len(g), n_detected=int(g.detected.sum()),
                effective_n=float(g[f"completeness_{geom}"].sum()),
                median=post.median, p16=post.quantile(0.16), p84=post.quantile(0.84),
                upper95=post.upper95)
    det = res[res.detected & np.isfinite(res.get("post_inc_deg_p50", pd.Series(np.nan, index=res.index)))]
    if len(det) >= 3:
        err = ((det.post_inc_deg_p84 - det.post_inc_deg_p16) / 2).clip(lower=0.5)
        i_b = det.pl_orbincl.fillna(90.0)
        post = mutual_inclination_posterior(det.post_inc_deg_p50, err, i_b)
        lo, md, hi = np.percentile(post.samples, [16, 50, 84])
        summary["mutual_inclination_rayleigh_sigma_deg"] = dict(n=len(det), median=md, p16=lo, p84=hi,
                                                                 upper95=float(np.percentile(post.samples, 95)))
    (args.out / "population_summary.json").write_text(json.dumps(summary, indent=1))
    for k, v in summary["occurrence"].items():
        print(f"{k:32s} n = {v['n_hosts']:3d}, detected {v['n_detected']:2d}, effective N "
              f"{v['effective_n']:6.1f}: f = {v['median']:.3f} (+{v['p84'] - v['median']:.3f} "
              f"-{v['median'] - v['p16']:.3f}), 95% < {v['upper95']:.3f}")
    if "mutual_inclination_rayleigh_sigma_deg" in summary:
        s = summary["mutual_inclination_rayleigh_sigma_deg"]
        print(f"mutual inclination Rayleigh sigma (n = {s['n']}): {s['median']:.1f} "
              f"(+{s['p84'] - s['median']:.1f} -{s['median'] - s['p16']:.1f}) deg, 95% < {s['upper95']:.1f}")


if __name__ == "__main__":
    main()
