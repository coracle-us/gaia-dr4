"""Compare posterior priors on recovered dry-run injections: flat Thiele-Innes vs Campbell.

Input : results/dryrun/{results,injections}.csv and results/dryrun/sim_epochs/<source_id>.xml
        (scripts/release_day.py --source simulate)
Output: results/posterior_bias.csv; summary of mass ratio (posterior median / truth) and 68%
        interval coverage for both samplers.
"""
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io.votable import parse_single_table

ROOT = Path(__file__).resolve().parents[1]
DRY = ROOT / "results/dryrun"


def run(args):
    from companions.epoch import EpochData
    from companions.posterior import companion_mass_mj, sample_orbit, sample_orbit_campbell
    from companions.search import search_orbit_ml
    sid, mstar, m_true, i_true = args
    ep = EpochData.from_table(parse_single_table(DRY / "sim_epochs" / f"{sid}.xml").to_table(), sid)
    ml = search_orbit_ml(ep, model="5p")
    out = dict(source_id=sid, mass_true=m_true, inc_true=i_true)
    for name, fn, kw in (("ti", sample_orbit, dict(n_walkers=32, n_steps=2000, burn=700)),
                         ("campbell", sample_orbit_campbell, dict(n_walkers=40, n_steps=3000, burn=1200))):
        post = fn(ep, ml, model="5p", seed=1, **kw)
        mass = companion_mass_mj(post.a0, post.linear[:, 2], post.period, mstar)
        lo, md, hi = np.percentile(mass, [16, 50, 84])
        ilo, imd, ihi = np.percentile(np.rad2deg(post.inc), [16, 50, 84])
        out.update({f"{name}_mass_p50": md, f"{name}_mass_in68": lo < m_true < hi,
                    f"{name}_inc_p50": imd, f"{name}_inc_in68": ilo < i_true < ihi})
    return out


def main():
    res = pd.read_csv(DRY / "results.csv")
    inj = pd.read_csv(DRY / "injections.csv").merge(res[["source_id", "classification", "st_mass"]],
                                                    on="source_id")
    inj = inj[inj.classification == "orbit"]
    jobs = [(int(r.source_id), r.st_mass if np.isfinite(r.st_mass) else 1.0, r.mass_true_mj,
             r.inc_true_deg) for r in inj.itertuples()]
    with Pool(min(30, len(jobs))) as pool:
        df = pd.DataFrame(pool.map(run, jobs))
    df.to_csv(ROOT / "results/posterior_bias.csv", index=False, float_format="%.5g")
    for name in ("ti", "campbell"):
        ratio = df[f"{name}_mass_p50"] / df.mass_true
        print(f"{name:9s}: mass ratio median {ratio.median():.3f} mean {ratio.mean():.3f}; "
              f"mass truth in 68%: {df[f'{name}_mass_in68'].mean():.2f}; "
              f"inc truth in 68%: {df[f'{name}_inc_in68'].mean():.2f} (n = {len(df)})")


if __name__ == "__main__":
    main()
