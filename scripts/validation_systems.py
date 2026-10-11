"""Expected DR4 signal for known outer companions of transiting giant-planet hosts.

Inputs : data/known_outer_companions.csv (hand-curated from the literature, refs per row;
         masses are M sin i unless mass_type == "true", so signals are lower limits there)
         data/targets.csv (scripts/build_targets.py: distance, sigma_fov, N_DR4, M_star)
Output : data/validation_known_companions.csv

dchi2_edge / dchi2_face: expected Delta-chi2 for an edge-on / face-on circular orbit,
N alpha^2 (1 + cos^2 i) / (4 sigma^2); detection threshold 50 (Lammers & Winn 2025).
Orbits with P > 4 yr are partly absorbed by the 5-parameter fit (flag p_gt_4yr); these
should appear in DR4 as accelerations / long-period trends rather than full orbits;
orbit_recoverable_edge = astrometry-only orbit expected (dchi2_edge > 50 and P <= 4 yr).
With the RV period as a prior (joint fit), P > 4 yr signals with large dchi2 are usable too.
"""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
MJUP_MSUN = 9.5458e-4

kc = pd.read_csv(ROOT / "data" / "known_outer_companions.csv")
t = (pd.read_csv(ROOT / "data" / "targets.csv", dtype={"source_id": "Int64"})
     .drop_duplicates("hostname"))
df = kc.merge(t[["hostname", "source_id", "pl_orbper", "dist_pc", "phot_g_mean_mag", "ruwe",
                 "mstar_msun", "sigma_fov_mas", "n_fov_dr4_est", "epoch_rv_expected"]],
              on="hostname", how="left")

mtot = df["mstar_msun"] + df["mass_mj"] * MJUP_MSUN
df["a_au"] = df["a_au"].fillna((mtot * (df["period_d"] / 365.25) ** 2) ** (1 / 3))
df["period_d"] = df["period_d"].fillna(365.25 * np.sqrt(df["a_au"] ** 3 / mtot))
df["p_gt_4yr"] = df["period_d"] > 4 * 365.25
alpha = df["mass_mj"] * MJUP_MSUN / df["mstar_msun"] * df["a_au"] / df["dist_pc"] * 1e3  # mas
df["alpha_uas"] = (alpha * 1e3).round(1)
base = df["n_fov_dr4_est"] * alpha**2 / (4 * df["sigma_fov_mas"] ** 2)
df["dchi2_edge"] = base.round(1)
df["dchi2_face"] = (2 * base).round(1)
df["orbit_recoverable_edge"] = (df["dchi2_edge"] > 50) & ~df["p_gt_4yr"]

out = ROOT / "data" / "validation_known_companions.csv"
df.to_csv(out, index=False, float_format="%.6g")
cols = ["hostname", "companion", "mass_mj", "mass_type", "period_d", "a_au", "dist_pc",
        "alpha_uas", "dchi2_edge", "dchi2_face", "p_gt_4yr", "orbit_recoverable_edge"]
print(df[cols].round(1).to_string(index=False))
print(f"not in targets: {sorted(df.loc[df['source_id'].isna(), 'hostname'])}")
