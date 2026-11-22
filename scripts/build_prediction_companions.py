"""Merge the archive pull (results/archive_pull, MacBook run 2026-11-22, commit 6bdc5d1) with the
literature-only companions in data/prediction_companions.csv into data/prediction_companions_v2.csv.

Rules
- Archive outer candidates (P > 100 d) on transiting giant-planet hosts; drop M < 0.1 MJ (small
  Kepler planets, not giant companions).
- Orbit from ONE published solution per companion (table ps): the default solution if it has a
  periastron time, else the most recent one with a periastron time, else the default solution.
  pscomppars mixes parameters from different papers, so it is not used for elements.
- Parc+2026 (CORALIE, arXiv 2610.25152) companions not yet in the archive, or with an updated orbit
  (WASP-8, -41, -47, -53, -81, -132), come from data/prediction_companions.csv (no periastron time).
- mass_type 'true' for archive 'Mass' (TTV/photodynamical or RV+astrometry), else 'msini'.
"""
from __future__ import annotations

from pathlib import Path
import re

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
AP = ROOT / "results/archive_pull"
PARC = {"WASP-8", "WASP-41", "WASP-47", "WASP-53", "WASP-81", "WASP-132", "WASP-22", "WASP-34",
        "WASP-67", "WASP-68", "WASP-99", "WASP-185"}
TRUE_MASS_OVERRIDE = {"WASP-8 c": (26.1, "Xiao+2026 arXiv:2609.01377 (RV+HG)"),
                      "WASP-53 c": (17.9, "Xiao+2026 arXiv:2609.01377 (RV+HG)"),
                      "TOI-2537 c": (8.9, "Xiao+2026 arXiv:2609.01377 (RV+HG)")}


def ref_text(s):
    m = re.search(r">\s*([^<]+?)\s*</a>", str(s))
    return m.group(1) if m else str(s)


def main():
    comp = pd.read_csv(AP / "companions_pscomppars.csv")
    ps = pd.read_csv(AP / "companions_ps_all_solutions.csv")
    cand = comp[comp.outer_candidate & (comp.pl_bmassj >= 0.1)]
    rows = []
    for _, c in cand.iterrows():
        if c.hostname in PARC:
            continue
        sol = ps[ps.pl_name == c.pl_name].copy()
        sol["has_tp"] = sol.pl_orbtper.notna()
        sol["ts"] = pd.to_datetime(sol.rowupdate, errors="coerce")
        pick = sol[(sol.default_flag == 1) & sol.has_tp]
        if pick.empty:
            pick = sol[sol.has_tp].sort_values("ts", ascending=False)
        if pick.empty:
            pick = sol[sol.default_flag == 1]
        s = pick.iloc[0]
        true = str(s.pl_bmassprov) == "Mass"
        rows.append(dict(hostname=c.hostname, companion=c.pl_name.split()[-1],
                         mass_mj=s.pl_bmassj, mass_type="true" if true else "msini",
                         msini_mj=s.pl_msinij if np.isfinite(s.pl_msinij) else (np.nan if true else s.pl_bmassj),
                         period_d=s.pl_orbper, a_au=np.nan, ecc=s.pl_orbeccen, omega_deg=s.pl_orblper,
                         tp_bjd=s.pl_orbtper, tp_err_d=s.pl_orbtpererr1, ref=ref_text(s.pl_refname),
                         discovery=c.discoverymethod, note="NASA Exoplanet Archive ps, one solution"))
    lit = pd.read_csv(ROOT / "data/prediction_companions.csv")
    lit = lit[lit.hostname.isin(PARC)].assign(omega_deg=np.nan, tp_bjd=np.nan, tp_err_d=np.nan,
                                               discovery="Radial Velocity")
    out = pd.concat([pd.DataFrame(rows), lit], ignore_index=True)
    for name, (m, ref) in TRUE_MASS_OVERRIDE.items():
        h, k = name.rsplit(" ", 1)
        sel = (out.hostname == h) & (out.companion == k)
        out.loc[sel, ["mass_mj", "mass_type"]] = [m, "true"]
        out.loc[sel, "note"] = out.loc[sel, "note"].fillna("") + f"; true mass {ref}"
    out = out.sort_values("hostname").reset_index(drop=True)
    out.to_csv(ROOT / "data/prediction_companions_v2.csv", index=False)
    print(out[["hostname", "companion", "mass_mj", "mass_type", "msini_mj", "period_d", "ecc",
               "tp_bjd", "ref"]].to_string())


if __name__ == "__main__":
    main()
