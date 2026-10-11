"""Tables for the RNAAS note from a release-day results.csv.

    python scripts/make_table.py results/dr4/results.csv --out paper/

Writes
  detections.tex : AASTeX deluxetable of orbit and acceleration detections (RNAAS allows one
                   table): host, G, d, P, e, companion mass and minimum mutual inclination
                   (posterior median with 16-84% range), detection statistic, vetting notes;
  hosts.mrt      : AAS machine-readable table of ALL hosts (classification, statistics, orbit
                   where detected, 50% / 90% completeness masses at P = 3 yr) for the supplement.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.table import Table

FLAG_NOTES = {"twin_excluded_lit_rv": "t", "impostor_not_excluded": "i",
              "stellar_companion_check": "s", "stability_confirmed": "c",
              "stability_not_confirmed": "u", "high_e_check": "e", "marginal": "m",
              "brown_dwarf_mass": "b", "stellar_mass": "*"}


def _pm(row, name, fmt="{:.1f}"):
    md, lo, hi = (row.get(f"post_{name}_p{q}", np.nan) for q in (50, 16, 84))
    if not np.isfinite(md):
        return "\\nodata"
    return f"${fmt.format(md)}^{{+{fmt.format(hi - md)}}}_{{-{fmt.format(md - lo)}}}$"


def _notes(flags) -> str:
    if not isinstance(flags, str):
        return ""
    return ",".join(FLAG_NOTES[f] for f in flags.split(";") if f in FLAG_NOTES)


def latex_detection_table(res: pd.DataFrame) -> str:
    det = res[res.classification.isin(["orbit", "acceleration"])].sort_values(
        ["classification", "dlnl2_orbit"], ascending=[False, False])
    lines = [r"\begin{deluxetable*}{lrrlrrrrl}",
             r"\tablecaption{Gaia DR4 astrometric companions to transiting giant-planet hosts\label{tab:det}}",
             r"\tablehead{\colhead{Host} & \colhead{$G$} & \colhead{$d$ (pc)} & \colhead{Type} & "
             r"\colhead{$P$ (yr)} & \colhead{$e$} & \colhead{$M_c$ ($M_J$)} & "
             r"\colhead{$\Delta I_{\min}$ (deg)} & \colhead{$\Delta\ln L$; notes}}",
             r"\startdata"]
    for _, r in det.iterrows():
        if r.classification == "orbit":
            p, e = f"{r.period_yr:.2f}", f"{r.e:.2f}"
            mass, dinc = _pm(r, "mass_mj"), _pm(r, "min_mutual_inc_deg")
            stat = r.dlnl2_orbit
        else:
            p = e = mass = dinc = r"\nodata"
            stat = max(r.dlnl2_acc7, r.dlnl2_acc9)
        lines.append(f"{r.hostname} & {r.phot_g_mean_mag:.1f} & {r.dist_pc:.0f} & {r.classification} & "
                     f"{p} & {e} & {mass} & {dinc} & {stat:.0f}; {_notes(r["flags"])} \\\\")
    lines += [r"\enddata",
              r"\tablecomments{Type: orbit (Keplerian, $P<4$ yr) or acceleration (long-period). "
              r"$\Delta\ln L$: $2\Delta\ln\mathcal{L}$ versus a single star. Notes: "
              + "; ".join(f"{v} = {k.replace('_', ' ')}" for k, v in FLAG_NOTES.items())
              + r". Upper limits for all other hosts are in the machine-readable table.}",
              r"\end{deluxetable*}"]
    return "\n".join(lines)


def write_mrt(res: pd.DataFrame, path: Path) -> None:
    t = Table()
    t["Host"] = np.array(res.hostname.astype(str), dtype=str)
    t["SourceID"] = res.source_id.astype("int64")
    t["Gmag"] = res.phot_g_mean_mag.astype(float).round(3)
    t["Dist"] = res.dist_pc.astype(float).round(1)
    t["Class"] = np.array(res.classification.fillna("unfitted").astype(str), dtype=str)
    t["dlnL_orb"] = res.dlnl2_orbit.astype(float).round(2)
    t["dlnL_acc"] = np.fmax(res.dlnl2_acc7.astype(float), res.dlnl2_acc9.astype(float)).round(2)
    t["Per"] = res.period_yr.astype(float).round(3)
    t["Mass"] = res.get("post_mass_mj_p50", pd.Series(np.nan, index=res.index)).astype(float).round(2)
    # completeness masses at the solver cap (500 M_J) mean "no useful sensitivity" -> blank
    t["m50_P3"] = res.m50_mj_P3.astype(float).where(res.m50_mj_P3 < 499).round(2)
    t["m90_P3"] = res.m90_mj_P3.astype(float).where(res.m90_mj_P3 < 499).round(2)
    t["Flags"] = np.array(res["flags"].fillna("-").replace("", "-").astype(str), dtype=str)
    from astropy.table import MaskedColumn
    for c in t.colnames:                                  # NaN -> masked (blank in the MRT)
        if t[c].dtype.kind == "f" and np.isnan(t[c]).any():
            t[c] = MaskedColumn(np.nan_to_num(t[c], nan=0.0), mask=np.isnan(t[c]))
    units = dict(Gmag="mag", Dist="pc", Per="yr")          # masses in M_Jup (see descriptions)
    for c, u in units.items():
        t[c].unit = u
    descr = dict(Host="Host star name", SourceID="Gaia source identifier", Gmag="Gaia G magnitude",
                 Dist="Distance", Class="Cascade classification (orbit/acceleration/none)",
                 dlnL_orb="2 Delta lnL, orbit vs single star", dlnL_acc="2 Delta lnL, best acceleration model",
                 Per="Orbital period (orbit detections)", Mass="Companion mass, posterior median (M_Jup)",
                 m50_P3="50% completeness mass (M_Jup), aligned orbit, P = 3 yr",
                 m90_P3="90% completeness mass (M_Jup), aligned orbit, P = 3 yr", Flags="Vetting flags")
    for c, d in descr.items():
        t[c].description = d
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    t.write(path, format="ascii.mrt", overwrite=True)


def main(argv=None) -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("results", type=Path)
    ap.add_argument("--out", type=Path, default=Path("paper"))
    args = ap.parse_args(argv)
    res = pd.read_csv(args.results)
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "detections.tex").write_text(latex_detection_table(res), encoding="utf-8")
    write_mrt(res, args.out / "hosts.mrt")
    print(f"wrote {args.out / 'detections.tex'} and {args.out / 'hosts.mrt'}")


if __name__ == "__main__":
    main()
