"""Pull orbital elements of known outer companions to transiting giant-planet hosts (MacBook run).

Cloud sessions cannot reach these archives, so this runs on Jon's MacBook. Read-only queries:
 1. NASA Exoplanet Archive pscomppars: every planet in a system that hosts one of the transiting
    giants in data/targets.csv, when that system has more than one planet.
 2. NASA Exoplanet Archive ps (all solutions) for those companions: periastron time, argument of
    periastron, eccentricity, msini and reference, so DR4 predictions can use the orbital phase.
 3. VizieR: RV-trend tables of Knutson+2014 (J/ApJ/785/126) and Bryan+2016 (J/ApJ/821/89).
    Table names are guessed from the bibcode; the script lists the catalogue's tables first and
    dumps all of them, so nothing is lost if the guess is wrong.
Outputs go to results/archive_pull/ (raw query results; never edit them by hand).
Usage: python scripts/pull_companion_params.py
"""
from __future__ import annotations

import io
from pathlib import Path

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "archive_pull"
ARCHIVE_TAP = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync"
COLS_COMP = ("pl_name,hostname,sy_pnum,discoverymethod,tran_flag,rv_flag,ttv_flag,pl_controv_flag,"
             "pl_orbper,pl_orbpererr1,pl_orbpererr2,pl_orbsmax,pl_bmassj,pl_bmassprov,pl_msinij,"
             "pl_orbeccen,pl_orbeccenerr1,pl_orblper,pl_orbtper,pl_orbtpererr1,pl_orbincl,pl_radj,"
             # pscomppars has no pl_refname (ORA-00904); it carries per-parameter reference links
             "pl_orbper_reflink,pl_bmassj_reflink,pl_orbeccen_reflink,pl_orbtper_reflink,disc_refname")
COLS_PS = ("pl_name,hostname,default_flag,pl_orbper,pl_orbpererr1,pl_orbeccen,pl_orbeccenerr1,"
           "pl_orblper,pl_orblpererr1,pl_orbtper,pl_orbtpererr1,pl_msinij,pl_bmassj,pl_bmassprov,"
           "pl_orbincl,pl_rvamp,pl_refname,rowupdate")


def tap(q: str) -> pd.DataFrame:
    r = requests.get(ARCHIVE_TAP, params={"query": q, "format": "csv"}, timeout=300)
    r.raise_for_status()
    return pd.read_csv(io.StringIO(r.text))


def quoted(names) -> str:
    return ",".join("'" + n.replace("'", "''") + "'" for n in names)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    hosts = sorted(pd.read_csv(ROOT / "data/targets.csv").hostname.dropna().unique())
    comp = []
    for i in range(0, len(hosts), 100):
        comp.append(tap(f"select {COLS_COMP} from pscomppars where sy_pnum > 1 "
                        f"and hostname in ({quoted(hosts[i:i + 100])})"))
    comp = pd.concat(comp, ignore_index=True)
    comp.to_csv(OUT / "pscomppars_multi.csv", index=False)
    # companions: not the transiting giant itself; keep everything, flag the likely outer ones
    giants = set(pd.read_csv(ROOT / "data/targets.csv").pl_name)
    outer = comp[~comp.pl_name.isin(giants)].copy()
    outer["outer_candidate"] = (outer.pl_orbper > 100) | outer.pl_orbper.isna()
    outer.to_csv(OUT / "companions_pscomppars.csv", index=False)
    print(f"{len(comp)} planets in {comp.hostname.nunique()} multi-planet systems; "
          f"{len(outer)} non-target planets, {int(outer.outer_candidate.sum())} with P > 100 d or no P")

    names = outer.loc[outer.outer_candidate, "pl_name"].tolist()
    ps = []
    for i in range(0, len(names), 100):
        ps.append(tap(f"select {COLS_PS} from ps where pl_name in ({quoted(names[i:i + 100])})"))
    ps = pd.concat(ps, ignore_index=True) if ps else pd.DataFrame()
    ps.to_csv(OUT / "companions_ps_all_solutions.csv", index=False)
    print(f"{len(ps)} published solutions for {ps.pl_name.nunique() if len(ps) else 0} companions")

    from astroquery.vizier import Vizier
    # CDS main host timed out from the MacBook on 2026-11-22; fall back to the CfA mirror
    for server in ("vizier.cds.unistra.fr", "vizier.cfa.harvard.edu"):
        try:
            v = Vizier(columns=["**"], row_limit=-1, vizier_server=server)
            for cat in ("J/ApJ/785/126", "J/ApJ/821/89"):
                tabs = v.get_catalogs(cat)
                print(cat, f"({server}) tables:", list(tabs.keys()))
                for key in tabs.keys():
                    fn = OUT / ("vizier_" + key.replace("/", "_") + ".csv")
                    tabs[key].to_pandas().to_csv(fn, index=False)
                    print("  wrote", fn.name, len(tabs[key]), "rows")
            break
        except Exception as exc:  # report, don't hide
            print(f"VizieR step failed on {server}:", repr(exc))


if __name__ == "__main__":
    main()
