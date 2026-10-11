"""Cross-match target hosts to pre-DR4 astrometric companion evidence.

Inputs : data/targets.csv (scripts/build_targets.py)
Queries (cached as raw results in data/raw/prior_astrometry/, never edited):
  * Gaia DR3 NSS: gaiadr3.nss_two_body_orbit, nss_acceleration_astro, nss_non_linear_spectro
  * HGCA EDR3 (Brandt 2021, ApJS 254, 42; VizieR J/ApJS/254/42/catalog)
  * Kervella et al. 2022 proper-motion anomaly (A&A 657, A7; VizieR J/A+A/657/A7/tablea1)
Output : data/target_prior_astrometry.csv, one row per host.

HGCA chi2 is for a constant proper motion (2 d.o.f.): chi2 > 11.8 is 3 sigma, > 28.7 is 5 sigma.
Kervella snrPMaH2EG3b is the S/N of the Gaia EDR3 vs long-term (Hip-EDR3) PM anomaly;
BinH2EG3b = 1 for S/N > 3. M25au is their companion mass for a 5 AU orbit (M_Jup).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pyvo
from astroquery.gaia import Gaia

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "prior_astrometry"
OUT = ROOT / "data" / "target_prior_astrometry.csv"
VIZIER_TAP = "https://tapvizier.cds.unistra.fr/TAPVizieR/tap"
CHI2_3SIG, CHI2_5SIG = 11.8, 28.7


def cached(name: str, fetch) -> pd.DataFrame:
    path = RAW / name
    if not path.exists():
        RAW.mkdir(parents=True, exist_ok=True)
        fetch().to_csv(path, index=False)
    return pd.read_csv(path)


def in_chunks(run, template: str, ids, size: int = 200) -> pd.DataFrame:
    parts = []
    for i in range(0, len(ids), size):
        id_list = ",".join(str(int(s)) for s in ids[i:i + size])
        parts.append(run(template.format(ids=id_list)))
    return pd.concat(parts, ignore_index=True)


def gaia(q: str) -> pd.DataFrame:
    return Gaia.launch_job(q).get_results().to_pandas()


def vizier(q: str) -> pd.DataFrame:
    return pyvo.dal.TAPService(VIZIER_TAP).search(q).to_table().to_pandas()


def main() -> None:
    # source_id must be read as nullable int: with NaNs present pandas would parse it as
    # float64, which silently corrupts 19-digit Gaia IDs.
    hosts = (pd.read_csv(ROOT / "data" / "targets.csv", dtype={"source_id": "Int64"})
             .dropna(subset=["source_id"]).drop_duplicates("hostname"))
    ids = hosts["source_id"].astype("int64").unique()

    orb = cached("dr3_nss_two_body_orbit.csv", lambda: in_chunks(gaia, """
        select source_id, nss_solution_type, period, period_error, eccentricity,
               significance, goodness_of_fit, a_thiele_innes, b_thiele_innes,
               f_thiele_innes, g_thiele_innes, semi_amplitude_primary
        from gaiadr3.nss_two_body_orbit where source_id in ({ids})""", ids))
    acc = cached("dr3_nss_acceleration_astro.csv", lambda: in_chunks(gaia, """
        select source_id, nss_solution_type, significance, accel_ra, accel_dec
        from gaiadr3.nss_acceleration_astro where source_id in ({ids})""", ids))
    nls = cached("dr3_nss_non_linear_spectro.csv", lambda: in_chunks(gaia, """
        select source_id, nss_solution_type, mean_velocity, first_deriv_velocity
        from gaiadr3.nss_non_linear_spectro where source_id in ({ids})""", ids))
    hgca = cached("hgca_edr3.csv", lambda: in_chunks(vizier, """
        select "HIP", "Gaia", "chi2" from "J/ApJS/254/42/catalog"
        where "Gaia" in ({ids})""", ids))
    kerv = cached("kervella2022_pma.csv", lambda: in_chunks(vizier, """
        select "HIP", "GaiaEDR3", "snrPMaH2EG3b", "BinH2EG3b", "dVt", "e_dVt",
               "M23au", "M25au", "M210au", "B"
        from "J/A+A/657/A7/tablea1" where "GaiaEDR3" in ({ids})""", ids))

    df = hosts[["hostname", "source_id", "pl_name", "pl_orbper", "pl_bmassj", "dist_pc",
                "phot_g_mean_mag", "ruwe", "non_single_star", "sy_snum", "sy_pnum",
                "snr_5MJ_3AU", "epoch_rv_expected"]].copy()
    df["source_id"] = df["source_id"].astype("int64")

    def agg(t, cols):
        t = t.copy()
        t["source_id"] = t["source_id"].astype("int64")
        return t.groupby("source_id").agg(cols)

    if len(orb):
        df = df.merge(agg(orb, {"nss_solution_type": lambda s: ";".join(sorted(s)),
                                "period": "first", "eccentricity": "first",
                                "significance": "max"})
                      .add_prefix("dr3_orbit_"), left_on="source_id", right_index=True, how="left")
    if len(acc):
        df = df.merge(agg(acc, {"nss_solution_type": "first", "significance": "max"})
                      .add_prefix("dr3_accel_"), left_on="source_id", right_index=True, how="left")
    if len(nls):
        df = df.merge(agg(nls, {"nss_solution_type": "first"}).add_prefix("dr3_rvtrend_"),
                      left_on="source_id", right_index=True, how="left")
    hgca = hgca.rename(columns={"Gaia": "source_id", "HIP": "hip", "chi2": "hgca_chi2"})
    df = df.merge(hgca[["source_id", "hip", "hgca_chi2"]], on="source_id", how="left")
    kerv = kerv.rename(columns={"GaiaEDR3": "source_id", "snrPMaH2EG3b": "pma_snr",
                                "BinH2EG3b": "pma_bin", "dVt": "pma_dvt_ms",
                                "e_dVt": "pma_dvt_err_ms", "M23au": "pma_m2_3au_mj",
                                "M25au": "pma_m2_5au_mj", "M210au": "pma_m2_10au_mj",
                                "B": "kervella_n_bound_cpm"}).drop(columns=["HIP"])
    df = df.merge(kerv, on="source_id", how="left")

    df["hgca_3sigma"] = df["hgca_chi2"] > CHI2_3SIG
    df["hgca_5sigma"] = df["hgca_chi2"] > CHI2_5SIG
    for c in ("dr3_orbit_nss_solution_type", "dr3_accel_nss_solution_type",
              "dr3_rvtrend_nss_solution_type"):
        if c not in df.columns:
            df[c] = np.nan
    flags = {
        "ruwe>1.4": df["ruwe"] > 1.4,
        "dr3_nss_orbit": df["dr3_orbit_nss_solution_type"].notna(),
        "dr3_nss_accel": df["dr3_accel_nss_solution_type"].notna(),
        "dr3_rv_trend": df["dr3_rvtrend_nss_solution_type"].notna(),
        "hgca_3sigma": df["hgca_chi2"] > CHI2_3SIG,
        "pma_snr>3": df["pma_snr"] > 3,
    }
    fl = pd.DataFrame(flags)
    df["prior_evidence"] = fl.apply(lambda r: ";".join(k for k, v in r.items() if v), axis=1)
    df["n_prior_evidence"] = fl.sum(axis=1)
    df = df.sort_values(["n_prior_evidence", "snr_5MJ_3AU"], ascending=[False, False])
    df.to_csv(OUT, index=False, float_format="%.6g")

    print(f"hosts: {len(df)}; in HGCA: {df['hgca_chi2'].notna().sum()}; "
          f"in Kervella+22: {df['pma_snr'].notna().sum()}")
    for k, v in flags.items():
        print(f"  {k:14s} {int(v.sum()):4d}")
    print(f"  any evidence   {(df['prior_evidence'] != '').sum():4d}")
    show = df[df["prior_evidence"] != ""]
    cols = ["hostname", "dist_pc", "phot_g_mean_mag", "ruwe", "hgca_chi2", "pma_snr",
            "pma_m2_5au_mj", "prior_evidence"]
    cols += [c for c in ("dr3_orbit_nss_solution_type", "dr3_orbit_period",
                         "dr3_accel_nss_solution_type") if c in df.columns]
    with pd.option_context("display.width", 250, "display.max_rows", 200):
        print(show[cols].round(2).to_string(index=False))


if __name__ == "__main__":
    main()
