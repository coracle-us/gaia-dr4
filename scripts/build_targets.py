"""Build data/targets.csv: transiting giant planets cross-matched to Gaia DR3,
with the expected DR4 astrometric signal of hypothetical outer companions.

Sources
-------
* NASA Exoplanet Archive, table ``pscomppars`` (TAP), transiting planets.
* Gaia DR3 ``gaiadr3.gaia_source`` (G, G_RVS, RUWE, parallax, matched transits)
  and ``external.gaiaedr3_distance`` (Bailer-Jones et al. 2021 geometric distance),
  via TAP upload on the ESA Gaia archive.

Raw query results are cached in data/raw/ (never edited); delete them to re-query.

Sensitivity model (rough, Phase 0 ranking only -- calibrate with injection-recovery)
-----------------------------------------------------------------------------------
* Photocentre semi-major axis of the star for a dark companion:
      alpha [mas] = (M_c / M_star) * a[AU] / d[pc] * 1000
* Per-FoV-transit along-scan precision sigma_fov(G) = sigma_ccd(G) / sqrt(8), with
  sigma_ccd the "EDR3 adjusted" per-CCD curve of Holl et al. 2023 (arXiv:2206.05439,
  Fig. 3) as tabulated in El-Badry's gaiamock (al_uncertainty_per_ccd_interp).
  For G < 13 we add 0.02 mas in quadrature (gaiamock draws U(0, 0.04) mas extra noise).
* N_DR4 = DR3 astrometric_matched_transits * (66 / 34) months, * 0.9 for rejected epochs.
* S/N = alpha * sqrt(N_DR4 / 2) / sigma_fov  (RMS of a circular orbit projected on
  random scan angles is alpha/sqrt(2)). Ignores absorption of the orbit into the
  5-parameter solution, which matters when P approaches or exceeds the 5.5-yr baseline;
  the P_<m>MJ_<a>AU_yr columns let you flag those.
"""
from __future__ import annotations

import io
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from astroquery.gaia import Gaia

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
OUT = ROOT / "data" / "targets.csv"

ARCHIVE_TAP = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync"
MJUP_MSUN = 9.5458e-4          # IAU nominal Jupiter / Sun mass ratio
DR3_MONTHS, DR4_MONTHS = 34.0, 66.0
DR4_BASELINE_YR = DR4_MONTHS / 12.0
COMPANION_MJ = (1, 5)
COMPANION_AU = (1, 3, 5)

# Holl et al. 2023 Fig. 3 "EDR3 adjusted" per-CCD AL uncertainty [mas], via gaiamock
_G_GRID = [4, 5, 6, 7, 8.2, 8.4, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
_SIG_CCD = [0.4, 0.35, 0.15, 0.17, 0.23, 0.13, 0.13, 0.135, 0.125, 0.13, 0.15, 0.23,
            0.36, 0.63, 1.05, 2.05, 4.1]


def sigma_fov_mas(g):
    sig = np.interp(g, _G_GRID, _SIG_CCD) / np.sqrt(8.0)
    return np.where(g < 13, np.hypot(sig, 0.02), sig)


def fetch_archive() -> pd.DataFrame:
    cache = RAW / "nea_pscomppars_transiting.csv"
    if not cache.exists():
        cols = ("pl_name,hostname,gaia_dr3_id,tran_flag,pl_bmassj,pl_bmassprov,pl_radj,"
                "pl_orbper,pl_orbsmax,st_mass,st_teff,sy_dist,sy_gaiamag,sy_snum,sy_pnum,ra,dec")
        q = f"select {cols} from pscomppars where tran_flag = 1"
        r = requests.get(ARCHIVE_TAP, params={"query": q, "format": "csv"}, timeout=300)
        r.raise_for_status()
        RAW.mkdir(parents=True, exist_ok=True)
        cache.write_text(r.text, encoding="utf-8")
    return pd.read_csv(cache)


def select_giants(nea: pd.DataFrame) -> pd.DataFrame:
    """Giant = M_p >= 0.3 M_J, or (no mass and R_p >= 0.7 R_J); drop M_p > 13 M_J."""
    m, r = nea["pl_bmassj"], nea["pl_radj"]
    giant = (m >= 0.3) | (m.isna() & (r >= 0.7))
    return nea[giant & ~(m > 13)].copy()


def fetch_gaia(source_ids: np.ndarray) -> pd.DataFrame:
    cache = RAW / "gaiadr3_targets.csv"
    if not cache.exists():
        # Chunked IN-list queries: TAP uploads were dropped by the archive (Oct 2026,
        # "archive in evolution" for DR4), plain sync queries work.
        ids = [int(s) for s in source_ids]
        parts = []
        for i in range(0, len(ids), 200):
            id_list = ",".join(str(s) for s in ids[i:i + 200])
            q = f"""
            select g.source_id, g.ra as gaia_ra, g.dec as gaia_dec, g.parallax,
                   g.parallax_error, g.phot_g_mean_mag, g.grvs_mag, g.rv_nb_transits, g.ruwe,
                   g.astrometric_matched_transits, g.astrometric_excess_noise,
                   g.astrometric_excess_noise_sig, g.non_single_star, g.ecl_lat,
                   d.r_med_geo, d.r_lo_geo, d.r_hi_geo
            from gaiadr3.gaia_source as g
            left join external.gaiaedr3_distance as d on d.source_id = g.source_id
            where g.source_id in ({id_list})
            """
            parts.append(Gaia.launch_job(q).get_results().to_pandas())
        pd.concat(parts).to_csv(cache, index=False)
    return pd.read_csv(cache)


def main() -> None:
    nea = fetch_archive()
    tg = select_giants(nea)
    tg["source_id"] = (tg["gaia_dr3_id"].str.extract(r"(\d+)\s*$")[0]
                       .astype("Int64"))
    ids = tg["source_id"].dropna().unique()
    gaia = fetch_gaia(ids)
    df = tg.merge(gaia, on="source_id", how="left")

    df["dist_pc"] = df["r_med_geo"].fillna(1000.0 / df["parallax"]).fillna(df["sy_dist"])
    df["mstar_msun"] = df["st_mass"].fillna(1.0)
    df["mstar_assumed"] = df["st_mass"].isna()
    g = df["phot_g_mean_mag"].fillna(df["sy_gaiamag"])
    df["sigma_fov_mas"] = sigma_fov_mas(g.to_numpy())
    df["n_fov_dr4_est"] = np.round(df["astrometric_matched_transits"] * DR4_MONTHS / DR3_MONTHS * 0.9)
    df["epoch_rv_expected"] = df["grvs_mag"] <= 14

    for mc in COMPANION_MJ:
        for a in COMPANION_AU:
            tag = f"{mc}MJ_{a}AU"
            alpha_mas = mc * MJUP_MSUN / df["mstar_msun"] * a / df["dist_pc"] * 1e3
            df[f"alpha_uas_{tag}"] = (alpha_mas * 1e3).round(2)
            df[f"snr_{tag}"] = (alpha_mas * np.sqrt(df["n_fov_dr4_est"] / 2)
                                / df["sigma_fov_mas"]).round(2)
            df[f"P_yr_{tag}"] = np.sqrt(a**3 / (df["mstar_msun"] + mc * MJUP_MSUN)).round(2)

    keep = ["pl_name", "hostname", "source_id", "ra", "dec", "pl_bmassj", "pl_bmassprov",
            "pl_radj", "pl_orbper", "pl_orbsmax", "st_mass", "mstar_msun", "mstar_assumed",
            "st_teff", "sy_snum", "sy_pnum", "dist_pc", "r_lo_geo", "r_hi_geo", "parallax",
            "parallax_error", "phot_g_mean_mag", "grvs_mag", "rv_nb_transits",
            "epoch_rv_expected", "ruwe", "astrometric_excess_noise",
            "astrometric_excess_noise_sig", "non_single_star", "astrometric_matched_transits",
            "n_fov_dr4_est", "sigma_fov_mas"]
    keep += [c for c in df.columns if c.startswith(("alpha_uas_", "snr_", "P_yr_"))]
    df = df[keep].sort_values("snr_5MJ_3AU", ascending=False)
    df.to_csv(OUT, index=False, float_format="%.6g")

    n = len(df)
    print(f"transiting giants: {n} planets, {df['hostname'].nunique()} hosts; "
          f"with DR3 match: {df['phot_g_mean_mag'].notna().sum()}")
    hosts = df.drop_duplicates("hostname")
    print(f"hosts with G_RVS<=14 (epoch RVs expected): {hosts['epoch_rv_expected'].sum()}")
    for tag in ("1MJ_1AU", "1MJ_3AU", "5MJ_1AU", "5MJ_3AU", "5MJ_5AU"):
        s = hosts[f"snr_{tag}"]
        print(f"  {tag:8s} hosts with S/N>=5: {(s >= 5).sum():4d}  >=10: {(s >= 10).sum():4d}"
              f"  >=20: {(s >= 20).sum():4d}")


if __name__ == "__main__":
    main()
