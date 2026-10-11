"""Impostor (luminous twin binary) analysis for the target hosts.

Inputs : data/raw/gaiaforecasts/DR4_mock_{exoplanet,planet_impostor}_catalog.csv
           (Lammers & Winn 2025 mocks, github.com/CalebLammers/GaiaForecasts, MIT)
         data/raw/mamajek/EEM_dwarf_UBVIJHK_colors_Teff_2022.04.16.txt
         results/target_sensitivity.csv, data/targets.csv
         queried + cached in data/raw/impostor/: NASA Exoplanet Archive pl_rvamp (literature RV
           semi-amplitude => RV follow-up exists) and Gaia DR3 RV statistics for the hosts.
Outputs: results/impostor_mock_purity.csv, results/impostor_vetting.csv; summary printed.

(a) Mock population: purity (planets / all detections) for FGK-like hosts in our range.
(b) Per host with m50(P = 3 yr) <= 13 M_J: the main-sequence twin that mimics the m50 signal at
    P = 2, 3, 4 yr (companions.impostor), its overluminosity, line splitting and blended RV.
(c) Evidence available now: literature RVs (high-resolution spectra would show the 10-30 km/s
    line splitting), DR3 RV variability (rv_chisq_pvalue < 0.01 & rv_renormalised_gof > 4 &
    rv_nb_transits >= 10; DR3 documentation), DR3 per-epoch RV scatter
    sigma_epoch = sqrt((err^2 - 0.11^2) 2N / pi), and overluminosity dM_G = M_G(obs) - M_G,MS(Teff).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import requests
from astroquery.gaia import Gaia

from companions.impostor import MainSequence, twin_mimicking_signal
from companions.posterior import MJUP_MSUN

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/raw/impostor"
GF = ROOT / "data/raw/gaiaforecasts"
RES = ROOT / "results"
SPECTRO_RESOLUTION_KMS = 6.0     # line splitting resolvable in R ~ 50,000 RV-follow-up spectra


def mock_purity() -> pd.DataFrame:
    pl = pd.read_csv(GF / "DR4_mock_exoplanet_catalog.csv")
    im = pd.read_csv(GF / "DR4_mock_planet_impostor_catalog.csv")
    sel_pl = lambda d: d[(d["Stellar mass [M_\\odot]"].between(0.7, 1.5))
                         & (d["G-band mag"] < 13.5) & (d["True distance [pc]"] < 300)]
    sel_im = lambda d: d[(d["Apparent stellar mass [M_\\odot]"].between(0.7, 1.5))
                         & (d["Apparent G-band mag"] < 13.5) & (d["True distance [pc]"] < 300)]
    rows = []
    for plo, phi in [(100, 365), (365, 730), (730, 1100), (1100, 1500)]:
        p = sel_pl(pl)
        p = p[p["Best-fit period [days]"].between(plo, phi)]
        i = sel_im(im)
        i = i[i["Best-fit period [days]"].between(plo, phi)]
        q = i["Secondary stellar mass [M_\\odot]"] / i["Primary stellar mass [M_\\odot]"]
        rows.append(dict(period_days=f"{plo}-{phi}", n_planets=len(p), n_impostors=len(i),
                         purity=len(p) / max(len(p) + len(i), 1),
                         impostor_q_median=float(q.median()) if len(i) else np.nan,
                         impostor_q_min=float(q.min()) if len(i) else np.nan))
    return pd.DataFrame(rows)


def query_cached(name, fetch) -> pd.DataFrame:
    path = RAW / name
    if not path.exists():
        RAW.mkdir(parents=True, exist_ok=True)
        fetch().to_csv(path, index=False)
    return pd.read_csv(path)


def fetch_rvamp() -> pd.DataFrame:
    q = "select pl_name, hostname, pl_rvamp, pl_bmassprov from pscomppars where tran_flag = 1"
    r = requests.get("https://exoplanetarchive.ipac.caltech.edu/TAP/sync",
                     params={"query": q, "format": "csv"}, timeout=300)
    r.raise_for_status()
    from io import StringIO
    return pd.read_csv(StringIO(r.text))


def fetch_dr3_rv(ids) -> pd.DataFrame:
    parts = []
    for i in range(0, len(ids), 200):
        id_list = ",".join(str(int(s)) for s in ids[i:i + 200])
        parts.append(Gaia.launch_job(f"""
            select source_id, radial_velocity_error, rv_nb_transits, rv_method_used,
                   rv_chisq_pvalue, rv_renormalised_gof, rv_amplitude_robust, grvs_mag
            from gaiadr3.gaia_source where source_id in ({id_list})""").get_results().to_pandas())
    return pd.concat(parts)


def main() -> None:
    purity = mock_purity()
    purity.to_csv(RES / "impostor_mock_purity.csv", index=False, float_format="%.4g")
    print("(a) GaiaForecasts DR4 mocks, hosts 0.7-1.5 M_sun, G < 13.5, d < 300 pc:")
    print(purity.round(3).to_string(index=False))

    ms = MainSequence.load()
    sens = pd.read_csv(RES / "target_sensitivity.csv", dtype={"source_id": "Int64"})
    tg = pd.read_csv(ROOT / "data/targets.csv", dtype={"source_id": "Int64"}).drop_duplicates("hostname")
    rvamp = query_cached("nea_rvamp.csv", fetch_rvamp)
    has_rv = rvamp.groupby("hostname")["pl_rvamp"].apply(lambda s: s.notna().any())
    dr3 = query_cached("dr3_rv_stats.csv", lambda: fetch_dr3_rv(sens.source_id.dropna().unique()))
    dr3["source_id"] = dr3["source_id"].astype("Int64")

    df = (sens[sens.m50_mj_P3 <= 13]
          .merge(tg[["hostname", "st_teff", "phot_g_mean_mag", "parallax"]], on="hostname", how="left")
          .merge(dr3, on="source_id", how="left"))
    df["lit_rv_followup"] = df.hostname.map(has_rv).fillna(False).astype(bool)
    teff_tab = np.array([(t, g) for t, g in zip(*_teff_mg_table(ms))])
    df["dM_G"] = (df.phot_g_mean_mag + 5 * np.log10(df.parallax / 100.0)
                  - np.interp(df.st_teff, teff_tab[:, 0], teff_tab[:, 1]))
    df["dr3_rv_variable"] = ((df.rv_chisq_pvalue < 0.01) & (df.rv_renormalised_gof > 4)
                             & (df.rv_nb_transits >= 10))
    n = df.rv_nb_transits.astype(float)
    err = df.radial_velocity_error.astype(float)
    df["dr3_sigma_epoch_kms"] = np.where(df.rv_method_used == 1, np.sqrt(
        np.clip(err**2 - 0.11**2, 0, None) * 2 * n / np.pi), np.nan)

    for p in (2.0, 3.0, 4.0):
        out = []
        for _, h in df.iterrows():
            plx, m1 = 1000.0 / h.dist_pc, h.mstar
            mc = h[f"m50_mj_P{p:g}"] * MJUP_MSUN
            a0 = ((m1 + mc) * p**2) ** (1 / 3) * mc / (m1 + mc) * plx
            out.append(twin_mimicking_signal(ms, m1, p, plx, a0) or {})
        tw = pd.DataFrame(out)
        df[f"twin_q_P{p:g}"] = tw["q"].to_numpy()
        df[f"twin_overlum_P{p:g}"] = tw["overluminosity_g"].to_numpy()
        df[f"twin_split_kms_P{p:g}"] = tw["max_line_separation"].to_numpy()
        df[f"twin_blended_rv_kms_P{p:g}"] = tw["blended_rv_amp_rp"].to_numpy()
        # fraction of an edge-on circular orbit with |dv| > resolution
        frac = 1 - 2 / np.pi * np.arcsin(np.clip(SPECTRO_RESOLUTION_KMS / tw["max_line_separation"], 0, 1))
        df[f"twin_frac_resolved_P{p:g}"] = frac.to_numpy()

    df["twin_excluded_by_lit_rv"] = df.lit_rv_followup & (df.twin_frac_resolved_P3 > 0.5)
    df.to_csv(RES / "impostor_vetting.csv", index=False, float_format="%.4g")

    print(f"\n(b,c) {len(df)} hosts with m50(P=3 yr) <= 13 M_J")
    print(f"  twin q (P=3 yr): median {df.twin_q_P3.median():.4f}, range {df.twin_q_P3.min():.4f}-"
          f"{df.twin_q_P3.max():.4f}; overluminosity median {df.twin_overlum_P3.median():.2f} mag")
    print(f"  line splitting (edge-on, P=2/3/4 yr) median "
          f"{df.twin_split_kms_P2.median():.1f}/{df.twin_split_kms_P3.median():.1f}/"
          f"{df.twin_split_kms_P4.median():.1f} km/s; blended RV amplitude median "
          f"{df.twin_blended_rv_kms_P3.median() * 1000:.0f} m/s (P=3 yr)")
    print(f"  literature RV follow-up: {df.lit_rv_followup.sum()} of {len(df)}; "
          f"twin excluded by it (resolved > 50% of orbit): {df.twin_excluded_by_lit_rv.sum()}")
    print(f"  DR3 RV-variable (chisq & gof & N>=10): {int(df.dr3_rv_variable.sum())} "
          f"(of {df.rv_chisq_pvalue.notna().sum()} with statistics); median DR3 per-epoch "
          f"scatter {np.nanmedian(df.dr3_sigma_epoch_kms):.2f} km/s")
    print(f"  overluminous dM_G < -0.5 mag (twin-like or evolved): {(df.dM_G < -0.5).sum()}; "
          f"|dM_G| < 0.3: {(df.dM_G.abs() < 0.3).sum()}")
    no_rv = df[~df.lit_rv_followup].sort_values("m50_mj_P3")
    print(f"  hosts WITHOUT literature RVs (twin not yet excluded): {len(no_rv)}")
    print(no_rv[["hostname", "dist_pc", "g_mag", "m50_mj_P3", "dM_G", "dr3_rv_variable",
                 "dr3_sigma_epoch_kms", "twin_split_kms_P3"]].head(15).round(2).to_string(index=False))


def _teff_mg_table(ms_unused):
    """Main-sequence Teff -> M_G from the same Pecaut & Mamajek table (sorted by Teff)."""
    from companions.impostor import TABLE
    lines = TABLE.read_text(encoding="utf-8").splitlines()
    header = next(l.split() for l in lines if l.startswith("#SpT"))
    it, ig = header.index("Teff"), header.index("M_G")
    rows = []
    for line in lines:
        parts = line.split()
        if line.startswith("#") or len(parts) != len(header) or not parts[0].endswith("V"):
            continue
        try:
            rows.append((float(parts[it]), float(parts[ig])))
        except ValueError:
            continue
    arr = np.array(sorted(set(rows)))
    return arr[:, 0], arr[:, 1]


if __name__ == "__main__":
    main()
