"""Euclid Q1 footprint (MER tiles) and the Gaia DR3 stars within 200 pc inside it.

Queries (cached in data/raw/euclid/):
  * ESA Euclid Science Archive (astroquery.esa.euclid, public 'PDR' environment):
    distinct VIS MER tiles in q1.mosaic_product (tile_index, ra, dec).
  * Gaia DR3: stars with parallax > 5 mas and parallax_over_error > 5 in circles around each
    Q1 field; kept if inside a tile core (MER tiles: 30' x 30' cores; assumption, see notes).
  * Coverage: tile cores over-count (edge tiles are partly empty), so 120 random hosts are checked
    for >= 1 MER source within 30" (sync cone queries; TAP uploads need a login), giving the
    covered fraction (cached: data/raw/euclid/q1_hosts_coverage_sample.csv).
Output: results/euclid_q1_hosts.csv; prints Q1 area and host counts, and scales to DR1 area.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/raw/euclid"
TILE_CORE_DEG = 0.5            # MER tile core side (30 arcmin)
DR1_AREA_DEG2 = 1900.0         # announced DR1 area (to be verified, notes/euclid_scoping.md)


def tiles() -> pd.DataFrame:
    path = RAW / "q1_vis_tiles.csv"
    if not path.exists():
        from astroquery.esa.euclid import Euclid
        RAW.mkdir(parents=True, exist_ok=True)
        r = Euclid.launch_job("select distinct tile_index, ra, dec from q1.mosaic_product "
                              "where filter_name = 'VIS'").get_results().to_pandas()
        r.to_csv(path, index=False)
    return pd.read_csv(path).drop_duplicates("tile_index")


def field_groups(t: pd.DataFrame, link_deg=3.0):
    """Group tiles into contiguous fields (single-linkage on tile centres)."""
    lab = -np.ones(len(t), int)
    xyz = np.c_[np.cos(np.radians(t.dec)) * np.cos(np.radians(t.ra)),
                np.cos(np.radians(t.dec)) * np.sin(np.radians(t.ra)), np.sin(np.radians(t.dec))]
    k = 0
    for i in range(len(t)):
        if lab[i] >= 0:
            continue
        stack, lab[i] = [i], k
        while stack:
            j = stack.pop()
            near = np.where((lab < 0) & (np.degrees(np.arccos(np.clip(xyz @ xyz[j], -1, 1))) < link_deg))[0]
            lab[near] = k
            stack.extend(near.tolist())
        k += 1
    return lab


def gaia_hosts(t: pd.DataFrame) -> pd.DataFrame:
    path = RAW / "q1_gaia_dr3_200pc.csv"
    if not path.exists():
        from astroquery.gaia import Gaia
        parts = []
        for g, sub in t.groupby("field"):
            ra0, dec0 = sub.ra.mean(), sub.dec.mean()
            sep = np.degrees(np.arccos(np.clip(
                np.sin(np.radians(dec0)) * np.sin(np.radians(sub.dec))
                + np.cos(np.radians(dec0)) * np.cos(np.radians(sub.dec)) * np.cos(np.radians(sub.ra - ra0)), -1, 1)))
            rad = sep.max() + 0.6
            q = f"""select source_id, ra, dec, parallax, parallax_error, pmra, pmdec, phot_g_mean_mag,
                           bp_rp, ruwe from gaiadr3.gaia_source
                    where 1 = contains(point('ICRS', ra, dec), circle('ICRS', {ra0}, {dec0}, {rad}))
                      and parallax > 5 and parallax_over_error > 5"""
            Gaia.ROW_LIMIT = -1                    # sync job (async queue is slow under load)
            parts.append(Gaia.launch_job(q).get_results().to_pandas().assign(field=g))
        pd.concat(parts).to_csv(path, index=False)
    return pd.read_csv(path)


def in_tile_core(stars: pd.DataFrame, t: pd.DataFrame) -> np.ndarray:
    keep = np.zeros(len(stars), bool)
    for _, tl in t.iterrows():
        dra = (stars.ra - tl.ra + 180) % 360 - 180
        keep |= ((np.abs(dra * np.cos(np.radians(tl.dec))) < TILE_CORE_DEG / 2)
                 & (np.abs(stars.dec - tl.dec) < TILE_CORE_DEG / 2))
    return keep


def coverage_fraction(stars: pd.DataFrame, n=120, seed=1) -> float:
    path = RAW / "q1_hosts_coverage_sample.csv"
    if not path.exists():
        from astroquery.esa.euclid import Euclid
        rows = []
        for r in stars.sample(n, random_state=seed).itertuples():
            q = (f"select count(*) as n from catalogue.mer_catalogue where 1=contains(point('ICRS', "
                 f"right_ascension, declination), circle('ICRS', {r.ra}, {r.dec}, 0.0083333))")
            try:
                rows.append((r.source_id, r.field, int(Euclid.launch_job(q).get_results()["n"][0])))
            except Exception:
                rows.append((r.source_id, r.field, -1))
        pd.DataFrame(rows, columns=["source_id", "field", "n_mer"]).to_csv(path, index=False)
    d = pd.read_csv(path)
    d = d[d.n_mer >= 0]
    return float((d.n_mer > 0).mean())


def main() -> None:
    t = tiles()
    t["field"] = field_groups(t)
    area = len(t) * TILE_CORE_DEG**2
    stars = gaia_hosts(t)
    stars = stars[in_tile_core(stars, t)].drop_duplicates("source_id")
    stars["dist_pc"] = 1000.0 / stars.parallax
    stars.to_csv(ROOT / "results/euclid_q1_hosts.csv", index=False, float_format="%.6g")
    print(f"Q1 VIS tiles: {len(t)} in {t.field.nunique()} fields -> area ~ {area:.1f} deg^2 "
          f"(tile cores {TILE_CORE_DEG * 60:.0f}')")
    for g, sub in t.groupby("field"):
        n = (stars.field == g).sum()
        print(f"  field {g}: {len(sub)} tiles, centre ({sub.ra.mean():.2f}, {sub.dec.mean():.2f}), "
              f"{n} Gaia stars < 200 pc")
    f = coverage_fraction(stars)
    eff_area = f * area
    print(f"covered fraction (MER source within 30\"): {f:.2f} -> effective area ~ {eff_area:.0f} deg^2")
    for d in (50, 100, 200):
        n_d = (stars.dist_pc < d).sum()
        print(f"Gaia DR3 stars within {d:3d} pc: {n_d} in tile cores, ~{f * n_d:.0f} covered; "
              f"density {f * n_d / eff_area:.2f}/deg^2 -> DR1 ({DR1_AREA_DEG2:.0f} deg^2) ~ "
              f"{f * n_d / eff_area * DR1_AREA_DEG2:,.0f}")
    print("(Q1 fields are high-latitude; DR1 wide patches have a broader latitude mix: approximate)")


if __name__ == "__main__":
    main()
