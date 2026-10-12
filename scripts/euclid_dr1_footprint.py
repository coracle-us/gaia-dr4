"""Which target hosts fall in the expected Euclid DR1 footprint.

Input : data/raw/euclid/cgv_map_dr1input_o13.fits.gz (ESA, euclid-dr1-expected-coverage: HEALPix
        nside 8192, NESTED, full sky, column T; ~2,100 deg^2 of DR1 *input* data, of which ~1,900
        deg^2 are expected in the release; "inclusion in DR1 is not guaranteed")
        data/targets.csv (transiting giant-planet hosts) + extra positions (e.g. TESS candidates)
Output: results/euclid_dr1_hosts.csv (host, ra, dec, coverage value T, in_dr1); summary printed:
        footprint area, Galactic-latitude mix of the footprint, hosts inside.
The 6.4 GB map is streamed row by row (1024 pixels per FITS row), never loaded whole.
"""
from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy_healpix import HEALPix
import astropy.units as u

ROOT = Path(__file__).resolve().parents[1]
MAP = ROOT / "data/raw/euclid/cgv_map_dr1input_o13.fits.gz"
NSIDE, PER_ROW = 8192, 1024
ROW_BYTES = PER_ROW * 8
HEADER_BYTES = 2 * 2880          # primary HDU header + BINTABLE header (one block each)
PIX_AREA_DEG2 = 41252.96125 / (12 * NSIDE**2)
EXTRA = {"TIC 159614322 (CTOI)": (220.00977, 65.48404),
         "TIC 147035057 (CTOI)": (104.27440, -16.47848)}


def stream_map(wanted_rows: set, lat_every: int = 4096):
    """Return ({row: values}, n_covered_pixels, sample of covered pixel indices)."""
    rows, n_cov, sample = {}, 0, []
    with gzip.open(MAP, "rb") as f:
        head = f.read(HEADER_BYTES)
        assert b"NSIDE   =                 8192" in head and b"NESTED" in head
        row, n_rows = 0, 12 * NSIDE**2 // PER_ROW           # trailing bytes are FITS padding
        while row < n_rows:
            buf = f.read(ROW_BYTES * min(256, n_rows - row))
            if not buf:
                break
            block = np.frombuffer(buf, dtype=">f8").reshape(-1, PER_ROW)
            nz = block > 0
            n_cov += int(nz.sum())
            for k in np.where(nz.any(axis=1))[0]:
                r = row + k
                if r in wanted_rows:
                    rows[r] = block[k].astype(float)
                if r % lat_every == 0:
                    sample.extend((r * PER_ROW + np.where(nz[k])[0]).tolist())
            for r in wanted_rows & set(range(row, row + len(block))):
                rows.setdefault(r, block[r - row].astype(float))
            row += len(block)
    return rows, n_cov, np.array(sample)


def main() -> None:
    hp = HEALPix(nside=NSIDE, order="nested")
    tg = pd.read_csv(ROOT / "data/targets.csv", dtype={"source_id": "Int64"})
    hosts = tg.drop_duplicates("hostname")[["hostname", "source_id", "ra", "dec", "pl_orbper",
                                            "dist_pc", "phot_g_mean_mag"]].copy()
    extra = pd.DataFrame([dict(hostname=k, ra=v[0], dec=v[1]) for k, v in EXTRA.items()])
    pts = pd.concat([hosts, extra], ignore_index=True)
    pix = hp.lonlat_to_healpix(pts.ra.to_numpy() * u.deg, pts.dec.to_numpy() * u.deg)
    pts["pix"] = pix
    rows, n_cov, sample = stream_map(set((pix // PER_ROW).tolist()))
    pts["coverage_T"] = [rows[p // PER_ROW][p % PER_ROW] if (p // PER_ROW) in rows else 0.0
                         for p in pts.pix]
    pts["in_dr1"] = pts.coverage_T > 0
    out = ROOT / "results/euclid_dr1_hosts.csv"
    pts.drop(columns="pix").to_csv(out, index=False, float_format="%.6g")

    area = n_cov * PIX_AREA_DEG2
    lon, lat = hp.healpix_to_lonlat(sample)
    b = SkyCoord(lon.to(u.deg), lat.to(u.deg), frame="icrs").galactic.b.deg
    print(f"DR1 input footprint: {area:,.0f} deg^2 (map T > 0); |b| of footprint pixels: "
          f"median {np.median(np.abs(b)):.0f} deg, 10-90% {np.percentile(np.abs(b), 10):.0f}-"
          f"{np.percentile(np.abs(b), 90):.0f} deg")
    h = pts[pts.hostname.isin(hosts.hostname)]
    ins = h[h.in_dr1]
    print(f"transiting giant-planet hosts in footprint: {len(ins)} of {len(h)}")
    for name in EXTRA:
        r = pts[pts.hostname == name].iloc[0]
        print(f"{name}: {'IN' if r.in_dr1 else 'not in'} DR1 footprint (T = {r.coverage_T:.3g})")
    if len(ins):
        sens = pd.read_csv(ROOT / "results/target_sensitivity.csv")
        ins = ins.merge(sens[["hostname", "m50_mj_P3"]], on="hostname", how="left")
        cols = ["hostname", "dist_pc", "phot_g_mean_mag", "pl_orbper", "m50_mj_P3"]
        print(ins.sort_values("dist_pc")[cols].round(2).to_string(index=False))


if __name__ == "__main__":
    main()
