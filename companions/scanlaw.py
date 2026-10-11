"""Predicted Gaia transits (times, scan angles) and along-scan parallax factors.

Times and scan angles: ESA GOST nominal scanning-law predictions on a HEALPix level-6 (nside 64,
nested) grid, as packaged in gaiascanlaw 0.2.0 (Z. Penoyre; PyPI). Only the data tables are used:
    pip download --no-deps gaiascanlaw==0.2.0 -d data/raw/gaiascanlaw   # 89 MB wheel
    then extract gaiascanlaw/data/*.fits from the wheel (zip) into that directory.
Known DR3 astrometric data gaps (cosmos.esa.int/web/gaia/dr3-data-gaps) are removed, as is the
Ecliptic Pole Scanning period (no astrometry). Gaps after DR3 are unknown until DR4 publishes
`gaia_gap`; the DR4 `commanded_scan_law` table will replace this module's source on release.
The prediction is for the pixel centre (~0.9 deg pixels), so a few transits near the field edges
differ from those of the exact position: the pixel (~0.9 deg) is comparable to the 0.7 deg
across-scan FoV width, so only ~40-70% of real transits are reproduced epoch-by-epoch, while the
number of transits, timing and scan angles are statistically representative (tests/test_scanlaw.py).
Good enough for injection-recovery statistics; not for fitting real data.

Parallax factor: computed for Gaia near L2 (0.01 AU beyond the Earth on the Sun-Earth line):
    f_AL = sin(psi) (X sin a - Y cos a) + cos(psi) (X cos a sin d + Y sin a sin d - Z cos d)
with (X, Y, Z) Gaia's barycentric ICRS position in AU.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from astropy.coordinates import get_body_barycentric
from astropy.table import Table, vstack
from astropy.time import Time
from astropy_healpix import HEALPix
import astropy.units as u

from companions.epoch import DAYS_PER_YEAR, JD_REF_EPOCH

ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data/raw/gaiascanlaw/gaiascanlaw/data"
CACHE = ROOT / "data/cache/gaiascanlaw_t_psi.npz"
L2_OFFSET_AU = 0.01


def _t_from_time(time: Time) -> np.ndarray:
    return (time.tcb.jd - JD_REF_EPOCH) / DAYS_PER_YEAR


T_EPSL_END = float(_t_from_time(Time("2014-08-22T21:01:25.599970336", scale="tcb")))
T_DR4_END = float(_t_from_time(Time("2020-01-20T22:01:30.250520158", scale="tcb")))


def parallax_factor_al(t, psi, ra_deg, dec_deg) -> np.ndarray:
    time = Time(JD_REF_EPOCH + np.asarray(t) * DAYS_PER_YEAR, format="jd", scale="tcb")
    earth = get_body_barycentric("earth", time).xyz.to_value(u.au)
    sun = get_body_barycentric("sun", time).xyz.to_value(u.au)
    d = earth - sun
    x, y, z = earth + L2_OFFSET_AU * d / np.linalg.norm(d, axis=0)
    a, dd = np.deg2rad(ra_deg), np.deg2rad(dec_deg)
    return (np.sin(psi) * (x * np.sin(a) - y * np.cos(a))
            + np.cos(psi) * (x * np.cos(a) * np.sin(dd) + y * np.sin(a) * np.sin(dd) - z * np.cos(dd)))


@dataclass
class ScanLaw:
    pix: np.ndarray        # sorted HEALPix (nside 64, nested)
    t: np.ndarray          # yr from J2017.5 (TCB)
    psi: np.ndarray        # rad
    gap_start: np.ndarray  # astrometry gaps, yr from J2017.5
    gap_end: np.ndarray

    @classmethod
    def load(cls, data_dir: Path = DATA_DIR) -> "ScanLaw":
        files = sorted(Path(data_dir).glob("gost_simple_*.fits"))
        if not files:
            raise FileNotFoundError(f"no gost_simple_*.fits in {data_dir}")
        gaps = Table.read(Path(data_dir) / "data_gaps.fits")
        gaps = gaps[gaps["type"] == "astrometry"]
        # gaiascanlaw "decimal year" values are Julian years (TCB): dy = 2017.5 at J2017.5.
        # (Checked against the DR4 pre-release epochs: astropy's calendar decimalyear is off
        # by hours; Julian years match to ~10-40 s.)
        conv = lambda dy: np.asarray(dy, float) - 2017.5
        if CACHE.exists():
            c = np.load(CACHE)
            pix, t, psi = c["pix"], c["t"], c["psi"]
        else:
            tab = vstack([Table.read(f) for f in files])
            order = np.argsort(np.asarray(tab["healpixel"]), kind="stable")
            pix = np.asarray(tab["healpixel"])[order]
            t = conv(tab["scan time [decimal year]"])[order]
            psi = np.asarray(tab["scan angle [radian]"], float)[order]
            CACHE.parent.mkdir(parents=True, exist_ok=True)
            np.savez(CACHE, pix=pix, t=t, psi=psi)
        return cls(pix, t, psi, conv(gaps["start [decimal year]"]), conv(gaps["end [decimal year]"]))

    def transits(self, ra_deg, dec_deg, t_end: float = T_DR4_END):
        """Predicted astrometric FoV transits in DR4: (t [yr from J2017.5], psi [rad]), sorted."""
        hp = HEALPix(nside=64, order="nested")
        p = int(hp.lonlat_to_healpix(ra_deg * u.deg, dec_deg * u.deg))
        lo, hi = np.searchsorted(self.pix, [p, p + 1])
        t, psi = self.t[lo:hi], self.psi[lo:hi]
        keep = (t > T_EPSL_END) & (t < t_end)
        for gs, ge in zip(self.gap_start, self.gap_end):
            keep &= ~((t > gs) & (t < ge))
        order = np.argsort(t[keep])
        return t[keep][order], psi[keep][order]
