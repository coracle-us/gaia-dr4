"""Luminous-binary impostors: twin binaries whose photocentre orbit mimics a planet.

For a binary with mass ratio q = M2/M1 and band flux ratio f = F2/F1, the photocentre moves on
an orbit of size a0 = a_rel |B - beta|, with B = q/(1+q) and beta = f/(1+f). A main-sequence
companion has f < q, so a0 > 0 and shrinks to 0 as q -> 1: near-twins mimic small planet-like
orbits. The flux-weighted ("blended") RV of an unresolved double-lined spectrum moves with the
same factor, v_blend = v_rel (B - beta_band), so twins hide in RV amplitude too. They show up as
(i) line splitting up to the full relative velocity K_rel (double-lined / broadened spectra),
(ii) overluminosity of up to 2.5 log10(2) = 0.75 mag.
Main-sequence M_G and G-RP versus mass: Pecaut & Mamajek (2013), table version 2022.04.16
(data/raw/mamajek/). RP is used as a proxy for the RVS band (845-872 nm).
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]
TABLE = ROOT / "data/raw/mamajek/EEM_dwarf_UBVIJHK_colors_Teff_2022.04.16.txt"
AU_PER_YR_KMS = 4.740470463


def photocentre_factor(q, flux_ratio):
    return q / (1 + q) - flux_ratio / (1 + flux_ratio)


@dataclass
class MainSequence:
    mass: np.ndarray
    m_g: np.ndarray
    g_rp: np.ndarray

    @classmethod
    def load(cls, path: Path = TABLE) -> "MainSequence":
        if not Path(path).exists():
            raise FileNotFoundError(path)
        lines = Path(path).read_text(encoding="utf-8").splitlines()
        header = next(l.split() for l in lines if l.startswith("#SpT"))
        i_m, i_g, i_grp = header.index("Msun"), header.index("M_G"), header.index("G-Rp")
        rows = []
        for line in lines:
            parts = line.split()
            if line.startswith("#") or len(parts) != len(header) or not parts[0].endswith("V"):
                continue
            try:
                rows.append((float(parts[i_m]), float(parts[i_g]), float(parts[i_grp])))
            except ValueError:                       # '...' = not tabulated
                continue
        arr = np.array(sorted(set(rows)))
        return cls(arr[:, 0], arr[:, 1], arr[:, 2])

    def abs_g(self, mass):
        return np.interp(mass, self.mass, self.m_g)

    def abs_rp(self, mass):
        return np.interp(mass, self.mass, self.m_g - self.g_rp)

    def flux_ratio(self, m1, q, band="G"):
        f = self.abs_g if band == "G" else self.abs_rp
        return 10 ** (-0.4 * (f(q * m1) - f(m1)))


def rv_signatures(m1, m2, period_yr, inc=np.pi / 2, e=0.0) -> dict:
    a_rel = ((m1 + m2) * period_yr**2) ** (1 / 3)
    k_rel = 2 * np.pi * a_rel * np.sin(inc) / (period_yr * np.sqrt(1 - e**2)) * AU_PER_YR_KMS
    return dict(a_rel_au=a_rel, k_rel=k_rel, k1=k_rel * m2 / (m1 + m2),
                k2=k_rel * m1 / (m1 + m2), max_line_separation=k_rel * (1 + e))


def twin_mimicking_signal(ms: MainSequence, m1, period_yr, parallax_mas, a0_target_mas,
                          inc=np.pi / 2) -> dict | None:
    """Main-sequence companion (q < 1) whose G-band photocentre orbit has size a0_target."""
    def a0(q):
        a_rel_mas = ((m1 * (1 + q)) * period_yr**2) ** (1 / 3) * parallax_mas
        return a_rel_mas * photocentre_factor(q, ms.flux_ratio(m1, q, "G"))

    q_lo = max(0.5, ms.mass.min() / m1)
    if a0(q_lo) < a0_target_mas:
        return None
    q = brentq(lambda x: a0(x) - a0_target_mas, q_lo, 1 - 1e-9, xtol=1e-12)
    f_g, f_rp = ms.flux_ratio(m1, q, "G"), ms.flux_ratio(m1, q, "RP")
    sig = rv_signatures(m1, q * m1, period_yr, inc)
    return dict(q=q, flux_ratio_g=f_g, flux_ratio_rp=f_rp,
                overluminosity_g=2.5 * np.log10(1 + f_g),
                blended_rv_amp_rp=sig["k_rel"] * photocentre_factor(q, f_rp), **sig)
