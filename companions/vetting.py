"""Release-day vetting inputs: DPAC NSS solutions for our hosts, and Gaia epoch-RV constancy.

summarize_dpac: per source_id, the DPAC non-single-star solution types found in DR4 TAP tables
(nss_two_body_orbit, nss_acceleration_astro, ...) and the orbital period if any. Any match must be
compared with our result and cited (CLAUDE.md: cite, don't claim).
epoch_rv_vetting: Gaia epoch RVs (companions.gaia_rv) minus the transiting planet's signal
(literature K = pl_rvamp, P = pl_orbper, T_c = pl_tranmid; circular), tested for constancy;
gaia_rv_variable if p < 1e-3. Hosts without epoch RVs or literature K get NaN / no correction.
rv_known_period_test: for named targets with a period known from photometry (rv_known_period_d,
e.g. TESS eclipses), the Gaia epoch-RV semi-amplitude at that period and its dchi2.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from companions.epoch import DAYS_PER_YEAR, JD_REF_EPOCH
from companions.gaia_rv import rv_constancy, rvdata_from_gaia_epoch

RV_VARIABLE_PVALUE = 1e-3


def summarize_dpac(tables: dict, source_ids) -> pd.DataFrame:
    rows = []
    for sid in source_ids:
        types, period = [], np.nan
        for df in tables.values():
            m = df[df["source_id"] == sid]
            if len(m):
                types += [str(x) for x in m["nss_solution_type"]]
                if "period" in m and np.isfinite(m["period"]).any():
                    period = float(m["period"].dropna().iloc[0])
        rows.append(dict(source_id=int(sid), dpac_solutions=";".join(dict.fromkeys(types)),
                         dpac_period_d=period))
    return pd.DataFrame(rows)


def epoch_rv_vetting(tables: dict, hosts: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, h in hosts.iterrows():
        sid = int(h["source_id"])
        row = dict(source_id=sid, gaia_rv_n=0, gaia_rv_pvalue=np.nan, gaia_rv_rms_ms=np.nan,
                   gaia_rv_variable=False)
        if sid in tables:
            rv = rvdata_from_gaia_epoch(tables[sid], sid)
            inner = []
            if all(np.isfinite(h.get(c, np.nan)) for c in ("pl_orbper", "pl_tranmid", "pl_rvamp")):
                inner = [(h["pl_orbper"] / DAYS_PER_YEAR,
                          (h["pl_tranmid"] - JD_REF_EPOCH) / DAYS_PER_YEAR, h["pl_rvamp"])]
            if len(rv.t) >= 3:
                c = rv_constancy(rv, inner=inner)
                row.update(gaia_rv_n=len(rv.t), gaia_rv_pvalue=c.pvalue, gaia_rv_rms_ms=c.rms_ms,
                           gaia_rv_variable=bool(c.pvalue < RV_VARIABLE_PVALUE))
                p_known = h.get("rv_known_period_d", np.nan)
                if np.isfinite(p_known) and len(rv.t) >= 5:
                    kp = rv_known_period_test(rv, p_known)
                    row.update(gaia_rv_known_period_k_ms=kp["k_ms"], gaia_rv_known_period_dchi2=kp["dchi2"])
        rows.append(row)
    return pd.DataFrame(rows)


def rv_known_period_test(rv, period_days: float) -> dict:
    """Circular orbit at a period known from photometry (e.g. TESS eclipses), fitted to Gaia
    epoch RVs: v = c + a sin(2 pi t / P) + b cos(2 pi t / P), weighted least squares.
    Returns K = hypot(a, b) [m/s] and dchi2 against a constant (2 extra parameters)."""
    p_yr = period_days / DAYS_PER_YEAR
    w = 1.0 / rv.err**2
    ph = 2 * np.pi * rv.t / p_yr
    design = np.column_stack([np.ones_like(rv.t), np.sin(ph), np.cos(ph)])
    a = design * np.sqrt(w)[:, None]
    coef, *_ = np.linalg.lstsq(a, rv.rv * np.sqrt(w), rcond=None)
    chi2_orbit = float(np.sum(w * (rv.rv - design @ coef) ** 2))
    mean = np.sum(w * rv.rv) / np.sum(w)
    chi2_const = float(np.sum(w * (rv.rv - mean) ** 2))
    return dict(n=len(rv.t), k_ms=float(np.hypot(coef[1], coef[2])), dchi2=chi2_const - chi2_orbit)
