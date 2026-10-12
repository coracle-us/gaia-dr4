"""Release-day vetting inputs: DPAC NSS solutions for our hosts, and Gaia epoch-RV constancy.

summarize_dpac: per source_id, the DPAC non-single-star solution types found in DR4 TAP tables
(nss_two_body_orbit, nss_acceleration_astro, ...) and the orbital period if any. Any match must be
compared with our result and cited (CLAUDE.md: cite, don't claim).
epoch_rv_vetting: Gaia epoch RVs (companions.gaia_rv) minus the transiting planet's signal
(literature K = pl_rvamp, P = pl_orbper, T_c = pl_tranmid; circular), tested for constancy;
gaia_rv_variable if p < 1e-3. Hosts without epoch RVs or literature K get NaN / no correction.
rv_known_period_test / known_period_rv_table: for named targets with periods known from photometry
(rv_known_period_d, rv_test_periods_d "P1;P2", reference epoch rv_t0_bjd;
e.g. TESS eclipses), the Gaia epoch-RV semi-amplitude at that period and its dchi2.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

from companions.epoch import DAYS_PER_YEAR, JD_REF_EPOCH
from companions.gaia_rv import rv_constancy, rvdata_from_gaia_epoch
from companions.orbit import solve_kepler

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


def rv_known_period_test(rv, period_days: float, t0_bjd: float = np.nan) -> dict:
    """Circular orbit at a period known from photometry (e.g. TESS eclipses), fitted to Gaia
    epoch RVs: v = c + a sin(2 pi t / P) + b cos(2 pi t / P), weighted least squares.
    Returns K = hypot(a, b) [m/s], its 1-sigma error (linearised; errors scaled by
    sqrt(chi2_red) when > 1), k_upper_ms = K + 3 sigma, and dchi2 against a constant (2 extra
    parameters). With t0_bjd (a primary eclipse or transit, companion in front): t0_offset_d,
    the time from T0 to the nearest fitted RV zero-crossing going from red- to blueshift, where
    the primary should be at conjunction (0 for a consistent eclipse; meaningful only when K is
    significant)."""
    p_yr = period_days / DAYS_PER_YEAR
    w = 1.0 / rv.err**2
    ph = 2 * np.pi * rv.t / p_yr
    design = np.column_stack([np.ones_like(rv.t), np.sin(ph), np.cos(ph)])
    a = design * np.sqrt(w)[:, None]
    coef, *_ = np.linalg.lstsq(a, rv.rv * np.sqrt(w), rcond=None)
    chi2_orbit = float(np.sum(w * (rv.rv - design @ coef) ** 2))
    mean = np.sum(w * rv.rv) / np.sum(w)
    chi2_const = float(np.sum(w * (rv.rv - mean) ** 2))
    cov = np.linalg.inv(a.T @ a) * max(1.0, chi2_orbit / max(len(rv.t) - 3, 1))
    k = float(np.hypot(coef[1], coef[2]))
    grad = np.array([0.0, coef[1], coef[2]]) / k if k > 0 else np.array([0.0, 1.0, 0.0])
    k_err = float(np.sqrt(grad @ cov @ grad))
    offset = np.nan
    if np.isfinite(t0_bjd):
        # v - c = K sin(theta + phi), theta = 2 pi t / P; decreasing zero-crossing at theta + phi = pi
        phi = np.arctan2(coef[2], coef[1])
        t_cross_d = (np.pi - phi) / (2 * np.pi) * period_days
        t0_d = t0_bjd - JD_REF_EPOCH                 # days from J2017.5, like rv.t (TDB vs TCB: s)
        offset = float((t_cross_d - t0_d + period_days / 2) % period_days - period_days / 2)
    return dict(n=len(rv.t), k_ms=k, k_err_ms=k_err, k_upper_ms=k + 3 * k_err,
                dchi2=chi2_const - chi2_orbit, t0_offset_d=offset)


def _true_anomaly(t_d, p_d, e, tp_d):
    big_e = solve_kepler(2 * np.pi * (t_d - tp_d) / p_d, e)
    return 2 * np.arctan2(np.sqrt(1 + e) * np.sin(big_e / 2), np.sqrt(1 - e) * np.cos(big_e / 2))


def rv_keplerian_fixed_period(rv, period_days: float, t0_bjd: float = np.nan, e_max: float = 0.9) -> dict:
    """Keplerian orbit at a fixed period fitted to Gaia epoch RVs:
    v = c + K [cos(nu + omega) + e cos(omega)], omega the primary's argument of periastron.
    For given (e, T_p) the model is linear in (c, K cos omega, -K sin omega); (e, T_p) come from a
    grid then least squares. Returns K [m/s], e, omega [deg], T_p [BJD], dchi2 against a
    constant, and with t0_bjd (primary eclipse / transit) t0_offset_d: the nearest primary
    conjunction (nu + omega = 90 deg, primary farthest) minus T0. For e ~ 0.5 the circular fit's
    phase is biased; use this one. Needs >= 6 epochs to have any freedom left."""
    t_d = rv.t * DAYS_PER_YEAR                                   # days from J2017.5
    w = 1.0 / rv.err**2
    sw = np.sqrt(w)

    def linear(e, tp):
        nu = _true_anomaly(t_d, period_days, e, tp)
        design = np.column_stack([np.ones_like(t_d), np.cos(nu), np.sin(nu)])
        coef, *_ = np.linalg.lstsq(design * sw[:, None], rv.rv * sw, rcond=None)
        return coef, (rv.rv - design @ coef) * sw

    best = (np.inf, 0.0, 0.0)
    for e in np.linspace(0, e_max, 19):
        for tp in np.linspace(0, period_days, 48, endpoint=False):
            chi2 = float(np.sum(linear(e, tp)[1] ** 2))
            if chi2 < best[0]:
                best = (chi2, e, tp)
    sol = least_squares(lambda x: linear(x[0], x[1])[1], [best[1], best[2]],
                        bounds=([0.0, best[2] - period_days], [e_max, best[2] + period_days]),
                        x_scale=[0.05, 0.02 * period_days])
    e, tp = float(sol.x[0]), float(sol.x[1])
    coef, resid = linear(e, tp)
    chi2 = float(np.sum(resid**2))
    mean = np.sum(w * rv.rv) / np.sum(w)
    chi2_const = float(np.sum(w * (rv.rv - mean) ** 2))
    k, omega = float(np.hypot(coef[1], coef[2])), float(np.arctan2(-coef[2], coef[1]))
    offset = np.nan
    if np.isfinite(t0_bjd):
        nu_c = np.pi / 2 - omega
        big_e = 2 * np.arctan(np.sqrt((1 - e) / (1 + e)) * np.tan(nu_c / 2))
        t_conj = tp + ((big_e - e * np.sin(big_e)) % (2 * np.pi)) / (2 * np.pi) * period_days
        offset = float((t_conj - (t0_bjd - JD_REF_EPOCH) + period_days / 2) % period_days - period_days / 2)
    return dict(n=len(rv.t), k_ms=k, e=e, omega_deg=float(np.rad2deg(omega) % 360),
                tp_bjd=tp + JD_REF_EPOCH, dchi2=chi2_const - chi2, t0_offset_d=offset)


def _periods(value) -> list[float]:
    if value is None or (isinstance(value, float) and not np.isfinite(value)):
        return []
    return [float(x) for x in str(value).split(";") if x.strip()]


def known_period_rv_table(tables: dict, hosts: pd.DataFrame) -> pd.DataFrame:
    """One row per (host, candidate period) for hosts with rv_test_periods_d ("P1;P2", days)
    and optional rv_t0_bjd: rv_known_period_test on the Gaia epoch RVs (n = 0 if none)."""
    rows = []
    for _, h in hosts.iterrows():
        sid = int(h["source_id"])
        rv = rvdata_from_gaia_epoch(tables[sid], sid) if sid in tables else None
        for p in _periods(h.get("rv_test_periods_d")):
            row = dict(source_id=sid, hostname=h.get("hostname"), period_d=p, n=0)
            if rv is not None and len(rv.t) >= 5:
                row.update(rv_known_period_test(rv, p, h.get("rv_t0_bjd", np.nan)))
            if rv is not None and len(rv.t) >= 6:
                kep = rv_keplerian_fixed_period(rv, p, h.get("rv_t0_bjd", np.nan))
                row.update({f"kep_{k}": v for k, v in kep.items() if k != "n"})
            rows.append(row)
    cols = ["source_id", "hostname", "period_d", "n", "k_ms", "k_err_ms", "k_upper_ms", "dchi2",
            "t0_offset_d", "kep_k_ms", "kep_e", "kep_omega_deg", "kep_tp_bjd", "kep_dchi2", "kep_t0_offset_d"]
    return pd.DataFrame(rows).reindex(columns=cols)
