"""Load ESA's DR4 pre-release epoch astrometry (12 sources, June 2026) and plot one star.

Input : data/raw/dr4_prerelease/epoch_astrometry/GAIA_DR4_PRERELEASE_EPOCH_ASTROMETRY_RAW.xml
        (from gaia-dr4-prerelease-epoch-astrometry_2026-06-26.zip, see SHA256SUMS.txt)
Output: figures/prerelease_<source_id>.png, and a per-source summary printed to stdout.

Top panel: raw along-scan centroid positions (local-plane w, relative to ra0/dec0) per CCD.
Bottom panel: residuals from ESA gaiasupdate's DR4-like 6-parameter single-star fit,
averaged per FoV transit. An unmodelled orbit shows up here.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io.votable import parse_single_table
from gaiasupdate.epoch_astrometry import GaiaEpochAstrometryArchive

ROOT = Path(__file__).resolve().parents[1]
XML = ROOT / "data/raw/dr4_prerelease/epoch_astrometry/GAIA_DR4_PRERELEASE_EPOCH_ASTROMETRY_RAW.xml"
FIG = ROOT / "figures"
NS_PER_YR = 365.25 * 86400e9
NAMES = {3937211745905473024: "HD 114762", 1457486023639239296: "Gaia-4",
         4318465066420528000: "Gaia BH3", 4181040337841125632: "HD 183633",
         1663617687609809280: "LSPM J1324+6124"}


def tcb_ns_to_year(t_ns):
    """obs_time_tcb is ns since J2010.0 (TCB); Julian-year approximation is fine for plots."""
    return 2010.0 + np.asarray(t_ns, dtype=float) / NS_PER_YR


def flat(col, fill=np.nan, dtype=float):
    """Concatenate a per-CCD array column (one masked array per FoV transit row)."""
    return np.concatenate([np.ma.filled(np.ma.asarray(x).astype(dtype), fill) for x in col])


def summarize(t) -> pd.DataFrame:
    rows = []
    for sid in np.unique(t["source_id"]):
        s = t[t["source_id"] == sid]
        used = flat(s["used_by_agis_al"], False, bool)
        err = flat(s["centroid_pos_error_al"])
        tt = flat(s["obs_time_tcb"])
        rows.append(dict(
            source_id=int(sid), name=NAMES.get(int(sid), ""), g_mag=float(np.median(s["g_mag"])),
            n_fov_transits=len(s), n_ccd_used_al=int(used.sum()),
            med_ccd_err_used_mas=float(np.nanmedian(err[used])),
            agis_source_excess_noise_mas=float(s["agis_source_excess_noise"][0]),
            first=float(tcb_ns_to_year(np.nanmin(tt))), last=float(tcb_ns_to_year(np.nanmax(tt)))))
    return pd.DataFrame(rows)


def plot_source(t, sid: int) -> Path:
    df = t[t["source_id"] == sid].to_pandas()
    fit = GaiaEpochAstrometryArchive.supdate(df, sid)   # DR4-like 6p_constrained_colour
    ts = fit["timestamps_of_used_observations"]
    res = pd.DataFrame({"t_ns": ts["obsTimeTcb"].to_numpy(dtype=float),
                        "res": fit["residuals"], "var": fit["total_variance"]})
    # group CCD measurements into FoV transits (CCDs of one transit lie within ~1 min)
    res = res.sort_values("t_ns")
    res["transit"] = (res["t_ns"].diff() > 120e9).cumsum()
    res["w"] = 1.0 / res["var"]
    per = res.groupby("transit").apply(
        lambda g: pd.Series({"t": tcb_ns_to_year(g["t_ns"].mean()),
                             "res": np.average(g["res"], weights=g["w"]),
                             "err": 1.0 / np.sqrt(g["w"].sum())}), include_groups=False)

    # raw per-CCD AL positions (all CCDs, used and unused)
    sel = t[t["source_id"] == sid]
    tt, w = flat(sel["obs_time_tcb"]), flat(sel["centroid_pos_al"])
    used = flat(sel["used_by_agis_al"], False, bool)

    plx, pmra, pmdec = fit["parameters"][2:5]
    name = NAMES.get(sid, str(sid))
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(8, 6.5), sharex=True,
                                 gridspec_kw={"height_ratios": [1, 1]})
    a1.plot(tcb_ns_to_year(tt[~used]), w[~used], "x", color="0.6", ms=4, label="not used by AGIS")
    a1.plot(tcb_ns_to_year(tt[used]), w[used], ".", color="C0", ms=3, label="used by AGIS (AL)")
    a1.set_ylabel("centroid_pos_al  w  [mas]")
    a1.set_title(f"{name}  (Gaia DR4 pre-release, source_id {sid})\n"
                 f"{len(df)} FoV transits, {used.sum()} CCD AL measurements used")
    a1.legend(loc="best", fontsize=8)
    a2.axhline(0, color="0.5", lw=0.8)
    a2.errorbar(per["t"], per["res"], yerr=per["err"], fmt="o", ms=3, color="C3", lw=0.8)
    a2.set_ylabel("AL residual per FoV transit [mas]")
    a2.set_xlabel("year (TCB)")
    a2.set_title(f"residuals from DR4-like 6-parameter fit (gaiasupdate): "
                 f"plx={plx:.3f} mas, pm=({pmra:.2f}, {pmdec:.2f}) mas/yr", fontsize=9)
    fig.tight_layout()
    FIG.mkdir(exist_ok=True)
    out = FIG / f"prerelease_{sid}.png"
    fig.savefig(out, dpi=130)
    rms = np.sqrt(np.average(per["res"] ** 2))
    print(f"{name}: plx={plx:.4f} mas pmra={pmra:.3f} pmdec={pmdec:.3f}; "
          f"per-transit residual RMS={rms:.3f} mas vs median per-transit error "
          f"{per['err'].median():.3f} mas -> {out.relative_to(ROOT)}")
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--source-id", type=int, default=1457486023639239296)  # Gaia-4
    args = ap.parse_args()
    t = parse_single_table(XML).to_table()
    with pd.option_context("display.width", 200, "display.max_columns", 20):
        print(summarize(t).round(3).to_string(index=False))
    plot_source(t, args.source_id)


if __name__ == "__main__":
    main()
