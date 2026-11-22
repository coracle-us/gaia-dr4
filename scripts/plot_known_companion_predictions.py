"""Figure for the pre-DR4 prediction note.

Detection probability in simulated DR4 epoch astrometry (results/known_companion_predictions.csv,
scripts/known_companion_predictions.py) against the naive expected full-orbit Delta chi^2,
with the injection-recovery completeness for P < 5.5 yr (results/completeness.csv) for reference.
Writes figures/known_companion_predictions.png and .pdf.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

ROOT = Path(__file__).resolve().parents[1]
RES, FIG = ROOT / "results", ROOT / "figures"
# label position in data coordinates (x, y, horizontal alignment); a thin line joins label and point
LABEL = {"TOI-2537": (70, 1.03, "right"), "WASP-132": (45, 0.95, "right"), "WASP-34": (35, 0.86, "right"),
         "K2-232": (32, 0.62, "right"), "KOI-12": (330, 0.55, "left"), "TOI-2295": (80, 0.38, "left"),
         "TOI-4562": (70, 0.06, "left"), "HAT-P-17": (260, 0.30, "right"), "WASP-68": (1150, 0.19, "left"),
         "HAT-P-2": (4200, 0.30, "left"), "HD 118203": (2300, 0.66, "left"), "WASP-53": (2300, 0.80, "left"),
         "WASP-99": (7500, 0.84, "left"), "WASP-8": (11000, 0.93, "left"), "TOI-201": (2300, 1.04, "left"),
         "WASP-81": (700, 1.04, "right")}


def main():
    s = pd.read_csv(RES / "known_companion_predictions.csv")
    c = pd.read_csv(RES / "completeness.csv")
    c = c[~c.aligned].copy()
    c["plo"] = c.pbin.str.extract(r"\(([\d.]+),").astype(float)
    c["dlo"] = c.dbin.str.extract(r"\(([\d.]+),").astype(float)
    c["dhi"] = c.dbin.str.extract(r", ([\d.]+)\]").astype(float)
    c["dmid"] = np.sqrt(c.dlo * c.dhi)
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    for sel, lab, ls in [(c.plo < 5.4, r"injections, $P<5.5$ yr", "-"),
                         (c.plo >= 5.4, r"injections, $5.5<P<8$ yr", "--")]:
        g = c[sel].groupby("dmid")[["frac_any", "n"]].apply(
            lambda d: np.average(d.frac_any, weights=d.n))
        ax.plot(g.index, g.values, ls, color="0.55", lw=1.2, label=lab, zorder=1)
    x = s.dchi2_exp_med.clip(lower=0.02)
    err = np.sqrt(s.p_any * (1 - s.p_any) / s.n_draws)
    ax.errorbar(x, s.p_any, yerr=err, fmt="none", ecolor="0.7", lw=0.8, zorder=2)
    sc = ax.scatter(x, s.p_any, c=s.period_yr, cmap="viridis", norm=LogNorm(0.3, 25), s=34,
                    edgecolor="k", lw=0.4, zorder=3)
    ph = s[np.isfinite(s.p_any_phase) & ((s.p_any_phase - s.p_any).abs() > 0.1)]
    for _, r in ph.iterrows():
        ax.annotate("", xy=(r.dchi2_exp_med, r.p_any_phase), xytext=(r.dchi2_exp_med, r.p_any),
                    arrowprops=dict(arrowstyle="->", color="C3", lw=1.1), zorder=4)
    for _, r in s.iterrows():
        if r.hostname in LABEL:
            lx, ly, ha = LABEL[r.hostname]
            ax.annotate(r.hostname, (r.dchi2_exp_med, r.p_any), xytext=(lx, ly), ha=ha, va="center",
                        fontsize=6.5, arrowprops=dict(arrowstyle="-", color="0.6", lw=0.5,
                                                      shrinkA=1, shrinkB=3), zorder=5)
    ax.plot([], [], color="C3", lw=1.1, label="with published periastron time")
    ax.set_xscale("log")
    ax.set_xlim(0.02, 3e4)
    ax.set_ylim(-0.04, 1.06)
    ax.set_xlabel(r"naive expected full-orbit $\Delta\chi^2$")
    ax.set_ylabel("detection probability in simulated DR4")
    cb = fig.colorbar(sc, ax=ax, pad=0.01)
    cb.set_label("companion period (yr)")
    ax.legend(fontsize=7, loc="upper left", frameon=False)
    fig.tight_layout()
    FIG.mkdir(exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"known_companion_predictions.{ext}", dpi=200)


if __name__ == "__main__":
    main()
