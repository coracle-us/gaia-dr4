"""Detection thresholds from the null calibration, and completeness from injection-recovery.

Inputs : results/calibration_null.csv (scripts/calibrate_null.py)
         results/injection_recovery.csv (scripts/injection_recovery.py), optional
Outputs: results/thresholds.csv          FAP per star -> threshold, per statistic
         results/completeness.csv        detection fraction in bins of (dchi2_exp, P, geometry)
         figures/completeness.png
Thresholds: empirical quantile where the sample allows it, else an exponential tail fitted
above the 95th percentile (P(>x) = 0.05 exp(-(x - x95) / lambda), lambda = mean excess).
Design FAP: 1e-4 per star -> ~0.08 expected false alarms over ~790 hosts (CLAUDE.md: FAP << 1/790).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
FAPS = (1e-2, 1e-3, 1e-4)
DESIGN_FAP = 1e-4
STATS = ("dlnl2_orbit", "dlnl2_acc7", "dlnl2_acc9")


def threshold(x: np.ndarray, fap: float) -> tuple[float, str]:
    x = np.sort(np.asarray(x))
    if fap * len(x) >= 10:
        return float(np.quantile(x, 1 - fap)), "empirical"
    x95 = np.quantile(x, 0.95)
    lam = np.mean(x[x > x95] - x95)
    return float(x95 + lam * np.log(0.05 / fap)), "exp-tail"


def null_thresholds(null: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for col in STATS:
        for fap in FAPS:
            thr, how = threshold(null[col], fap)
            rows.append(dict(statistic=col, fap_per_star=fap, threshold=thr, method=how,
                             n_null=len(null), max_null=null[col].max()))
    # stability across magnitude: design-FAP orbit threshold per G tercile
    for lo, hi in [(0, 11), (11, 13), (13, 25)]:
        sub = null[(null.g_mag >= lo) & (null.g_mag < hi)]
        if len(sub) > 200:
            thr, how = threshold(sub["dlnl2_orbit"], DESIGN_FAP)
            rows.append(dict(statistic=f"dlnl2_orbit[G {lo}-{hi}]", fap_per_star=DESIGN_FAP,
                             threshold=thr, method=how, n_null=len(sub),
                             max_null=sub["dlnl2_orbit"].max()))
    return pd.DataFrame(rows)


def completeness(inj: pd.DataFrame, thr_orbit: float, thr_acc: float) -> pd.DataFrame:
    inj = inj.copy()
    inj["det_orbit"] = inj["dlnl2_orbit"] > thr_orbit
    inj["det_acc"] = np.maximum(inj["dlnl2_acc7"], inj["dlnl2_acc9"]) > thr_acc
    inj["period_ok"] = np.abs(inj["period_fit"] / inj["period_true"] - 1) < 0.1
    d_bins = [3, 10, 30, 100, 300, 1000, 3000]
    p_bins = [0.2, 1.0, 2.0, 3.0, 4.0, 5.5, 8.0]
    inj["dbin"] = pd.cut(inj["dchi2_exp_actual"], d_bins)
    inj["pbin"] = pd.cut(inj["period_true"], p_bins)
    g = inj.groupby(["aligned", "pbin", "dbin"], observed=True)
    out = g.agg(n=("det_orbit", "size"), frac_orbit=("det_orbit", "mean"),
                frac_orbit_and_period=("period_ok", lambda s: np.mean(
                    s & inj.loc[s.index, "det_orbit"])),
                frac_acc=("det_acc", "mean"),
                frac_any=("det_orbit", lambda s: np.mean(s | inj.loc[s.index, "det_acc"]))
                ).reset_index()
    return out, inj


def plot(inj: pd.DataFrame, path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
    d_edges = np.logspace(np.log10(3), np.log10(3000), 13)
    d_mid = np.sqrt(d_edges[1:] * d_edges[:-1])
    for ax, aligned in zip(axes, (True, False)):
        sub = inj[inj.aligned == aligned]
        for (plo, phi), c in zip([(0.2, 2), (2, 4), (4, 5.5), (5.5, 8)],
                                 ["C0", "C1", "C2", "C3"]):
            s = sub[(sub.period_true >= plo) & (sub.period_true < phi)]
            idx = np.digitize(s.dchi2_exp_actual, d_edges) - 1
            frac = [s.det_orbit[idx == k].mean() if np.sum(idx == k) >= 5 else np.nan
                    for k in range(len(d_mid))]
            ax.semilogx(d_mid, frac, "o-", color=c, label=f"P {plo}-{phi} yr (orbit)")
            if plo >= 4:
                frac_a = [(s.det_orbit | s.det_acc)[idx == k].mean() if np.sum(idx == k) >= 5
                          else np.nan for k in range(len(d_mid))]
                ax.semilogx(d_mid, frac_a, "s--", color=c, alpha=0.6,
                            label=f"P {plo}-{phi} yr (orbit or accel.)")
        ax.axvline(50, color="0.6", lw=0.8, ls=":")
        ax.set_xlabel(r"expected full-orbit $\Delta\chi^2$  $N a_0^2(1+\cos^2 i)/4\sigma^2$")
        ax.set_title("aligned (i = 90 +- 3 deg)" if aligned else "isotropic inclination")
        ax.grid(alpha=0.3)
    axes[0].set_ylabel(f"detection fraction (FAP {DESIGN_FAP:g} per star)")
    axes[1].legend(fontsize=7, loc="lower right")
    fig.tight_layout()
    path.parent.mkdir(exist_ok=True)
    fig.savefig(path, dpi=130)


def main() -> None:
    null = pd.read_csv(RES / "calibration_null.csv")
    thr = null_thresholds(null)
    thr.to_csv(RES / "thresholds.csv", index=False, float_format="%.4g")
    print(thr.to_string(index=False))
    inj_path = RES / "injection_recovery.csv"
    if not inj_path.exists():
        return
    pick = lambda s: float(thr[(thr.statistic == s) & (thr.fap_per_star == DESIGN_FAP)].threshold.iloc[0])
    t_orb = pick("dlnl2_orbit")
    t_acc = max(pick("dlnl2_acc7"), pick("dlnl2_acc9"))   # conservative for "either"
    comp, inj = completeness(pd.read_csv(inj_path), t_orb, t_acc)
    comp.to_csv(RES / "completeness.csv", index=False, float_format="%.4g")
    with pd.option_context("display.width", 200, "display.max_rows", 200):
        print(comp.round(2).to_string(index=False))
    plot(inj, ROOT / "figures" / "completeness.png")


if __name__ == "__main__":
    main()
