"""Static HTML dashboard for a release-day run directory (scripts/release_day.py --out DIR).

    build_dashboard(DIR) -> DIR/dashboard/index.html (+ plots/*.png)

Reads whatever the run has written so far, so it works mid-run and after:
  status.json      stage (fetch / fit / vetting / done), counts, timings
  hosts.csv        every host in the run
  progress.jsonl   one cascade row per fitted host, appended as hosts finish
  results.csv      final table incl. DPAC / epoch-RV vetting (when stage = done)
  diagnostics/<source_id>.json   plot data for detections and named targets
                                 (companions.diagnostics); PNGs are cached in dashboard/plots/
                                 and embedded in the page, so index.html is self-contained
  injections.csv, population_summary.json   when present
Manual vetting decisions (independent re-derivation, impostor checks done by hand) come from a
CSV with columns hostname, check, status (pass / warn / fail / pending), note; they override the
automatic status of that check. Pages auto-refresh until the run is done.
"""
from __future__ import annotations

import base64
import html
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PREDICTIONS = ROOT / "results/dr4_predictions.csv"
PREDICTIONS_NOTE = ROOT / "notes/dr4_predictions.md"
REFRESH_S = 60
PLOTS = ("period_profile", "residuals_time", "residuals_phase", "sky_orbit")
CHECKS = ("significance", "stability", "impostor_twin", "stellar_companion", "gaia_epoch_rv",
          "dpac", "mass", "eccentricity", "rederivation")


# ---------------------------------------------------------------- vetting checklist

def _num(row, key):
    try:
        v = float(row.get(key, np.nan))
    except (TypeError, ValueError):
        return np.nan
    return v


def _flags(row) -> set:
    f = row.get("flags")
    return set(str(f).split(";")) - {""} if isinstance(f, str) else set()


def vetting_checklist(row: dict, thr: dict, manual: dict) -> list[dict]:
    """One entry per CLAUDE.md section-5 check: {check, status, detail}; status is pass / warn /
    fail / pending (needs work before a claim) / info / na. manual: {(hostname, check): (status,
    note)} overrides."""
    cls, flags = row.get("classification"), _flags(row)
    out = {}
    d_orb = _num(row, "dlnl2_orbit")
    if cls == "orbit":
        txt = f"dlnl2 = {d_orb:.1f} ({d_orb / thr['orbit']:.1f}x the FAP 1e-4 threshold {thr['orbit']:.1f})"
        out["significance"] = ("warn", txt + "; marginal (< 1.25x)") if "marginal" in flags else ("pass", txt)
    elif cls == "acceleration":
        a7, a9 = _num(row, "dlnl2_acc7"), _num(row, "dlnl2_acc9")
        out["significance"] = ("pass", f"dlnl2 acc7 = {a7:.1f} (T {thr['acc7']:.1f}), acc9 = {a9:.1f} "
                                       f"(T {thr['acc9']:.1f}), orbit = {d_orb:.1f}")
    else:
        out["significance"] = ("na", "not a detection")
    if "stability_confirmed" in flags:
        out["stability"] = ("pass", "orbit re-found in all five drop-one-fifth subsets")
    elif "stability_not_confirmed" in flags:
        out["stability"] = ("warn", "not re-found in every subset (weak evidence against: "
                                    "44-69% of genuine detections pass)")
    elif cls == "orbit":
        out["stability"] = ("pending", "stability check not run")
    else:
        out["stability"] = ("na", "orbit-only check")
    fdl = _num(row, "fraction_double_lined")
    single = "; single-lined in DR4 RVS" if fdl == 0 else ""
    if np.isfinite(fdl) and fdl > 0:
        out["impostor_twin"] = ("fail", f"double-lined in DR4 RVS (fraction_double_lined = {fdl:.2f}): "
                                        "a stellar binary")
    elif "twin_excluded_lit_rv" in flags:
        out["impostor_twin"] = ("pass", "literature RVs exclude a near-twin binary" + single)
    elif "impostor_not_excluded" in flags:
        out["impostor_twin"] = ("pending", "no literature RVs: needs DR4 SB2 / line-width check" + (
            single + " (does not exclude twins split by less than ~25 km/s)" if single else ""))
    else:
        out["impostor_twin"] = ("na", "")
    if "stellar_companion_check" in flags:
        out["stellar_companion"] = ("warn", "RUWE > 1.4 or prior HGCA / PMa / DR3 NSS evidence: "
                                            "check for a known stellar companion")
    else:
        out["stellar_companion"] = ("pass", "RUWE <= 1.4, no prior acceleration evidence")
    n_rv = _num(row, "gaia_rv_n")
    if "gaia_rv_n" not in row:
        out["gaia_epoch_rv"] = ("na", "epoch-RV vetting not run (simulated run, or stage not reached)")
    elif not np.isfinite(n_rv) or n_rv == 0:
        out["gaia_epoch_rv"] = ("pending", "no Gaia epoch RVs (yet)")
    elif str(row.get("gaia_rv_variable")).lower() in ("true", "1", "1.0"):
        out["gaia_epoch_rv"] = ("warn", f"RV-variable after removing the transiting planet "
                                        f"(p = {_num(row, 'gaia_rv_pvalue'):.2g}, N = {n_rv:.0f})")
    else:
        out["gaia_epoch_rv"] = ("pass", f"constant (p = {_num(row, 'gaia_rv_pvalue'):.2g}, N = {n_rv:.0f})")
    if "dpac_solutions" not in row:
        out["dpac"] = ("na", "not cross-matched (simulated run)")
    elif isinstance(row["dpac_solutions"], str) and row["dpac_solutions"]:
        out["dpac"] = ("info", f"DPAC NSS: {row['dpac_solutions']}; compare and cite")
    else:
        out["dpac"] = ("pass", "no DPAC NSS solution")
    if "stellar_mass" in flags:
        out["mass"] = ("fail", "dark-companion mass > 80 M_J: a star, not a planet")
    elif "brown_dwarf_mass" in flags:
        out["mass"] = ("info", "13-80 M_J (brown dwarf)")
    elif cls == "orbit":
        out["mass"] = ("pass", f"{_num(row, 'mass_mj'):.1f} M_J")
    else:
        out["mass"] = ("na", "")
    if "high_e_check" in flags:
        out["eccentricity"] = ("warn", "e > 0.8: noise solutions concentrate at e > 0.9")
    elif cls == "orbit":
        out["eccentricity"] = ("pass", f"e = {_num(row, 'e'):.2f}")
    else:
        out["eccentricity"] = ("na", "")
    out["rederivation"] = ("pending", "independent re-derivation (fresh script or separate agent) "
                                      "not yet recorded")
    host = row.get("hostname")
    checks = []
    for name in CHECKS:
        status, detail = out[name]
        if (host, name) in manual:
            status, note = manual[(host, name)]
            detail = f"manual: {note}"
        checks.append(dict(check=name, status=status, detail=detail))
    return checks


def overall_status(checks) -> str:
    statuses = {c["status"] for c in checks}
    for s in ("fail", "pending", "warn"):
        if s in statuses:
            return s
    return "pass"


def load_manual(path) -> dict:
    if path is None or not Path(path).exists():
        return {}
    df = pd.read_csv(path).fillna("")
    return {(r.hostname, r.check): (str(r.status).strip().lower(), str(r.note)) for r in df.itertuples()}


# ---------------------------------------------------------------- run state

def _read_json(path):
    return json.loads(Path(path).read_text()) if Path(path).exists() else None


def load_run(out: Path) -> dict:
    out = Path(out)
    status = _read_json(out / "status.json")
    summary = _read_json(out / "run_summary.json")
    if status is None:
        status = dict(stage="done" if summary else "not started")
        if summary:
            status.update(source=summary.get("source"), n_hosts=summary.get("n_hosts"),
                          git_commit=summary.get("git_commit"), thresholds=summary.get("thresholds"))
    hosts = pd.read_csv(out / "hosts.csv") if (out / "hosts.csv").exists() else None
    if status.get("stage") == "done" and (out / "results.csv").exists():
        res = pd.read_csv(out / "results.csv")
    else:
        rows = []
        if (out / "progress.jsonl").exists():
            for line in (out / "progress.jsonl").read_text().splitlines():
                if line.strip():
                    rows.append(json.loads(line))
        res = pd.DataFrame(rows)
        if len(res) and hosts is not None:
            extra = [c for c in hosts.columns if c not in res.columns or c == "source_id"]
            res = res.merge(hosts[extra], on="source_id", how="left")
    if "named_target" not in res.columns:
        res["named_target"] = False
    res["named_target"] = res["named_target"].fillna(False).astype(bool)
    for col in ("classification", "status", "flags"):
        if col not in res.columns:
            res[col] = np.nan
    thr = status.get("thresholds") or (summary or {}).get("thresholds")
    return dict(out=out, status=status, summary=summary, hosts=hosts, results=res, thresholds=thr,
                injections=pd.read_csv(out / "injections.csv") if (out / "injections.csv").exists() else None,
                population=_read_json(out / "population_summary.json"),
                named_rv=pd.read_csv(out / "named_rv.csv") if (out / "named_rv.csv").exists() else None)


def fetch_progress(status: dict, hosts) -> tuple[int, int] | None:
    """(n downloaded, n wanted) from the DataLink manifest while fetching."""
    path = status.get("fetch_manifest")
    if not path or hosts is None or not Path(path).exists():
        return None
    got = {int(s) for v in json.loads(Path(path).read_text()).values() for s in v}
    want = set(hosts["source_id"].astype("int64"))
    return len(got & want), len(want)


# ---------------------------------------------------------------- plots

def _stale(png: Path, src: Path) -> bool:
    return not png.exists() or png.stat().st_mtime < src.stat().st_mtime


def make_host_plots(diag_path: Path, plot_dir: Path, sid: int, label: str) -> list[str]:
    names = [f"{sid}_{p}.png" for p in PLOTS]
    if not any(_stale(plot_dir / n, diag_path) for n in names):
        return names
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    d = json.loads(diag_path.read_text())
    tr = {k: np.array(v, dtype=float) for k, v in d["transits"].items()}
    thr, orb, acc = d["thresholds"], d["orbit"], d["accel"]
    plot_dir.mkdir(parents=True, exist_ok=True)

    def save(fig, name):
        fig.tight_layout()
        fig.savefig(plot_dir / name, dpi=110)
        plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 3.2))
    pp = d.get("period_profile") or {}
    if pp:
        ax.plot(np.array(pp["period_yr"], float) * 365.25, np.array(pp["dchi2"], float), lw=0.8)
        # grid Delta-chi2 is at the single-star noise, not the final dlnl2: no threshold line
        ax.axvline(orb["period_yr"] * 365.25, color="C2", lw=0.8, label=f"best P = {orb['period_yr'] * 365.25:.1f} d")
        if pp["period_yr"][-1] / pp["period_yr"][0] > 3 or pp["period_yr"][0] / pp["period_yr"][-1] > 3:
            ax.set_xscale("log")
        ax.legend(fontsize=7)
    ax.set(xlabel="period [d]", ylabel="Δχ² (grid, single-star noise)",
           title=f"{label}: period search (final dlnl2 = {orb['dlnl2']:.1f}, threshold {thr['orbit']:.1f})")
    ax.title.set_fontsize(9)
    save(fig, names[0])

    fig, axes = plt.subplots(2, 1, figsize=(6, 4.6), sharex=True)
    for ax, key, txt in ((axes[0], "orbit", f"orbit fit, dlnl2 = {orb['dlnl2']:.1f}"),
                         (axes[1], "accel", f"{acc['order']}p acceleration fit, dlnl2 = {acc['dlnl2']:.1f}")):
        ax.errorbar(tr["t"] + 2017.5, tr[f"{key}_resid"], tr["err"], fmt=".", ms=3, lw=0.6, color="0.3")
        ax.plot(tr["t"] + 2017.5, tr[f"{key}_model"], "o", mfc="none", ms=4, color="C1", label="model")
        ax.axhline(0, color="0.6", lw=0.5)
        ax.set(ylabel="AL [mas]", title=txt)
        ax.title.set_fontsize(9)
    axes[0].legend(fontsize=7)
    axes[1].set_xlabel("epoch [yr]")
    save(fig, names[1])

    fig, ax = plt.subplots(figsize=(6, 3.2))
    for shift in (0, 1):
        ax.errorbar(tr["phase"] + shift, tr["orbit_resid"], tr["err"], fmt=".", ms=3, lw=0.6, color="0.3")
        ax.plot(tr["phase"] + shift, tr["orbit_model"], "o", mfc="none", ms=4, color="C1")
    ax.axhline(0, color="0.6", lw=0.5)
    ax.set(xlabel=f"orbital phase (P = {orb['period_yr'] * 365.25:.1f} d, e = {orb['e']:.2f})",
           ylabel="AL [mas]", title=f"{label}: per-transit AL signal (model = open circles)")
    ax.title.set_fontsize(9)
    save(fig, names[2])

    fig, ax = plt.subplots(figsize=(4, 4))
    sky = d["sky_orbit"]
    ax.plot(sky["ra_mas"], sky["dec_mas"], color="C0")
    ax.plot([0], [0], "k+")
    ax.set_aspect("equal")
    ax.invert_xaxis()
    ax.set(xlabel="Δα* [mas]", ylabel="Δδ [mas]", title="fitted photocentre orbit")
    ax.title.set_fontsize(9)
    save(fig, names[3])
    return names


def make_limits_plot(res: pd.DataFrame, plot_dir: Path) -> str | None:
    col = "m50_mj_P3"
    if col not in res or not np.isfinite(res[col]).any():
        return None
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    x = res.loc[res.classification == "none", col].dropna()
    x = x[x > 0]
    if not len(x):
        return None
    plot_dir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(6, 3.2))
    ax.hist(np.log10(x), bins=40)
    ax.set(xlabel="log10 50%-completeness mass at P = 3 yr [M_J]", ylabel="hosts",
           title=f"upper limits for {len(x)} non-detections")
    fig.tight_layout()
    fig.savefig(plot_dir / "limits_hist.png", dpi=110)
    plt.close(fig)
    return "limits_hist.png"


# ---------------------------------------------------------------- HTML

def esc(x) -> str:
    return html.escape("" if x is None or (isinstance(x, float) and np.isnan(x)) else str(x))


def fmt(v, nd=1) -> str:
    try:
        v = float(v)
    except (TypeError, ValueError):
        return esc(v)
    return "" if not np.isfinite(v) else f"{v:.{nd}f}"


def img(plot_dir: Path, name: str) -> str:
    """<img> with the PNG embedded, so the page is one self-contained file."""
    data = base64.b64encode((plot_dir / name).read_bytes()).decode("ascii")
    return f'<img alt="{esc(name)}" src="data:image/png;base64,{data}">'


def pill(status: str) -> str:
    return f'<span class="pill {esc(status)}">{esc(status)}</span>'


def table(headers, rows, table_id="") -> str:
    head = "".join(f"<th>{h}</th>" for h in headers)
    body = "".join("<tr>" + "".join(f"<td>{c}</td>" for c in r) + "</tr>" for r in rows)
    return f'<div class="scroll"><table{f' id="{table_id}"' if table_id else ''}><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></div>'


CSS = """
:root{--bg:#f7f7f5;--fg:#1d1d1f;--muted:#6b6b70;--card:#fff;--line:#e2e2e0;--accent:#2b59c3;
--pass:#1f7a3d;--warn:#a86400;--fail:#b3261e;--pending:#4a5a78;--info:#2b59c3;--na:#8a8a8f}
@media (prefers-color-scheme: dark){:root:not([data-theme="light"]){--bg:#161618;--fg:#ececee;
--muted:#9a9aa0;--card:#202023;--line:#333338;--accent:#7aa2ff;--pass:#5cc47e;--warn:#e0a040;
--fail:#ff7a70;--pending:#9fb0d0;--info:#7aa2ff;--na:#77777c}}
body{background:var(--bg);color:var(--fg);font:14px/1.45 system-ui,-apple-system,Segoe UI,sans-serif;
margin:0;padding:16px}
main{max-width:1200px;margin:0 auto}
h1{font-size:22px;margin:0 0 4px}h2{font-size:17px;margin:28px 0 8px}h3{font-size:15px;margin:0 0 6px}
.muted{color:var(--muted)}a{color:var(--accent)}
.card{background:var(--card);border:1px solid var(--line);border-radius:8px;padding:12px 14px;margin:10px 0}
.stats{display:flex;flex-wrap:wrap;gap:10px}.stat{background:var(--card);border:1px solid var(--line);
border-radius:8px;padding:8px 14px;min-width:110px}.stat b{display:block;font-size:20px}
.bar{height:10px;background:var(--line);border-radius:5px;overflow:hidden;margin:6px 0}
.bar>div{height:100%;background:var(--accent)}
.scroll{overflow-x:auto}table{border-collapse:collapse;width:100%;font-size:13px}
th,td{border-bottom:1px solid var(--line);padding:4px 8px;text-align:left;vertical-align:top;white-space:nowrap}
td.wrap{white-space:normal}th{font-weight:600;color:var(--muted)}
.pill{display:inline-block;padding:0 8px;border-radius:10px;font-size:12px;color:#fff;background:var(--na)}
.pill.pass{background:var(--pass)}.pill.warn{background:var(--warn)}.pill.fail{background:var(--fail)}
.pill.pending{background:var(--pending)}.pill.info{background:var(--info)}
.plots{display:grid;grid-template-columns:repeat(auto-fit,minmax(300px,1fr));gap:8px}
.plots img{width:100%;height:auto;background:#fff;border-radius:4px}
input[type=search]{padding:6px 8px;width:100%;max-width:320px;border:1px solid var(--line);
border-radius:6px;background:var(--card);color:var(--fg)}
"""

FILTER_JS = """
document.getElementById('hostfilter').addEventListener('input',function(e){
var q=e.target.value.toLowerCase();document.querySelectorAll('#allhosts tbody tr').forEach(function(r){
r.style.display=r.textContent.toLowerCase().indexOf(q)>=0?'':'none';});});
"""


def _progress_section(run) -> str:
    st, res = run["status"], run["results"]
    stage = st.get("stage", "?")
    n_jobs = st.get("n_jobs") or st.get("n_hosts") or 0
    n_done = int((res["status"].notna()).sum()) if len(res) else 0
    n_err = int(res["status"].astype(str).str.startswith("error").sum()) if len(res) else 0
    parts = [f'<div class="stats"><div class="stat"><span class="muted">stage</span><b>{esc(stage)}</b></div>',
             f'<div class="stat"><span class="muted">hosts fitted</span><b>{n_done} / {n_jobs}</b></div>',
             f'<div class="stat"><span class="muted">errors</span><b>{n_err}</b></div>']
    fp = fetch_progress(st, run["hosts"])
    if fp:
        parts.append(f'<div class="stat"><span class="muted">downloaded</span><b>{fp[0]} / {fp[1]}</b></div>')
    eta = ""
    try:
        t0, t1 = pd.Timestamp(st["fit_started"]), pd.Timestamp(st["updated"])
        if stage == "fit" and 0 < n_done < n_jobs:
            rate = n_done / max((t1 - t0).total_seconds(), 1)
            eta = f" · about {(n_jobs - n_done) / rate / 60:.0f} min to go"
    except (KeyError, ValueError, TypeError):
        pass
    parts.append("</div>")
    frac = n_done / n_jobs if n_jobs else 0
    parts.append(f'<div class="bar"><div style="width:{100 * frac:.1f}%"></div></div>'
                 f'<div class="muted">updated {esc(st.get("updated", ""))}{esc(eta)}</div>')
    errs = res[res["status"].astype(str).str.startswith("error")] if len(res) else res
    if len(errs):
        parts.append("<h3>Errors</h3>" + table(["host", "source_id", "error"],
                     [[esc(r.hostname), esc(r.source_id), f'<span class="wrap">{esc(r.status)}</span>']
                      for r in errs.itertuples()]))
    return "".join(parts)


def _post(row, key, nd=1):
    md = _num(row, f"post_{key}_p50")
    if not np.isfinite(md):
        return fmt(row.get(key), nd)
    return f"{md:.{nd}f} <span class=\"muted\">[{_num(row, f'post_{key}_p16'):.{nd}f}, {_num(row, f'post_{key}_p84'):.{nd}f}]</span>"


def build_dashboard(out, manual_path=None, predictions_path=PREDICTIONS) -> Path:
    run = load_run(Path(out))
    out = run["out"]
    dash = out / "dashboard"
    plot_dir = dash / "plots"
    dash.mkdir(parents=True, exist_ok=True)
    res, st = run["results"], run["status"]
    thr = run["thresholds"] or dict(orbit=np.nan, acc7=np.nan, acc9=np.nan)
    if manual_path is None:
        manual_path = out / "manual_vetting.csv"
    manual = load_manual(manual_path)
    done = st.get("stage") == "done"
    body = []

    body.append(f"<h1>Gaia DR4 companions: release-day run</h1><div class='muted'>{esc(out)} · source "
                f"{esc(st.get('source', '?'))} · commit {esc(str(st.get('git_commit', ''))[:10])} · "
                f"thresholds orbit {fmt(thr['orbit'])} / acc7 {fmt(thr['acc7'])} / acc9 {fmt(thr['acc9'])}"
                f" (FAP 1e-4 per star) · manual vetting: {esc(manual_path)}</div>")
    body.append("<h2>Progress</h2>" + _progress_section(run))

    counts = res["classification"].value_counts().to_dict() if len(res) else {}
    body.append("<div class='stats'>" + "".join(
        f"<div class='stat'><span class='muted'>{esc(k)}</span><b>{counts.get(k, 0)}</b></div>"
        for k in ("orbit", "acceleration", "none")) + "</div>")

    cand = res[res["classification"].isin(["orbit", "acceleration"])].copy() if len(res) else res
    show_plots = {}
    for r in pd.concat([cand, res[res["named_target"]]]).drop_duplicates("source_id").itertuples() if len(res) else []:
        diag = out / "diagnostics" / f"{int(r.source_id)}.json"
        if diag.exists():
            try:
                show_plots[int(r.source_id)] = make_host_plots(diag, plot_dir, int(r.source_id), str(r.hostname))
            except Exception as err:                  # a bad plot must not kill the page
                show_plots[int(r.source_id)] = f"plot failed: {err}"

    body.append("<h2>Candidates</h2>")
    if not len(cand):
        body.append("<div class='card muted'>No orbit or acceleration detections yet.</div>")
    else:
        cand["_score"] = np.fmax(cand["dlnl2_orbit"] / thr["orbit"],
                                 np.fmax(cand["dlnl2_acc7"] / thr["acc7"], cand["dlnl2_acc9"] / thr["acc9"]))
        cand = cand.sort_values("_score", ascending=False)
        rows, details = [], []
        for r in cand.to_dict("records"):
            checks = vetting_checklist(r, thr, manual)
            ov = overall_status(checks)
            sid = int(r["source_id"])
            rows.append([f'<a href="#host-{sid}">{esc(r["hostname"])}</a>', esc(r["classification"]), pill(ov),
                         fmt(r["dlnl2_orbit"]), f'{fmt(r["dlnl2_acc7"])} / {fmt(r["dlnl2_acc9"])}',
                         _post(r, "period_yr", 2), _post(r, "e", 2), _post(r, "mass_mj", 1),
                         _post(r, "min_mutual_inc_deg", 0), esc(r.get("flags"))])
            plots = show_plots.get(sid)
            if isinstance(plots, list):
                imgs = "".join(img(plot_dir, p) for p in plots)
            else:
                imgs = f"<div class='muted'>{esc(plots or 'no diagnostics file')}</div>"
            details.append(
                f'<div class="card" id="host-{sid}"><h3>{esc(r["hostname"])} {pill(ov)}</h3>'
                f'<div class="muted">source_id {sid} · G {fmt(r.get("phot_g_mean_mag"), 2)} · '
                f'{fmt(r.get("dist_pc"), 0)} pc · {esc(r["classification"])} · P {fmt(r.get("period_yr"), 3)} yr · '
                f'a0 {fmt(r.get("a0_mas"), 3)} mas · i {fmt(r.get("inc_deg"), 0)}° · '
                f'jitter {fmt(r.get("s_ccd"), 3)} / {fmt(r.get("s_tr"), 3)} mas</div>'
                + table(["check", "status", "detail"],
                        [[esc(c["check"]), pill(c["status"]), f'<span class="wrap">{esc(c["detail"])}</span>']
                         for c in checks])
                + f'<div class="plots">{imgs}</div></div>')
        body.append(table(["host", "class", "vetting", "dlnl2 orbit", "acc7 / acc9", "P [yr]", "e",
                           "mass [M_J]", "min ΔI [°]", "flags"], rows))
        body.append("<p class='muted'>Posterior medians with 16-84% ranges where sampled. A claim needs "
                    "every check at pass (or info), including the manual independent re-derivation.</p>")
        body.extend(details)

    named = res[res["named_target"]] if len(res) else res
    if len(named):
        body.append("<h2>Named targets vs pre-registered predictions</h2>")
        note = os.path.relpath(PREDICTIONS_NOTE, dash) if PREDICTIONS_NOTE.exists() else None
        if note:
            body.append(f'<p class="muted">Predictions written 2026-10-12 before release: '
                        f'<a href="{esc(Path(note).as_posix())}">notes/dr4_predictions.md</a>.</p>')
        preds = pd.read_csv(predictions_path) if predictions_path and Path(predictions_path).exists() else None
        for r in named.to_dict("records"):
            sid = int(r["source_id"])
            obs = [("classification", esc(r.get("classification"))),
                   ("status", esc(r.get("status"))),
                   ("dlnl2 orbit (period window)", fmt(r.get("dlnl2_orbit"))),
                   ("dlnl2 acc7 / acc9", f'{fmt(r.get("dlnl2_acc7"))} / {fmt(r.get("dlnl2_acc9"))}'),
                   ("Gaia epoch RVs: N, p(constant), rms [m/s]",
                    f'{fmt(r.get("gaia_rv_n"), 0)}, {fmt(r.get("gaia_rv_pvalue"), 3)}, {fmt(r.get("gaia_rv_rms_ms"), 0)}'),
                   ("RV at known period: K [m/s], Δχ²",
                    f'{fmt(r.get("gaia_rv_known_period_k_ms"), 0)}, {fmt(r.get("gaia_rv_known_period_dchi2"))}'),
                   ("DPAC NSS", esc(r.get("dpac_solutions"))),
                   ("DR4 RVS double-lined fraction", fmt(r.get("fraction_double_lined"), 2))]
            diag_path = out / "diagnostics" / f"{sid}.json"
            if diag_path.exists():
                o = json.loads(diag_path.read_text())["orbit"]
                sig = "" if _num(r, "dlnl2_orbit") > thr["orbit"] else " (not significant)"
                obs.insert(3, (f"best orbit in the search window{sig}: P [d], a0 [mas], e, i [°]",
                               f'{fmt(o["period_yr"] * 365.25, 2)}, {fmt(o.get("a0_mas"), 3)}, '
                               f'{fmt(o["e"], 2)}, {fmt(o.get("inc_deg"), 0)}'))
            block = [f'<div class="card"><h3>{esc(r["hostname"])}</h3>',
                     table(["observed in DR4", "value"], [[k, v] for k, v in obs])]
            nrv = run["named_rv"]
            if nrv is not None and len(nrv[nrv["source_id"] == sid]):
                block.append("<div class='muted' style='margin-top:8px'>Gaia epoch RVs at the photometric "
                             "period(s): circular fit; T0 offset = fitted conjunction minus the eclipse/transit "
                             "epoch (0 if consistent)</div>" + table(
                    ["P [d]", "N", "K [m/s]", "K 3σ upper [m/s]", "Δχ² vs constant", "T0 offset [d]"],
                    [[fmt(x["period_d"], 3), fmt(x["n"], 0),
                      f'{fmt(x.get("k_ms"), 0)} ± {fmt(x.get("k_err_ms"), 0)}', fmt(x.get("k_upper_ms"), 0),
                      fmt(x.get("dchi2"), 1), fmt(x.get("t0_offset_d"), 2)]
                     for x in nrv[nrv["source_id"] == sid].to_dict("records")]))
            if preds is not None:
                p = preds[preds["target"] == r["hostname"]]
                if len(p):
                    cols = [c for c in ("hypothesis", "period_d", "m2_mjup", "a0_mas", "dchi2_exp", "k1_kms",
                                        "line_split_kms") if c in p.columns]
                    block.append("<div class='muted' style='margin-top:8px'>predicted</div>" + table(
                        cols, [[esc(x["hypothesis"])] + [fmt(x[c], 3) for c in cols[1:]]
                               for x in p.to_dict("records")]))
            plots = show_plots.get(sid)
            if isinstance(plots, list):
                block.append('<div class="plots">' + "".join(img(plot_dir, p) for p in plots) + "</div>")
            block.append("</div>")
            body.append("".join(block))

    if run["injections"] is not None and len(res):
        inj = run["injections"].merge(res, on="source_id", how="left")
        ok = (inj["classification"] == "orbit") & ((inj["period_yr"] / inj["period_true_yr"] - 1).abs() < 0.1)
        body.append(f"<h2>Injected planets (simulated run)</h2><p>{int(ok.sum())} of {len(inj)} recovered "
                    f"as orbits with the period within 10%.</p>")
        body.append(table(["host", "P true [yr]", "P fit [yr]", "mass true [M_J]", "mass fit [M_J]", "class", "recovered"],
                          [[esc(x["hostname"]), fmt(x["period_true_yr"], 2), fmt(x.get("period_yr"), 2),
                            fmt(x["mass_true_mj"], 1), fmt(x.get("mass_mj"), 1), esc(x.get("classification")),
                            pill("pass" if k else "fail")] for x, k in zip(inj.to_dict("records"), ok)]))

    if done:
        body.append("<h2>Population</h2>")
        lim = make_limits_plot(res, plot_dir)
        if run["population"]:
            occ = run["population"].get("occurrence", {})
            pop = run["population"].get("population", {})
            body.append(f"<p class='muted'>Companions {pop.get('mass_mj')} M_J, {pop.get('period_yr')} yr.</p>")
            body.append(table(["sample", "hosts", "detected", "effective N", "occurrence (median)", "16-84%", "95% upper"],
                              [[esc(k), v["n_hosts"], v["n_detected"], fmt(v["effective_n"]),
                                fmt(v["median"], 3), f'{fmt(v["p16"], 3)}-{fmt(v["p84"], 3)}', fmt(v["upper95"], 3)]
                               for k, v in occ.items()]))
        else:
            body.append("<p class='muted'>Run scripts/population_analysis.py on results.csv to add occurrence rates.</p>")
        if lim:
            body.append(f'<div class="plots">{img(plot_dir, lim)}</div>')

    if len(res):
        body.append("<h2>All hosts</h2><input type='search' id='hostfilter' placeholder='filter hosts, flags, status'>")
        allh = res.copy()
        allh["_score"] = allh["dlnl2_orbit"].astype(float) / thr["orbit"]
        allh = allh.sort_values("_score", ascending=False)
        cols = ["phot_g_mean_mag", "dist_pc", "m50_mj_P3"]
        body.append(table(["host", "G", "d [pc]", "status", "class", "dlnl2 orbit", "acc7", "acc9",
                           "m50 at 3 yr [M_J]", "flags"],
                          [[esc(x["hostname"]), fmt(x.get(cols[0]), 2), fmt(x.get(cols[1]), 0),
                            f'<span class="wrap">{esc(x.get("status"))}</span>', esc(x.get("classification")),
                            fmt(x.get("dlnl2_orbit")), fmt(x.get("dlnl2_acc7")), fmt(x.get("dlnl2_acc9")),
                            fmt(x.get(cols[2])), esc(x.get("flags"))] for x in allh.to_dict("records")],
                          table_id="allhosts"))

    refresh = "" if done else f'<meta http-equiv="refresh" content="{REFRESH_S}">'
    page = (f'<!doctype html><html lang="en"><head><meta charset="utf-8">'
            f'<meta name="viewport" content="width=device-width,initial-scale=1">{refresh}'
            f"<title>DR4 run status</title><style>{CSS}</style></head><body><main>"
            + "".join(body) + f"</main><script>{FILTER_JS if len(res) else ''}</script></body></html>")
    path = dash / "index.html"
    tmp = dash / "index.html.tmp"
    tmp.write_text(page, encoding="utf-8")
    os.replace(tmp, path)                             # never serve a half-written page
    return path
