"""Unattended Gaia DR4 release-day driver. Each call is one 'tick'; launchd runs it every 10 min
(config/launchd/us.coracle.gaia-dr4-watch.plist) and the Claude check-in tasks run it too.

    python scripts/dr4_watch.py                 # one tick (does nothing before the window opens)
    python scripts/dr4_watch.py --status        # print the state, change nothing
    python scripts/dr4_watch.py --clear         # release a 'blocked' stage after a fix (retry probe)
    python scripts/dr4_watch.py --post          # (internal) the post-run chain, launched detached

Stages, kept in results/dr4_watch/state.json:
  polling      TAP has no `gaiadr4` schema yet (archive offline or DR4 not out).
  probing      DR4 tables visible; fetch Gaia-4 by DataLink, check columns, run the cascade on it
               (expect an orbit with P in probe_period_d). DataLink errors are retried every tick.
  blocked      The probe fetched data but the columns or the Gaia-4 sanity check failed, or tracked
               files are modified. Nothing launches until someone fixes it and runs --clear.
  main_running scripts/release_day.py --source datalink (detached; resumes if relaunched).
               Relaunched if it dies or leaves > max_missing_frac hosts without data.
  post_running calibrate_null_real, population_analysis, make_table (into results/dr4/), dashboard.
  done         Nothing more to do automatically; vetting (docs/RELEASE_DAY.md) is by hand.

Automatic: polling, the probe, downloads, fits, the post-run chain. Never automatic: editing
config or code, choosing thresholds_real, git commits or pushes, anything in paper/.
"""
import argparse
import datetime as dt
import json
import os
import shlex
import subprocess
import sys
import urllib.parse
import urllib.request
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
CONFIG = ROOT / "config/release_day.json"
PY = sys.executable


def utcnow() -> dt.datetime:
    return dt.datetime.now(dt.timezone.utc)


def stamp() -> str:
    return utcnow().strftime("%Y-%m-%dT%H:%M:%SZ")


def alive(pid) -> bool:
    if not pid:
        return False
    try:
        os.kill(int(pid), 0)
    except (ProcessLookupError, PermissionError, ValueError):
        return False
    try:                       # a zombie child of this process is not alive
        return os.waitpid(int(pid), os.WNOHANG) == (0, 0)
    except ChildProcessError:
        return True


class Watch:
    def __init__(self, args):
        self.args = args
        self.cfg = json.loads(args.config.read_text())
        self.ac = self.cfg["autorun"]
        self.dir = ROOT / self.ac["state_dir"]
        self.dir.mkdir(parents=True, exist_ok=True)
        self.path = self.dir / "state.json"
        self.s = json.loads(self.path.read_text()) if self.path.exists() else dict(stage="polling", history=[])
        self.out = ROOT / self.ac["out"]
        self.procs = args.procs or max(1, (os.cpu_count() or 8) - 2)

    # ---- bookkeeping ------------------------------------------------------------------------
    def log(self, msg, stage=None, notify=False):
        if stage:
            self.s["stage"] = stage
        self.s["message"] = msg
        self.s["history"].append(dict(t=stamp(), stage=self.s["stage"], msg=msg))
        self.s["history"] = self.s["history"][-200:]
        with open(self.dir / "watch.log", "a", encoding="utf-8") as f:
            f.write(f"{stamp()} [{self.s['stage']}] {msg}\n")
        if notify:
            self.notify(msg)

    def notify(self, msg):
        if self.args.quiet:
            return
        script = f'display notification {json.dumps(msg[:200])} with title "Gaia DR4 watch"'
        subprocess.run(["osascript", "-e", script], capture_output=True)

    def save(self):
        self.s["last_tick"] = stamp()
        self.s["procs"] = self.procs
        tmp = self.path.with_suffix(".tmp")
        tmp.write_text(json.dumps(self.s, indent=1, default=str))
        tmp.replace(self.path)
        lines = [f"# Gaia DR4 watch: {self.s['stage']}", "", f"Last tick: {self.s['last_tick']}",
                 f"Message: {self.s.get('message', '')}", ""]
        for k in ("attention", "probe", "main_launches", "main_pid", "post"):
            if self.s.get(k):
                lines.append(f"- {k}: {json.dumps(self.s[k], default=str)}")
        lines += ["", "Recent history:"] + [f"- {h['t']} [{h['stage']}] {h['msg']}"
                                            for h in self.s["history"][-15:]]
        (self.dir / "STATUS.md").write_text("\n".join(lines) + "\n")

    def attention(self, msg):
        att = self.s.setdefault("attention", [])
        if msg not in att:
            att.append(msg)
            self.log("ATTENTION: " + msg, notify=True)

    def spawn(self, cmd, logname):
        logf = open(self.dir / logname, "a", encoding="utf-8")
        logf.write(f"\n===== {stamp()} {shlex.join(cmd)}\n")
        logf.flush()
        env = dict(os.environ, PYTHONPATH=os.pathsep.join(
            [str(ROOT)] + [x for x in os.environ.get("PYTHONPATH", "").split(os.pathsep) if x]))
        p = subprocess.Popen(cmd, cwd=ROOT, stdout=logf, stderr=subprocess.STDOUT, env=env,
                             stdin=subprocess.DEVNULL, start_new_session=True)
        return p.pid

    def keep_awake(self):
        """One caffeinate for the release window, so the Mac does not sleep mid-download."""
        if sys.platform != "darwin" or alive(self.s.get("caffeinate_pid")) or self.args.no_caffeinate:
            return
        hours = float(self.ac.get("keep_awake_hours", 72))
        self.s["caffeinate_pid"] = self.spawn(["caffeinate", "-i", "-s", "-t", str(int(hours * 3600))],
                                              "caffeinate.log")
        self.log(f"caffeinate started for {hours:.0f} h (pid {self.s['caffeinate_pid']})")

    # ---- checks -----------------------------------------------------------------------------
    def tap_tables(self, schema):
        """Table names in `schema` from TAP_SCHEMA, or None if the archive does not answer."""
        q = f"select table_name from tap_schema.tables where schema_name = '{schema}'"
        url = self.ac["tap_url"] + "?" + urllib.parse.urlencode(
            dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", QUERY=q))
        try:
            with urllib.request.urlopen(url, timeout=90) as r:
                rows = r.read().decode().splitlines()[1:]
        except Exception as err:
            return None, f"{type(err).__name__}: {err}"
        return sorted(x.strip().strip('"') for x in rows if x.strip()), None

    def check_config_tables(self, tables):
        want = set(self.cfg.get("dpac_tables", [])) | set(self.cfg["control"]["nss_tables"])
        want |= {self.cfg["control"]["gaia_source_table"], self.cfg["rv_summary"]["table"]}
        have = {t.lower() for t in tables}
        missing = sorted(t for t in want if t.lower() not in have)
        if missing:
            self.attention(f"config tables not in TAP (vetting/control steps will fail until fixed "
                           f"in config/release_day.json): {missing}")
        return missing

    def fetch_probe(self):
        """Gaia-4 epoch astrometry as an astropy Table, trying each release string. Returns
        (table, release, errors); table is None if every attempt failed."""
        sid = int(self.ac["probe_source_id"])
        if self.args.probe_file:                     # offline test of the sanity path
            from astropy.io.votable import parse_single_table
            return parse_single_table(self.args.probe_file).to_table(), "file", []
        from companions.archive import DataLinkSource
        errors = []
        for rel in [self.cfg["data_release"]] + [r for r in self.ac.get("release_candidates", [])
                                                 if r != self.cfg["data_release"]]:
            raw = ROOT / self.ac["raw_dir"] / "probe" / rel.replace(" ", "_")
            src = DataLinkSource(rel, raw, retrieval_type=self.cfg["retrieval_type"],
                                 data_structure=self.cfg["data_structure"], batch_size=1,
                                 retries=2, backoff_s=20)
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                tables = src.fetch([sid])
            if sid in tables:
                return tables[sid], rel, errors
            errors.append(f"{rel!r}: " + "; ".join(str(x.message)[:300] for x in w) or "no rows")
        return None, None, errors

    def sanity(self, table):
        from companions.cascade import HostInfo, Thresholds, run_cascade
        from companions.epoch import EpochData
        sid = int(self.ac["probe_source_id"])
        ep = EpochData.from_table(table, sid)
        thr = Thresholds.from_csv(ROOT / self.cfg["thresholds_file"], fap=self.cfg["fap_per_star"])
        model = self.args.probe_model or self.cfg["model"]
        row = run_cascade(ep, HostInfo(hostname="Gaia-4", mstar=float(self.ac["probe_mstar"])), thr,
                          model=model, do_posterior=False).row()
        p_d = float(row.get("period_yr", float("nan"))) * 365.25
        lo, hi = self.ac["probe_period_d"]
        ok = row.get("classification") == "orbit" and lo <= p_d <= hi
        info = dict(classification=row.get("classification"), period_d=round(p_d, 1),
                    mass_mj=row.get("mass_mj"), dlnl2_orbit=row.get("dlnl2_orbit"), n_obs=len(ep.w),
                    model=model)
        return ok, info

    def tracked_changes(self):
        r = subprocess.run(["git", "status", "--porcelain", "--untracked-files=no"], cwd=ROOT,
                           capture_output=True, text=True)
        return r.stdout.strip()

    # ---- stages -----------------------------------------------------------------------------
    def tick(self):
        (self.dir / "heartbeat.txt").write_text(f"{stamp()} stage={self.s['stage']} pid={os.getpid()}\n")
        stage = self.s["stage"]
        if stage in ("done", "blocked"):
            return
        start = dt.datetime.fromisoformat(self.ac["window_start_utc"].replace("Z", "+00:00"))
        if utcnow() < start and not self.args.force:
            return                                    # before the window: heartbeat only
        self.keep_awake()
        if stage == "polling":
            self.poll()
        if self.s["stage"] == "probing":
            self.probe()
        if self.s["stage"] == "main_running":
            self.watch_main()
        if self.s["stage"] == "post_running":
            self.watch_post()
        self.save()

    def poll(self):
        tables, err = self.tap_tables(self.args.schema or self.ac["schema"])
        if not tables:
            self.s["polls"] = self.s.get("polls", 0) + 1
            self.log(f"no {self.args.schema or self.ac['schema']} tables yet "
                     f"({err or 'schema empty'}); poll #{self.s['polls']}")
            return
        (self.dir / "tap_tables.txt").write_text("\n".join(tables) + "\n")
        self.s["n_tap_tables"] = len(tables)
        self.log(f"DR4 visible in TAP: {len(tables)} tables (results/dr4_watch/tap_tables.txt)",
                 stage="probing", notify=True)
        self.check_config_tables(tables)

    def probe(self):
        table, rel, errors = self.fetch_probe()
        if table is None:
            n = self.s["probe_failures"] = self.s.get("probe_failures", 0) + 1
            self.log(f"DataLink probe failed (#{n}): {' | '.join(errors)[:600]}")
            if n == int(self.ac.get("probe_failures_before_alert", 12)):
                self.attention(f"DataLink probe has failed {n} times; check data_release / "
                               f"release_candidates and the archive status")
            return
        cols = list(table.colnames)
        try:
            ok, info = self.sanity(table)
        except Exception as err:
            ok, info = False, dict(error=f"{type(err).__name__}: {err}")
        self.s["probe"] = dict(release=rel, columns=cols, **info)
        if rel not in ("file", self.cfg["data_release"]):
            self.attention(f"DataLink worked with data_release={rel!r}, not the configured "
                           f"{self.cfg['data_release']!r}: set it in config/release_day.json")
            self.log("probe OK only with another release string; blocked until config is fixed",
                     stage="blocked", notify=True)
            return
        if not ok:
            self.log(f"Gaia-4 sanity check FAILED: {info}", stage="blocked", notify=True)
            return
        changed = self.tracked_changes()
        if changed and not self.args.allow_dirty:
            self.log(f"tracked files modified; commit before the main run: {changed[:300]}",
                     stage="blocked", notify=True)
            return
        self.log(f"Gaia-4 sanity OK: {info}", notify=True)
        self.launch_main()

    def main_cmd(self):
        extra = shlex.split(self.args.main_args) if self.args.main_args else ["--source", "datalink"]
        return [PY, "scripts/release_day.py", *extra, "--out", str(self.out),
                "--procs", str(self.procs), "--raw-dir", str(ROOT / self.ac["raw_dir"])]

    def launch_main(self):
        self.s["main_launches"] = self.s.get("main_launches", 0) + 1
        self.s["main_started"] = stamp()
        self.s["main_pid"] = self.spawn(self.main_cmd(), "release_day.log")
        self.log(f"main run launched (#{self.s['main_launches']}, pid {self.s['main_pid']}): "
                 f"{shlex.join(self.main_cmd())}", stage="main_running", notify=True)

    def watch_main(self):
        if alive(self.s.get("main_pid")):
            st = self.out / "status.json"
            if st.exists():
                try:
                    d = json.loads(st.read_text())
                    self.s["message"] = (f"main run: stage {d.get('stage')}, fetched "
                                         f"{d.get('n_fetched', '?')}, fitted {d.get('n_done', '?')}/"
                                         f"{d.get('n_jobs', '?')}")
                except json.JSONDecodeError:
                    pass
            return
        summ_path = self.out / "run_summary.json"
        started = dt.datetime.fromisoformat(self.s["main_started"].replace("Z", "+00:00")).timestamp()
        summ = (json.loads(summ_path.read_text())
                if summ_path.exists() and summ_path.stat().st_mtime >= started else None)
        retry = self.s["main_launches"] < int(self.ac.get("max_launches", 6))
        if summ is None:
            if retry:
                self.log("main run ended without a run_summary.json; relaunching (resumes)")
                self.launch_main()
            else:
                self.log("main run keeps dying; see results/dr4_watch/release_day.log",
                         stage="blocked", notify=True)
            return
        frac = summ["n_missing"] / max(summ["n_hosts"], 1)
        self.s["main_summary"] = {k: summ[k] for k in ("n_hosts", "n_fitted", "n_missing", "n_errors",
                                                       "classifications")}
        if frac > float(self.ac.get("max_missing_frac", 0.05)) and retry:
            self.log(f"main run finished with {summ['n_missing']} hosts missing data "
                     f"({frac:.0%}); relaunching to fetch the rest")
            self.launch_main()
            return
        msg = (f"main run done: {summ['n_fitted']} fitted, {summ['n_missing']} missing, "
               f"{summ['n_errors']} errors, {summ['classifications']}")
        if summ["n_missing"]:
            self.attention(f"{summ['n_missing']} hosts still have no epoch data after "
                           f"{self.s['main_launches']} launches")
        self.log(msg, notify=True)
        cmd = [PY, "scripts/dr4_watch.py", "--post", "--config", str(self.args.config),
               "--procs", str(self.procs)]
        if self.args.main_args:
            cmd += ["--main-args", self.args.main_args]
        self.s["post_pid"] = self.spawn(cmd, "post.log")
        self.log(f"post-run chain launched (pid {self.s['post_pid']})", stage="post_running")

    def watch_post(self):
        if alive(self.s.get("post_pid")):
            return
        p = self.dir / "post_status.json"
        post = json.loads(p.read_text()) if p.exists() else {}
        self.s["post"] = post
        bad = [k for k, v in post.items() if v != 0]
        if bad or not post:
            self.attention(f"post-run steps failed: {bad or 'no post_status.json'} (results/dr4_watch/post.log)")
        self.log("automatic part finished; vetting is next (docs/RELEASE_DAY.md)", stage="done",
                 notify=True)

    def post(self):
        """Run after the main run, sequentially; record each step's return code."""
        res = str(self.out / "results.csv")
        steps = {
            "calibrate_null_real": [PY, "scripts/calibrate_null_real.py", "--k", "5",
                                    "--procs", str(self.procs), "--config", str(self.args.config)],
            "population_analysis": [PY, "scripts/population_analysis.py", res, "--out", str(self.out)],
            "make_table": [PY, "scripts/make_table.py", res, "--out", str(self.out / "table")],
            "dashboard": [PY, "scripts/dashboard.py", str(self.out)],
        }
        if self.args.main_args and "simulate" in self.args.main_args:
            steps.pop("calibrate_null_real")         # needs real control-star data
        status = {}
        for name, cmd in steps.items():
            print(f"===== {stamp()} {name}: {shlex.join(cmd)}", flush=True)
            status[name] = subprocess.run(cmd, cwd=ROOT).returncode
            (self.dir / "post_status.json").write_text(json.dumps(status, indent=1))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", type=Path, default=CONFIG)
    ap.add_argument("--status", action="store_true", help="print state.json and exit")
    ap.add_argument("--clear", action="store_true", help="blocked -> probing (after a fix)")
    ap.add_argument("--post", action="store_true", help=argparse.SUPPRESS)
    ap.add_argument("--procs", type=int, default=None)
    ap.add_argument("--quiet", action="store_true", help="no macOS notifications")
    ap.add_argument("--no-caffeinate", action="store_true")
    # testing only
    ap.add_argument("--force", action="store_true", help="ignore window_start_utc")
    ap.add_argument("--schema", default=None, help="test: poll this schema instead (e.g. gaiadr3)")
    ap.add_argument("--probe-file", type=Path, default=None, help="test: Gaia-4 VOTable instead of DataLink")
    ap.add_argument("--probe-model", default=None, help="test: cascade model for the probe")
    ap.add_argument("--main-args", default=None, help="test: release_day.py args instead of --source datalink")
    ap.add_argument("--allow-dirty", action="store_true")
    args = ap.parse_args(argv)
    w = Watch(args)
    if args.status:
        print(json.dumps({k: v for k, v in w.s.items() if k != "history"}, indent=1, default=str))
        for h in w.s["history"][-10:]:
            print(f"{h['t']} [{h['stage']}] {h['msg']}")
        return
    if args.clear:
        w.s.pop("attention", None)
        w.s["probe_failures"] = 0
        w.log("cleared by hand", stage="probing" if w.s["stage"] == "blocked" else w.s["stage"])
        w.save()
        return
    if args.post:
        w.post()
        return
    w.tick()


if __name__ == "__main__":
    main()
