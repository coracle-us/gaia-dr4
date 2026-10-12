"""Build (or keep rebuilding) the HTML dashboard for a release-day run directory.

    python scripts/dashboard.py results/dr4                 # once
    python scripts/dashboard.py results/dr4 --watch 60      # every 60 s until the run is done
    python scripts/dashboard.py results/dr4 --manual results/dr4/manual_vetting.csv

Open results/dr4/dashboard/index.html in a browser; it reloads itself while the run is going.
scripts/release_day.py already rebuilds the page as hosts finish; use this to rebuild after
editing manual vetting decisions or after scripts/population_analysis.py has run.
Manual vetting CSV columns: hostname, check, status, note (checks: see companions.dashboard.CHECKS).
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

from companions.dashboard import build_dashboard


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out", type=Path, help="release_day.py --out directory")
    ap.add_argument("--manual", type=Path, default=None,
                    help="manual vetting CSV (default OUT/manual_vetting.csv)")
    ap.add_argument("--watch", type=float, default=0, help="rebuild every N seconds until done")
    args = ap.parse_args(argv)
    while True:
        path = build_dashboard(args.out, manual_path=args.manual)
        print(f"{time.strftime('%H:%M:%S')} wrote {path}")
        status = args.out / "status.json"
        stage = json.loads(status.read_text()).get("stage") if status.exists() else "done"
        if not args.watch or stage == "done":
            break
        time.sleep(args.watch)


if __name__ == "__main__":
    main()
