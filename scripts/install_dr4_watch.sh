#!/bin/zsh
# Install (or with --uninstall remove) the launchd agent that runs scripts/dr4_watch.py every 10 min.
set -e
LABEL=us.coracle.gaia-dr4-watch
SRC="${0:A:h:h}/config/launchd/$LABEL.plist"
DST="$HOME/Library/LaunchAgents/$LABEL.plist"
launchctl bootout "gui/$(id -u)/$LABEL" 2>/dev/null || true
if [[ "$1" == "--uninstall" ]]; then rm -f "$DST"; echo "removed $LABEL"; exit 0; fi
plutil -lint "$SRC"
cp "$SRC" "$DST"
launchctl bootstrap "gui/$(id -u)" "$DST"
launchctl print "gui/$(id -u)/$LABEL" | grep -E "state|run interval|last exit" || true
echo "installed $LABEL; heartbeat: results/dr4_watch/heartbeat.txt, log: ~/Library/Logs/gaia-dr4-watch.log"
