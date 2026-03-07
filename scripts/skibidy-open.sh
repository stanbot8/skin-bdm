#!/usr/bin/env bash
# Open a .skibidy project file in the SkiBiDy Dashboard.
# Usage: skibidy-open.sh path/to/study.skibidy
#   or:  double-click a .skibidy file (via desktop file association)

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DASHBOARD="$SCRIPT_DIR/dashboard.py"

if [ ! -f "$DASHBOARD" ]; then
    echo "Error: dashboard.py not found at $DASHBOARD"
    exit 1
fi

exec python3 "$DASHBOARD" "$@"
