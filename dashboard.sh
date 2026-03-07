#!/usr/bin/env bash
# Launch the SkiBiDy dashboard
cd "$(dirname "$0")" || exit 1
source ~/biodynamo/build/bin/thisbdm.sh
python3 scripts/dashboard.py "$@"
