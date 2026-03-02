#!/bin/bash
# Run experiment scenarios
#
# Usage:
#   ./studies/run-scenarios.sh                                              # all scenarios
#   ./studies/run-scenarios.sh studies/diabetic-wound/scenarios/wound_size.toml  # one scenario
#   ./studies/run-scenarios.sh --quick                                      # all with 2 runs
#   ./studies/run-scenarios.sh --runs=3 studies/diabetic-wound/scenarios/biofilm* # specific count
#
# Output lands in output/scenarios/<scenario_name>/

STUDY="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$STUDY/.." && pwd)"
cd "$REPO"

RUNS=""
SCENARIOS=()

# --- Parse arguments ---
for arg in "$@"; do
  case "$arg" in
    --quick)     RUNS="--runs=2" ;;
    --runs=*)    RUNS="$arg" ;;
    --help|-h)
      sed -n '2,9p' "$0" | sed 's/^# \?//'
      echo ""
      echo "Available scenarios:"
      for f in studies/*/scenarios/*.toml; do
        study=$(basename "$(dirname "$(dirname "$f")")")
        name=$(basename "$f" .toml)
        desc=$(head -5 "$f" | grep "^#" | tail -1 | sed 's/^# //')
        printf "  %-12s %-25s %s\n" "$study" "$name" "$desc"
      done
      exit 0 ;;
    *)
      if [ -f "$arg" ]; then
        SCENARIOS+=("$arg")
      else
        echo "Unknown argument or file not found: $arg"
        exit 1
      fi ;;
  esac
done

# Default: all scenario TOMLs across all studies
if [ ${#SCENARIOS[@]} -eq 0 ]; then
  SCENARIOS=(studies/*/scenarios/*.toml)
fi

# --- Check BioDynaMo ---
if [ -z "$BDMSYS" ]; then
  echo "ERROR: BioDynaMo not sourced. Run: source <path>/bin/thisbdm.sh"
  exit 1
fi

set -eo pipefail

# --- Build if needed ---
BINARY=build/skibidy
if [ ! -f "$BINARY" ] || [ -n "$(find src/ CMakeLists.txt -newer "$BINARY" 2>/dev/null | head -1)" ]; then
  echo "Building..."
  mkdir -p build
  cd build && cmake .. -DCMAKE_BUILD_TYPE=Release > /dev/null 2>&1 && make -j"$(nproc)" 2>&1 | tail -3
  cd "$REPO"
  if [ ! -f "$BINARY" ]; then
    echo "Build failed."
    exit 1
  fi
else
  echo "Binary up to date."
fi

# --- Run scenarios ---
n=${#SCENARIOS[@]}
echo ""
echo "============================================"
echo "  Scenario Runner"
echo "  ${n} scenario(s) to run"
echo "============================================"
echo ""

t0=$SECONDS
for scenario in "${SCENARIOS[@]}"; do
  echo ">>> $(basename "$scenario" .toml)"
  python3 -u scripts/study/scenario_runner.py "$scenario" $RUNS
  echo ""
done

elapsed=$(( SECONDS - t0 ))
mins=$(( elapsed / 60 ))
secs=$(( elapsed % 60 ))

echo "============================================"
echo "  All scenarios completed in ${mins}m ${secs}s"
echo "  Results in: studies/scenarios/"
echo "============================================"
