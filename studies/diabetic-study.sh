#!/bin/bash
# Diabetic wound treatment study
#
# Runs a diabetic wound baseline + therapeutic interventions, then generates
# an Excel workbook comparing outcomes.  Pre-generated results ship as
# example.xlsx so you can inspect the output without running the sim.
#
# Usage:
#   ./studies/diabetic-study.sh                        # all 8 treatments (~3-4h)
#   ./studies/diabetic-study.sh --quick                # baseline + 2 (~45min)
#   ./studies/diabetic-study.sh --treatments=hbo,npwt  # pick specific ones
#   ./studies/diabetic-study.sh --combos=pairs         # singles + all pairwise combos
#   ./studies/diabetic-study.sh --combos=triples       # singles + 2- and 3-combos
#   ./studies/diabetic-study.sh --combos=all           # singles + all 2^N-1 subsets
#   ./studies/diabetic-study.sh --no-excel             # skip Excel generation
#
# Output lands in studies/diabetic-results/

echo "Loading..."

STUDY="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$STUDY/.." && pwd)"
cd "$REPO"

RESULTS="$STUDY/diabetic-results"
TREATMENTS="all"
COMBOS=""
EXCEL=true

# --- Parse arguments ---
for arg in "$@"; do
  case "$arg" in
    --quick)        TREATMENTS="growth_factor,npwt" ;;
    --treatments=*) TREATMENTS="${arg#*=}" ;;
    --combos=*)     COMBOS="${arg#*=}" ;;
    --no-excel)     EXCEL=false ;;
    --help|-h)
      sed -n '2,17p' "$0" | sed 's/^# \?//'
      echo ""
      echo "Available treatments:"
      for f in treatments/*.toml; do
        name=$(basename "$f" .toml)
        desc=$(head -1 "$f" | sed 's/^# Treatment: //')
        printf "  %-20s %s\n" "$name" "$desc"
      done
      exit 0 ;;
    *) echo "Unknown argument: $arg"; exit 1 ;;
  esac
done

# --- Source BioDynaMo (before strict mode -- thisbdm.sh isn't compatible) ---
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

# --- Count treatments ---
if [ "$TREATMENTS" = "all" ]; then
  n_treat=$(ls treatments/*.toml | wc -l)
else
  n_treat=$(echo "$TREATMENTS" | tr ',' '\n' | wc -l)
fi

echo ""
echo "============================================"
echo "  Diabetic Treatment Study"
echo "  ${n_treat} treatments (+ baseline)"
if [ -n "$COMBOS" ]; then
  echo "  Combinations: ${COMBOS}"
fi
echo "============================================"
echo ""

# --- Run treatment study ---
t0=$SECONDS
COMBO_FLAG=""
if [ -n "$COMBOS" ]; then
  COMBO_FLAG="--combos=${COMBOS}"
fi
python3 -u scripts/study/treatment_study.py --treatments="${TREATMENTS}" ${COMBO_FLAG}
elapsed=$(( SECONDS - t0 ))
mins=$(( elapsed / 60 ))
secs=$(( elapsed % 60 ))
echo ""
echo "Study completed in ${mins}m ${secs}s"

# --- Collect results into study dir ---
SRC="output/treatment_study"
mkdir -p "$RESULTS"
cp "$SRC"/*.csv "$RESULTS"/ 2>/dev/null || true

# --- Generate Excel workbook ---
if $EXCEL; then
  if python3 -c "import openpyxl" 2>/dev/null; then
    echo ""
    echo "Generating Excel workbook..."
    python3 scripts/study/generate_study.py "$RESULTS" -o "$RESULTS" --name "diabetic_study" 2>&1 || true
  else
    echo ""
    echo "Skipping Excel (pip install openpyxl to enable)"
  fi
fi

# --- Generate plots for baseline ---
if [ -f "$RESULTS/metrics_baseline.csv" ]; then
  echo ""
  echo "Generating baseline plots..."
  PLOT_DIR="$RESULTS/plots"
  mkdir -p "$PLOT_DIR"
  python3 scripts/analysis/plot_metrics.py "$RESULTS/metrics_baseline.csv" -o "$PLOT_DIR" 2>&1 || true
fi

# --- Resilience analysis ---
if [ -f scripts/analysis/resilience_analysis.py ]; then
  echo ""
  echo "Running resilience analysis..."
  python3 scripts/analysis/resilience_analysis.py "$RESULTS" 2>&1 || true
fi

# --- Summary ---
n_csv=$(ls "$RESULTS"/*.csv 2>/dev/null | wc -l)
n_xlsx=$(ls "$RESULTS"/*.xlsx 2>/dev/null | wc -l)
n_plots=$(ls "$RESULTS"/plots/*.png 2>/dev/null | wc -l)

echo ""
echo "============================================"
echo "  Study complete"
echo "============================================"
echo "  Results in: studies/diabetic-results/"
echo "  CSVs:  ${n_csv}"
echo "  Excel: ${n_xlsx}"
echo "  Plots: ${n_plots}"
echo ""
ls -1 "$RESULTS"/ 2>/dev/null | sed 's/^/    /'
echo ""

# Open folder if desktop is available
if [ -n "${DISPLAY:-}" ]; then
  xdg-open "$RESULTS" 2>/dev/null &
fi
