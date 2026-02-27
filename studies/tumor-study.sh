#!/bin/bash
# Tumor growth study
#
# Runs tumor and tumor+wound scenarios with consensus batches, then generates
# comparison plots and validation.
#
# Usage:
#   ./studies/tumor-study.sh              # all scenarios (10-run consensus each)
#   ./studies/tumor-study.sh --quick      # 5-run consensus
#   ./studies/tumor-study.sh -n 20        # 20-run consensus
#
# Output lands in studies/tumor-results/

echo "Loading..."

STUDY="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$STUDY/.." && pwd)"
cd "$REPO"

RESULTS="$STUDY/tumor-results"
N_RUNS=10

# --- Parse arguments ---
for arg in "$@"; do
  case "$arg" in
    --quick)  N_RUNS=5 ;;
    -n=*)     N_RUNS="${arg#*=}" ;;
    -n)       shift_next=true ;;
    --help|-h)
      sed -n '2,14p' "$0" | sed 's/^# \?//'
      exit 0 ;;
    *)
      if [ "$shift_next" = true ]; then
        N_RUNS="$arg"
        shift_next=false
      else
        echo "Unknown argument: $arg"; exit 1
      fi
      ;;
  esac
done

# --- Source BioDynaMo ---
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

echo ""
echo "============================================"
echo "  Tumor Growth Study"
echo "  ${N_RUNS}-run consensus per scenario"
echo "============================================"
echo ""

mkdir -p "$RESULTS"
t0=$SECONDS

# --- Scenario 1: Tumor only (no wound) ---
echo "[1/3] Tumor only (normal skin)"
python3 -u batch/batch.py -n "$N_RUNS" --skin normal --preset tumor --validate 2>&1 | tee "$RESULTS/log_tumor.txt"
LATEST_TUMOR=$(ls -td batch/results/consensus_normal_tumor_* 2>/dev/null | head -1)
if [ -n "$LATEST_TUMOR" ]; then
  cp "$LATEST_TUMOR"/consensus.csv "$RESULTS/consensus_tumor.csv" 2>/dev/null || true
  cp "$LATEST_TUMOR"/metrics.csv "$RESULTS/metrics_tumor.csv" 2>/dev/null || true
fi
echo ""

# --- Scenario 2: Tumor + wound ---
echo "[2/3] Tumor with wound (normal skin)"
python3 -u batch/batch.py -n "$N_RUNS" --skin normal --preset tumor_wound --validate 2>&1 | tee "$RESULTS/log_tumor_wound.txt"
LATEST_TW=$(ls -td batch/results/consensus_normal_tumor_wound_* 2>/dev/null | head -1)
if [ -n "$LATEST_TW" ]; then
  cp "$LATEST_TW"/consensus.csv "$RESULTS/consensus_tumor_wound.csv" 2>/dev/null || true
  cp "$LATEST_TW"/metrics.csv "$RESULTS/metrics_tumor_wound.csv" 2>/dev/null || true
fi
echo ""

# --- Scenario 3: Aged skin + tumor ---
echo "[3/3] Tumor with aged skin"
python3 -u batch/batch.py -n "$N_RUNS" --skin aged --preset tumor --validate 2>&1 | tee "$RESULTS/log_aged_tumor.txt"
LATEST_AT=$(ls -td batch/results/consensus_aged_tumor_* 2>/dev/null | head -1)
if [ -n "$LATEST_AT" ]; then
  cp "$LATEST_AT"/consensus.csv "$RESULTS/consensus_aged_tumor.csv" 2>/dev/null || true
  cp "$LATEST_AT"/metrics.csv "$RESULTS/metrics_aged_tumor.csv" 2>/dev/null || true
fi
echo ""

# --- Generate plots ---
for csv_file in "$RESULTS"/metrics_*.csv; do
  if [ -f "$csv_file" ]; then
    name=$(basename "$csv_file" .csv | sed 's/^metrics_//')
    echo "Plotting: $name"
    PLOT_DIR="$RESULTS/plots"
    mkdir -p "$PLOT_DIR"
    python3 scripts/analysis/plot_metrics.py "$csv_file" -o "$PLOT_DIR" --prefix "${name}_" 2>&1 || true
  fi
done

# --- Compare scenarios ---
if [ -f "$RESULTS/metrics_tumor.csv" ] && [ -f "$RESULTS/metrics_tumor_wound.csv" ]; then
  echo ""
  echo "Comparing tumor vs tumor+wound..."
  python3 scripts/analysis/compare_sims.py \
    "$RESULTS/metrics_tumor.csv" "$RESULTS/metrics_tumor_wound.csv" \
    --labels "tumor_only,tumor+wound" 2>&1 || true
fi

elapsed=$(( SECONDS - t0 ))
mins=$(( elapsed / 60 ))
secs=$(( elapsed % 60 ))

# --- Summary ---
n_csv=$(ls "$RESULTS"/*.csv 2>/dev/null | wc -l)
n_plots=$(ls "$RESULTS"/plots/*.png 2>/dev/null | wc -l)

echo ""
echo "============================================"
echo "  Tumor study complete (${mins}m ${secs}s)"
echo "============================================"
echo "  Results in: studies/tumor-results/"
echo "  CSVs:  ${n_csv}"
echo "  Plots: ${n_plots}"
echo ""
ls -1 "$RESULTS"/ 2>/dev/null | sed 's/^/    /'
echo ""
