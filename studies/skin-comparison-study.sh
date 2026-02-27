#!/bin/bash
# Skin type comparison study
#
# Runs wound healing with normal, aged, and diabetic skin profiles using
# consensus batches, then compares healing trajectories.
#
# Usage:
#   ./studies/skin-comparison-study.sh           # all profiles (10-run consensus)
#   ./studies/skin-comparison-study.sh --quick   # 5-run consensus
#   ./studies/skin-comparison-study.sh -n 20     # 20-run consensus
#
# Output lands in studies/skin-comparison-results/

echo "Loading..."

STUDY="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$STUDY/.." && pwd)"
cd "$REPO"

RESULTS="$STUDY/skin-comparison-results"
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
echo "  Skin Type Comparison Study"
echo "  ${N_RUNS}-run consensus per profile"
echo "============================================"
echo ""

mkdir -p "$RESULTS"
t0=$SECONDS

# --- Profile 1: Normal wound ---
echo "[1/3] Normal skin + wound"
python3 -u batch/batch.py -n "$N_RUNS" --skin normal --preset wound --validate 2>&1 | tee "$RESULTS/log_normal.txt"
LATEST=$(ls -td batch/results/consensus_normal_wound_* 2>/dev/null | head -1)
if [ -n "$LATEST" ]; then
  cp "$LATEST"/consensus.csv "$RESULTS/consensus_normal.csv" 2>/dev/null || true
  cp "$LATEST"/metrics.csv "$RESULTS/metrics_normal.csv" 2>/dev/null || true
fi
echo ""

# --- Profile 2: Aged wound ---
echo "[2/3] Aged skin + wound"
python3 -u batch/batch.py -n "$N_RUNS" --skin aged --preset wound --validate 2>&1 | tee "$RESULTS/log_aged.txt"
LATEST=$(ls -td batch/results/consensus_aged_wound_* 2>/dev/null | head -1)
if [ -n "$LATEST" ]; then
  cp "$LATEST"/consensus.csv "$RESULTS/consensus_aged.csv" 2>/dev/null || true
  cp "$LATEST"/metrics.csv "$RESULTS/metrics_aged.csv" 2>/dev/null || true
fi
echo ""

# --- Profile 3: Diabetic wound ---
echo "[3/3] Diabetic skin + wound"
python3 -u batch/batch.py -n "$N_RUNS" --skin diabetic --preset diabetic_wound --validate 2>&1 | tee "$RESULTS/log_diabetic.txt"
LATEST=$(ls -td batch/results/consensus_diabetic_diabetic_wound_* 2>/dev/null | head -1)
if [ -n "$LATEST" ]; then
  cp "$LATEST"/consensus.csv "$RESULTS/consensus_diabetic.csv" 2>/dev/null || true
  cp "$LATEST"/metrics.csv "$RESULTS/metrics_diabetic.csv" 2>/dev/null || true
fi
echo ""

# --- Generate plots for each ---
for csv_file in "$RESULTS"/metrics_*.csv; do
  if [ -f "$csv_file" ]; then
    name=$(basename "$csv_file" .csv | sed 's/^metrics_//')
    echo "Plotting: $name"
    PLOT_DIR="$RESULTS/plots"
    mkdir -p "$PLOT_DIR"
    python3 scripts/analysis/plot_metrics.py "$csv_file" -o "$PLOT_DIR" --prefix "${name}_" 2>&1 || true
  fi
done

# --- Cross-profile comparison ---
csvs=""
labels=""
for profile in normal aged diabetic; do
  f="$RESULTS/metrics_${profile}.csv"
  if [ -f "$f" ]; then
    csvs="$csvs $f"
    labels="${labels:+$labels,}$profile"
  fi
done

if [ -n "$csvs" ]; then
  echo ""
  echo "Comparing profiles..."
  python3 scripts/analysis/compare_sims.py $csvs --labels "$labels" 2>&1 || true
fi

elapsed=$(( SECONDS - t0 ))
mins=$(( elapsed / 60 ))
secs=$(( elapsed % 60 ))

# --- Summary ---
n_csv=$(ls "$RESULTS"/*.csv 2>/dev/null | wc -l)
n_plots=$(ls "$RESULTS"/plots/*.png 2>/dev/null | wc -l)

echo ""
echo "============================================"
echo "  Skin comparison complete (${mins}m ${secs}s)"
echo "============================================"
echo "  Results in: studies/skin-comparison-results/"
echo "  CSVs:  ${n_csv}"
echo "  Plots: ${n_plots}"
echo ""
ls -1 "$RESULTS"/ 2>/dev/null | sed 's/^/    /'
echo ""
