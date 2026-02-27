#!/bin/bash
# Adaptive surrogate-guided treatment combination search
#
# Efficiently explores 2^N-1 treatment combinations using adaptive consensus
# (early stopping when variance is low), an additive surrogate model (predict
# combos from singles), and guided selection (only run promising combos).
#
# Usage:
#   ./studies/adaptive-study.sh                      # full search (~30-45min)
#   ./studies/adaptive-study.sh --quick              # top-5, 2-3 runs (~10min)
#   ./studies/adaptive-study.sh --top-k=10           # top-10 combos
#   ./studies/adaptive-study.sh --treatments=hbo,npwt,growth_factor
#
# Output lands in studies/adaptive-results/

echo "Loading..."

STUDY="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$STUDY/.." && pwd)"
cd "$REPO"

RESULTS="$STUDY/adaptive-results"
TOP_K=20
MIN_RUNS=3
MAX_RUNS=10
CV_THRESH=0.05
TREATMENTS="all"
EXCEL=true

# --- Parse arguments ---
for arg in "$@"; do
  case "$arg" in
    --quick)
      TOP_K=5; MIN_RUNS=2; MAX_RUNS=3 ;;
    --top-k=*)    TOP_K="${arg#*=}" ;;
    --min-runs=*) MIN_RUNS="${arg#*=}" ;;
    --max-runs=*) MAX_RUNS="${arg#*=}" ;;
    --cv=*)       CV_THRESH="${arg#*=}" ;;
    --treatments=*) TREATMENTS="${arg#*=}" ;;
    --no-excel)   EXCEL=false ;;
    --help|-h)
      sed -n '2,14p' "$0" | sed 's/^# \?//'
      echo ""
      echo "Available treatments:"
      for f in treatments/*.toml; do
        name=$(basename "$f" .toml)
        [ "$name" = "combination" ] && continue
        desc=$(head -1 "$f" | sed 's/^# Treatment: //')
        printf "  %-20s %s\n" "$name" "$desc"
      done
      exit 0 ;;
    *) echo "Unknown argument: $arg"; exit 1 ;;
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
  cd build && cmake .. -DCMAKE_BUILD_TYPE=Release > /dev/null 2>&1
  make -j"$(nproc)" 2>&1 | tail -3
  if [ ${PIPESTATUS[0]} -ne 0 ]; then
    echo "Build failed."
    exit 1
  fi
  cd "$REPO"
else
  echo "Binary up to date."
fi

echo ""
echo "============================================"
echo "  Adaptive Combination Search"
echo "  Top-K: ${TOP_K}, Runs: ${MIN_RUNS}-${MAX_RUNS}"
echo "  CV threshold: ${CV_THRESH}"
echo "============================================"
echo ""

# --- Run adaptive study ---
t0=$SECONDS
python3 -u scripts/study/adaptive_study.py \
  --treatments="${TREATMENTS}" \
  --top-k="${TOP_K}" \
  --min-runs="${MIN_RUNS}" \
  --max-runs="${MAX_RUNS}" \
  --cv-threshold="${CV_THRESH}"
elapsed=$(( SECONDS - t0 ))
mins=$(( elapsed / 60 ))
secs=$(( elapsed % 60 ))
echo ""
echo "Study completed in ${mins}m ${secs}s"

# --- Collect results ---
SRC="output/adaptive_study"
mkdir -p "$RESULTS"
cp "$SRC"/*.csv "$RESULTS"/ 2>/dev/null || true

# --- Generate Excel workbook ---
if $EXCEL; then
  if python3 -c "import openpyxl" 2>/dev/null; then
    echo ""
    echo "Generating Excel workbook..."
    python3 scripts/study/generate_study.py "$SRC" -o "$RESULTS" --name "adaptive_study" 2>&1 || true
  else
    echo ""
    echo "Skipping Excel (pip install openpyxl to enable)"
  fi
fi

# --- Summary ---
n_csv=$(ls "$RESULTS"/*.csv 2>/dev/null | wc -l)
n_xlsx=$(ls "$RESULTS"/*.xlsx 2>/dev/null | wc -l)

echo ""
echo "============================================"
echo "  Adaptive study complete"
echo "============================================"
echo "  Results in: studies/adaptive-results/"
echo "  CSVs:  ${n_csv}"
echo "  Excel: ${n_xlsx}"
echo ""
ls -1 "$RESULTS"/ 2>/dev/null | sed 's/^/    /'
echo ""
