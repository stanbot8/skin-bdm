#!/usr/bin/env bash
# A/B performance comparison: BDM master vs dev/stanbot8
# Builds BDM from each branch, rebuilds skibidy, runs batch timing.
#
# Usage:
#   ./scripts/perf_ab.sh              # defaults: master vs dev/stanbot8, 5 runs
#   ./scripts/perf_ab.sh 10           # 10 runs per branch
#   ./scripts/perf_ab.sh 5 master feat/my-branch
#
# Output: scripts/perf_ab_results.txt

set -euo pipefail

RUNS="${1:-5}"
BRANCH_A="${2:-master}"
BRANCH_B="${3:-dev/stanbot8}"
BDM_DIR="$HOME/biodynamo"
SKI_DIR="$HOME/skibidy"
RESULTS="$SKI_DIR/scripts/perf_ab_results.txt"

log() { echo "[perf_ab] $*"; }

build_bdm() {
    local branch="$1"
    log "Checking out BDM branch: $branch"
    cd "$BDM_DIR"
    git checkout "$branch" --quiet
    git submodule update --init --recursive --quiet 2>/dev/null || true

    log "Building BDM ($branch)..."
    cd build
    cmake -DCMAKE_BUILD_TYPE=Release "$BDM_DIR" > /dev/null 2>&1
    make -j"$(nproc)" 2>&1 | tail -1
    cd "$SKI_DIR"
}

build_skibidy() {
    log "Rebuilding skibidy..."
    cd "$SKI_DIR"
    source "$BDM_DIR/build/bin/thisbdm.sh" > /dev/null 2>&1
    rm -rf build
    mkdir build && cd build
    cmake -DCMAKE_BUILD_TYPE=Release .. > /dev/null 2>&1
    make -j"$(nproc)" 2>&1 | tail -1
    cd "$SKI_DIR"
}

run_batch() {
    local label="$1"
    local n="$2"
    log "Running $n batch simulations ($label)..."
    cd "$SKI_DIR"
    source "$BDM_DIR/build/bin/thisbdm.sh" > /dev/null 2>&1

    local times=()
    for i in $(seq 1 "$n"); do
        local start
        start=$(date +%s%N)
        python3 batch/batch.py -n 1 --skin normal --study wound > /dev/null 2>&1
        local end
        end=$(date +%s%N)
        local elapsed
        elapsed=$(echo "scale=2; ($end - $start) / 1000000000" | bc)
        times+=("$elapsed")
        log "  Run $i/$n: ${elapsed}s"
    done

    # Compute mean and stddev
    local sum=0
    for t in "${times[@]}"; do
        sum=$(echo "$sum + $t" | bc)
    done
    local mean
    mean=$(echo "scale=2; $sum / $n" | bc)

    local sumsq=0
    for t in "${times[@]}"; do
        local diff
        diff=$(echo "$t - $mean" | bc)
        sumsq=$(echo "$sumsq + $diff * $diff" | bc)
    done
    local stddev
    stddev=$(echo "scale=2; sqrt($sumsq / $n)" | bc)

    echo "$label: ${mean}s +/- ${stddev}s (n=$n, times: ${times[*]})"
}

collect_perf_stat() {
    local label="$1"
    log "Collecting perf stat ($label)..."
    cd "$SKI_DIR"
    source "$BDM_DIR/build/bin/thisbdm.sh" > /dev/null 2>&1

    # Single run with perf stat
    python3 -c "
import subprocess, os, sys, shutil
sys.path.insert(0, os.path.join('batch', os.pardir))
from batch import lib
lib.apply_study('wound')
lib.override_param('skin.headless', True)
lib.override_param('visualization.export', False)
lib.merge_config()
" 2>/dev/null

    perf stat -e task-clock,cycles,instructions,cache-references,cache-misses,branches,branch-misses \
        ./build/skibidy 2>&1 | grep -E '(task-clock|cycles|instructions|cache-ref|cache-miss|branch)' || \
        log "  perf stat not available (try: sudo sysctl kernel.perf_event_paranoid=1)"
}

# Warmup run (prime filesystem cache)
warmup() {
    log "Warmup run..."
    cd "$SKI_DIR"
    source "$BDM_DIR/build/bin/thisbdm.sh" > /dev/null 2>&1
    python3 batch/batch.py -n 1 --skin normal --study wound > /dev/null 2>&1
}

# --- Main ---

log "A/B Performance Comparison"
log "Branch A: $BRANCH_A"
log "Branch B: $BRANCH_B"
log "Runs per branch: $RUNS"
log "BDM dir: $BDM_DIR"
log "Skibidy dir: $SKI_DIR"
echo

{
    echo "=== A/B Performance Comparison ==="
    echo "Date: $(date)"
    echo "Branch A: $BRANCH_A"
    echo "Branch B: $BRANCH_B"
    echo "Runs: $RUNS"
    echo "CPU: $(lscpu | grep 'Model name' | sed 's/.*: *//')"
    echo "Cores: $(nproc)"
    echo

    # --- Branch A ---
    build_bdm "$BRANCH_A"
    build_skibidy
    warmup
    echo "--- Branch A ($BRANCH_A) ---"
    run_batch "$BRANCH_A" "$RUNS"
    echo
    echo "perf stat ($BRANCH_A):"
    collect_perf_stat "$BRANCH_A" 2>&1 || true
    echo

    # --- Branch B ---
    build_bdm "$BRANCH_B"
    build_skibidy
    warmup
    echo "--- Branch B ($BRANCH_B) ---"
    run_batch "$BRANCH_B" "$RUNS"
    echo
    echo "perf stat ($BRANCH_B):"
    collect_perf_stat "$BRANCH_B" 2>&1 || true
    echo

    # Restore dev branch
    cd "$BDM_DIR" && git checkout "$BRANCH_B" --quiet
    echo "=== Done ==="

} 2>&1 | tee "$RESULTS"

log "Results saved to $RESULTS"
