#!/bin/bash
VIEW=true
VALIDATE=true
STUDY=""
SKIN=""
TREATMENT=""
COMPARE=false

# --- Interactive menu (no arguments) ---
if [ $# -eq 0 ]; then
  echo "=== skibidy ==="
  echo ""

  # Skin profile picker (normal first, then alphabetical)
  echo "Skin profile:"
  printf "  0) %-12s no overlay\n" "none"
  skins=("normal")
  for f in profiles/*.toml; do
    name=$(basename "$f" .toml)
    [ "$name" = "TEMPLATE" ] || [ "$name" = "normal" ] && continue
    skins+=("$name")
  done
  for i in "${!skins[@]}"; do
    name="${skins[$i]}"
    desc=$(head -1 "profiles/${name}.toml" | sed 's/^# Profile: //')
    printf "  %d) %-12s %s\n" "$((i+1))" "$name" "$desc"
  done
  printf "Pick skin [1]: "
  read -r skin_choice
  skin_choice="${skin_choice:-1}"
  if [ "$skin_choice" -ge 1 ] && [ "$skin_choice" -le "${#skins[@]}" ] 2>/dev/null; then
    SKIN="${skins[$((skin_choice-1))]}"
  fi
  echo ""

  # Study picker (ordered by complexity)
  echo "Study:"
  printf "  0) %-14s no overlay\n" "none"
  studies=("wound" "diabetic-wound" "tumor" "tumor-wound" "baseline")
  for i in "${!studies[@]}"; do
    name="${studies[$i]}"
    desc=$(head -1 "studies/${name}/preset.toml" | sed 's/^# Preset: //')
    printf "  %d) %-14s %s\n" "$((i+1))" "$name" "$desc"
  done
  printf "Pick study [1]: "
  read -r study_choice
  study_choice="${study_choice:-1}"
  if [ "$study_choice" -ge 1 ] && [ "$study_choice" -le "${#studies[@]}" ] 2>/dev/null; then
    STUDY="${studies[$((study_choice-1))]}"
  fi
  echo ""

  # Summary
  echo "Running: skin=${SKIN:-none} study=${STUDY:-none}"
  echo "---"
else
  for arg in "$@"; do
    case "$arg" in
      --no-view)     VIEW=false ;;
      --no-validate) VALIDATE=false ;;
      --study=*)     STUDY="${arg#--study=}" ;;
      --skin=*)      SKIN="${arg#--skin=}" ;;
      --diabetic)    STUDY="diabetic-wound" ;;
      --treatment=*) TREATMENT="${arg#--treatment=}" ;;
      --compare)     COMPARE=true ;;
      --list-treatments)
        echo "Available treatments:"
        for f in treatments/*.toml; do
          name=$(basename "$f" .toml)
          desc=$(head -1 "$f" | sed 's/^# Treatment: //')
          printf "  %-20s %s\n" "$name" "$desc"
        done
        exit 0 ;;
      --list-skins)
        echo "Available skin profiles:"
        for f in profiles/*.toml; do
          name=$(basename "$f" .toml)
          [ "$name" = "TEMPLATE" ] && continue
          desc=$(head -1 "$f" | sed 's/^# Profile: //')
          keys=$(grep -c '^\s*[a-z].*=' "$f" 2>/dev/null || echo 0)
          printf "  %-14s %s (%d overrides)\n" "$name" "$desc" "$keys"
        done
        echo ""
        echo "Create custom profiles from profiles/TEMPLATE.toml"
        exit 0 ;;
      --list-studies)
        echo "Available studies:"
        for d in studies/*/; do
          name=$(basename "$d")
          f="${d}preset.toml"
          [ -f "$f" ] || continue
          desc=$(head -1 "$f" | sed 's/^# Preset: //')
          printf "  %-14s %s\n" "$name" "$desc"
        done
        exit 0 ;;
      --help|-h)
        echo "Usage: ./run.sh [OPTIONS]"
        echo "       ./run.sh                (interactive menu)"
        echo ""
        echo "Options:"
        echo "  --skin=NAME     apply skin profile (e.g. normal, aged)"
        echo "  --study=NAME    apply study config (e.g. wound, diabetic-wound, tumor)"
        echo "  --diabetic      shorthand for --study=diabetic-wound"
        echo "  --compare       run normal + diabetic-wound back-to-back and compare"
        echo "  --list-skins    show available skin profiles"
        echo "  --list-studies  show available studies"
        echo "  --no-view       skip biodynamo view after simulation"
        echo "  --no-validate   skip validation scripts after simulation"
        exit 0 ;;
    esac
  done

  # Handle --study <name> and --skin <name> (space-separated) for convenience
  args=("$@")
  for i in "${!args[@]}"; do
    if [ "${args[$i]}" = "--study" ] && [ -n "${args[$((i+1))]}" ]; then
      STUDY="${args[$((i+1))]}"
    fi
    if [ "${args[$i]}" = "--skin" ] && [ -n "${args[$((i+1))]}" ]; then
      SKIN="${args[$((i+1))]}"
    fi
  done
fi

if [ -z "$BDMSYS" ]; then
  echo "ERROR: BioDynaMo not sourced. Run: source ~/biodynamo/build/bin/thisbdm.sh"
  exit 1
fi

# Detect BDM version change and force rebuild
BDM_STAMP="build/.bdm_source"
CURRENT_BDM="$(cd "$BDMSYS/.." && git rev-parse HEAD 2>/dev/null || echo unknown)"
if [ -f "$BDM_STAMP" ] && [ "$(cat "$BDM_STAMP")" != "$CURRENT_BDM" ]; then
  echo "BDM source changed, forcing rebuild..."
  rm -rf build
fi
cd "$(dirname "$0")"

# --- Compare mode: run normal + diabetic-wound, then diff ---
if [ "$COMPARE" = true ]; then
  echo "=== Compare mode: normal vs diabetic-wound ==="
  rm -rf output
  mkdir -p output/compare

  # Normal run
  echo "--- Running normal wound ---"
  rm -f bdm.toml
  python3 scripts/config/merge_config.py || exit 1
  python3 scripts/config/apply_preset.py studies/wound/preset.toml bdm.toml || exit 1
  python3 scripts/config/apply_preset.py profiles/normal.toml bdm.toml || exit 1
  ./build/skibidy 2>&1 | grep -v 'pvbatch\|SIGABRT\|Aborted\|ParaviewAdaptor\|Error in'
  cp output/skibidy/metrics.csv output/compare/metrics_normal.csv
  python3 scripts/analysis/plot_metrics.py || true
  python3 literature/validate_all.py output/skibidy/metrics.csv || true
  mkdir -p output/compare/plots_normal
  cp output/plots/*.png output/compare/plots_normal/ 2>/dev/null
  rm -rf output/skibidy output/plots

  # Diabetic run
  echo "--- Running diabetic wound ---"
  rm -f bdm.toml
  python3 scripts/config/merge_config.py || exit 1
  python3 scripts/config/apply_preset.py studies/diabetic-wound/preset.toml bdm.toml || exit 1
  ./build/skibidy 2>&1 | grep -v 'pvbatch\|SIGABRT\|Aborted\|ParaviewAdaptor\|Error in'
  cp output/skibidy/metrics.csv output/compare/metrics_diabetic.csv
  python3 scripts/analysis/plot_metrics.py || true
  python3 literature/validate_all.py output/skibidy/metrics.csv || true
  mkdir -p output/compare/plots_diabetic
  cp output/plots/*.png output/compare/plots_diabetic/ 2>/dev/null

  # Compare
  echo ""
  python3 scripts/analysis/compare_sims.py output/compare/metrics_normal.csv output/compare/metrics_diabetic.csv
  exit 0
fi

# --- Single run ---

# Build runtime config: merge bdm.core.toml + modules/ -> bdm.toml
# bdm.core.toml is committed to git; bdm.toml is gitignored and regenerated.
rm -f bdm.toml
python3 scripts/config/merge_config.py || exit 1

# Phase 1: skin profile (biology parameters)
if [ -n "$SKIN" ]; then
  SKIN_FILE="profiles/${SKIN}.toml"
  if [ ! -f "$SKIN_FILE" ]; then
    echo "ERROR: skin profile '$SKIN' not found. Use --list-skins to see available profiles."
    exit 1
  fi
  python3 scripts/config/apply_preset.py "$SKIN_FILE" bdm.toml || exit 1
fi

# Phase 2: study config (scenario parameters)
if [ -n "$STUDY" ]; then
  STUDY_FILE="studies/${STUDY}/preset.toml"
  if [ ! -f "$STUDY_FILE" ]; then
    echo "ERROR: study '$STUDY' not found. Use --list-studies to see available studies."
    exit 1
  fi
  python3 scripts/config/apply_preset.py "$STUDY_FILE" bdm.toml || exit 1
fi

# Phase 3: treatment overlay (applied on top of profile + study config)
if [ -n "$TREATMENT" ]; then
  TREATMENT_FILE="treatments/${TREATMENT}.toml"
  if [ ! -f "$TREATMENT_FILE" ]; then
    echo "ERROR: treatment '$TREATMENT' not found. Use --list-treatments to see available treatments."
    exit 1
  fi
  python3 scripts/config/apply_preset.py "$TREATMENT_FILE" bdm.toml || exit 1
fi

# Build (skip if binary is newer than all sources)
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
  echo "$CURRENT_BDM" > "$BDM_STAMP"
  cd ..
else
  echo "Binary up to date, skipping build."
fi

# Run simulation directly (no biodynamo run wrapper)
rm -rf output
./build/skibidy || exit 1

# Post-processing: patch pvsm only if pvbatch succeeded (has DISPLAY)
if [ -n "$DISPLAY" ] && [ -f output/skibidy/*.pvsm ] 2>/dev/null; then
  python3 scripts/viz/patch_pvsm.py
fi

# Plot metrics
python3 scripts/analysis/plot_metrics.py || true

# Validation
if [ "$VALIDATE" = true ]; then
  CSV="output/skibidy/metrics.csv"
  if [ -f "$CSV" ]; then
    echo "=== Running validation ==="
    python3 literature/validate_all.py "$CSV"
    echo "=== Validation complete ==="
  else
    echo "Warning: $CSV not found, skipping validation."
  fi
fi

# View (only if requested and DISPLAY is available)
if [ "$VIEW" = true ] && [ -n "$DISPLAY" ]; then
  biodynamo view || true
fi
