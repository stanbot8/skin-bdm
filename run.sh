#!/bin/bash
VIEW=true
VALIDATE=true
PRESET=""
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

  # Scenario picker (ordered by complexity)
  echo "Scenario:"
  printf "  0) %-14s no overlay\n" "none"
  presets=("wound" "diabetic_wound" "tumor" "tumor_wound" "nothing")
  for i in "${!presets[@]}"; do
    name="${presets[$i]}"
    desc=$(head -1 "presets/${name}.toml" | sed 's/^# Preset: //')
    printf "  %d) %-14s %s\n" "$((i+1))" "$name" "$desc"
  done
  printf "Pick scenario [1]: "
  read -r preset_choice
  preset_choice="${preset_choice:-1}"
  if [ "$preset_choice" -ge 1 ] && [ "$preset_choice" -le "${#presets[@]}" ] 2>/dev/null; then
    PRESET="${presets[$((preset_choice-1))]}"
  fi
  echo ""

  # Summary
  echo "Running: skin=${SKIN:-none} scenario=${PRESET:-none}"
  echo "---"
else
  for arg in "$@"; do
    case "$arg" in
      --no-view)     VIEW=false ;;
      --no-validate) VALIDATE=false ;;
      --preset=*)    PRESET="${arg#--preset=}" ;;
      --skin=*)      SKIN="${arg#--skin=}" ;;
      --diabetic)    PRESET="diabetic_wound" ;;
      --treatment=*) TREATMENT="${arg#--treatment=}" ;;
      --compare)     COMPARE=true ;;
      --study)       python3 scripts/study/treatment_study.py; exit 0 ;;
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
      --list-presets)
        echo "Available presets:"
        for f in presets/*.toml; do
          name=$(basename "$f" .toml)
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
        echo "  --preset=NAME   apply preset (e.g. wound, tumor, diabetic_wound)"
        echo "  --diabetic      shorthand for --preset=diabetic_wound"
        echo "  --compare       run normal + diabetic back-to-back and compare"
        echo "  --list-skins    show available skin profiles"
        echo "  --list-presets  show available presets"
        echo "  --no-view       skip biodynamo view after simulation"
        echo "  --no-validate   skip validation scripts after simulation"
        exit 0 ;;
    esac
  done

  # Handle --preset <name> and --skin <name> (space-separated) for convenience
  args=("$@")
  for i in "${!args[@]}"; do
    if [ "${args[$i]}" = "--preset" ] && [ -n "${args[$((i+1))]}" ]; then
      PRESET="${args[$((i+1))]}"
    fi
    if [ "${args[$i]}" = "--skin" ] && [ -n "${args[$((i+1))]}" ]; then
      SKIN="${args[$((i+1))]}"
    fi
  done
fi

if [ -z "$BDMSYS" ]; then
  echo "ERROR: BioDynaMo not sourced. Run: source <path>/bin/thisbdm.sh"
  exit 1
fi
cd "$(dirname "$0")"

# --- Compare mode: run normal + diabetic, then diff ---
if [ "$COMPARE" = true ]; then
  echo "=== Compare mode: normal vs diabetic ==="
  rm -rf output
  mkdir -p output/compare

  # Normal run
  echo "--- Running normal wound ---"
  rm -f bdm.toml
  python3 scripts/config/merge_config.py || exit 1
  python3 scripts/config/apply_preset.py presets/wound.toml bdm.toml || exit 1
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
  python3 scripts/config/apply_preset.py presets/diabetic_wound.toml bdm.toml || exit 1
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

# Phase 2: preset (scenario parameters)
if [ -n "$PRESET" ]; then
  PRESET_FILE="presets/${PRESET}.toml"
  if [ ! -f "$PRESET_FILE" ]; then
    echo "ERROR: preset '$PRESET' not found. Use --list-presets to see available presets."
    exit 1
  fi
  python3 scripts/config/apply_preset.py "$PRESET_FILE" bdm.toml || exit 1
fi

# Phase 3: treatment overlay (applied on top of profile + preset)
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
