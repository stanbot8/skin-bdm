"""Shared utilities for the batch/sweep system.

Provides config merging, parameter override, simulation running,
CSV aggregation, and outcome extraction.
"""

import csv
import math
import os
import re
import shutil
import subprocess
import sys
import time

# Project root (one level up from batch/)
ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), os.pardir))


# ---------------------------------------------------------------------------
# Config management
# ---------------------------------------------------------------------------

def merge_config():
    """Run merge_config.py to produce bdm.toml from bdm.core.toml + modules."""
    rc = subprocess.run(
        ["python3", os.path.join(ROOT, "scripts", "config", "merge_config.py")],
        cwd=ROOT, capture_output=True, text=True)
    if rc.returncode != 0:
        print(f"Config merge failed: {rc.stderr}")
        sys.exit(1)
    return rc.stdout.strip()


def apply_overlay(overlay_path):
    """Apply a TOML overlay (profile or preset) to bdm.toml."""
    rc = subprocess.run(
        ["python3", os.path.join(ROOT, "scripts", "config", "apply_preset.py"),
         overlay_path, "bdm.toml"],
        cwd=ROOT, capture_output=True, text=True)
    if rc.returncode != 0:
        print(f"Overlay failed ({overlay_path}): {rc.stderr}")
        sys.exit(1)


def apply_profile(name):
    """Apply a skin profile by name."""
    path = os.path.join(ROOT, "profiles", f"{name}.toml")
    if not os.path.isfile(path):
        print(f"ERROR: skin profile '{name}' not found.")
        sys.exit(1)
    apply_overlay(path)


def apply_preset(name):
    """Apply a preset by name."""
    path = os.path.join(ROOT, "presets", f"{name}.toml")
    if not os.path.isfile(path):
        print(f"ERROR: preset '{name}' not found.")
        sys.exit(1)
    apply_overlay(path)


def override_param(param_path, value):
    """Override a single dotted TOML parameter in bdm.toml.

    param_path: e.g. "skin.immune.cytokine_rate"
    value: the value to set (number, bool, or string)
    """
    bdm_toml = os.path.join(ROOT, "bdm.toml")
    with open(bdm_toml) as f:
        lines = f.readlines()

    # Parse dotted path: "skin.immune.cytokine_rate" -> section "[skin.immune]", key "cytokine_rate"
    parts = param_path.split(".")
    key = parts[-1]
    section = ".".join(parts[:-1])
    section_header = f"[{section}]"

    # Find the section, then the key within it
    in_section = False
    found = False
    for i, line in enumerate(lines):
        stripped = line.strip()
        if stripped.startswith("["):
            in_section = (stripped == section_header)
            if found:
                break  # past our section
        elif in_section:
            # Match key = value (ignoring comments)
            m = re.match(rf"^(\s*){re.escape(key)}\s*=\s*", line)
            if m:
                indent = m.group(1)
                if isinstance(value, bool):
                    val_str = "true" if value else "false"
                elif isinstance(value, (int, float)):
                    val_str = str(value)
                else:
                    val_str = f'"{value}"'
                lines[i] = f"{indent}{key} = {val_str}\n"
                found = True

    if not found:
        print(f"WARNING: param '{param_path}' not found in bdm.toml")

    with open(bdm_toml, "w") as f:
        f.writelines(lines)

    return found


# ---------------------------------------------------------------------------
# Build
# ---------------------------------------------------------------------------

def build_if_needed():
    """Build the binary if sources are newer."""
    binary = os.path.join(ROOT, "build", "skibidy")
    if os.path.isfile(binary):
        bin_mtime = os.path.getmtime(binary)
        needs_build = False
        for root, _, files in os.walk(os.path.join(ROOT, "src")):
            for f in files:
                if os.path.getmtime(os.path.join(root, f)) > bin_mtime:
                    needs_build = True
                    break
            if needs_build:
                break
        if os.path.getmtime(os.path.join(ROOT, "CMakeLists.txt")) > bin_mtime:
            needs_build = True
        if not needs_build:
            return
    print("Building...")
    rc = subprocess.run(
        "cd build && cmake .. -DCMAKE_BUILD_TYPE=Release > /dev/null 2>&1 "
        "&& make -j$(nproc) 2>&1 | tail -3",
        shell=True, cwd=ROOT)
    if rc.returncode != 0:
        print("Build failed.")
        sys.exit(1)


# ---------------------------------------------------------------------------
# Simulation execution
# ---------------------------------------------------------------------------

def run_simulation():
    """Run one simulation. Returns (success, elapsed_seconds).

    Success is determined by metrics.csv existing (not exit code),
    because ParaView viz export can crash headless without affecting
    simulation results.
    """
    output_dir = os.path.join(ROOT, "output")
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir, ignore_errors=True)

    # Disable viz export for batch runs (headless pvbatch crashes, wastes time)
    override_param("visualization.export", False)

    t0 = time.time()
    subprocess.run(
        [os.path.join(ROOT, "build", "skibidy")],
        cwd=ROOT, capture_output=True, text=True)
    elapsed = time.time() - t0

    return os.path.isfile(get_metrics_path()), elapsed


def get_metrics_path():
    """Return path to the metrics CSV from the last run."""
    return os.path.join(ROOT, "output", "skibidy", "metrics.csv")


# ---------------------------------------------------------------------------
# CSV loading and aggregation
# ---------------------------------------------------------------------------

def load_csv(path):
    """Load a CSV into a dict of lists (numeric where possible).

    Strips NUL bytes from corrupted files (partial writes from crashed sims).
    """
    with open(path, errors="replace") as f:
        clean = (row.replace("\x00", "") for row in f if not row.startswith("#"))
        reader = csv.DictReader(clean)
        rows = [r for r in reader if all(v is not None for v in r.values())]
    if not rows:
        return {}
    data = {}
    for key in rows[0]:
        try:
            data[key] = [float(r[key]) for r in rows]
        except (ValueError, TypeError):
            data[key] = [r[key] for r in rows]
    return data


def aggregate_csvs(csv_paths, min_length_fraction=0.5):
    """Compute mean and std across multiple CSV files.

    Runs shorter than min_length_fraction of the median length are
    excluded as likely early terminations.  For remaining runs, the
    consensus extends to the longest run; rows with fewer contributors
    average over whatever is available.

    Returns (mean_data, std_data) where each is a dict of lists.
    """
    all_data = [d for d in (load_csv(p) for p in csv_paths) if d]
    if not all_data:
        return {}, {}

    columns = [k for k in all_data[0] if isinstance(all_data[0][k][0], float)]
    if not columns:
        return {}, {}

    # Filter out abnormally short runs
    lengths = [len(d[columns[0]]) for d in all_data]
    lengths_sorted = sorted(lengths)
    median_len = lengths_sorted[len(lengths_sorted) // 2]
    min_len = max(1, int(median_len * min_length_fraction))

    filtered = []
    n_dropped = 0
    for d, length in zip(all_data, lengths):
        if length >= min_len:
            filtered.append(d)
        else:
            n_dropped += 1

    if n_dropped:
        print(f"  Dropped {n_dropped} short run(s) "
              f"(< {min_len} rows, median {median_len})")

    if not filtered:
        filtered = all_data  # fallback: keep everything

    n_rows = max(len(d[columns[0]]) for d in filtered)

    mean_data = {}
    std_data = {}
    for col in columns:
        means = []
        stds = []
        for i in range(n_rows):
            vals = [d[col][i] for d in filtered if i < len(d[col])]
            if not vals:
                break
            m = sum(vals) / len(vals)
            v = sum((x - m) ** 2 for x in vals) / len(vals)
            means.append(m)
            stds.append(math.sqrt(v))
        mean_data[col] = means
        std_data[col] = stds

    return mean_data, std_data


def write_csv(data, path, columns=None):
    """Write a dict-of-lists to CSV."""
    if columns is None:
        columns = list(data.keys())
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(columns)
        n_rows = len(data[columns[0]])
        for i in range(n_rows):
            writer.writerow([f"{data[col][i]:.6g}" for col in columns])


# ---------------------------------------------------------------------------
# Outcome extraction
# ---------------------------------------------------------------------------

def extract_outcome(data, column, measure="final"):
    """Extract a scalar outcome from a metrics timeseries.

    Measures:
      final     - last value
      peak      - maximum value
      auc       - area under curve (trapezoidal, using time_h)
      time_to_N - first time value exceeds N (e.g. time_to_90)
    """
    values = data.get(column, [])
    if not values:
        return float("nan")

    if measure == "final":
        return values[-1]
    elif measure == "peak":
        return max(values)
    elif measure == "auc":
        times = data.get("time_h", list(range(len(values))))
        total = 0
        for i in range(1, len(values)):
            dt = times[i] - times[i - 1]
            total += 0.5 * (values[i] + values[i - 1]) * dt
        return total
    elif measure.startswith("time_to_"):
        threshold = float(measure.split("_")[-1])
        times = data.get("time_h", list(range(len(values))))
        for t, v in zip(times, values):
            if v >= threshold:
                return t
        return float("nan")  # never reached
    else:
        raise ValueError(f"Unknown measure: {measure}")


# ---------------------------------------------------------------------------
# Validation integration
# ---------------------------------------------------------------------------

def run_validation(metrics_path):
    """Run validate_all.py on a metrics CSV. Returns returncode."""
    rc = subprocess.run(
        ["python3", os.path.join(ROOT, "literature", "validate_all.py"),
         metrics_path],
        cwd=ROOT)
    return rc.returncode
