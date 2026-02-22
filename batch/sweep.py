#!/usr/bin/env python3
"""Parameter sweep runner.

Reads a TOML sweep config, runs simulations across parameter values
with replicates, and produces summary statistics.

Usage:
    python3 batch/sweep.py batch/configs/cytokine_rate.toml
    python3 batch/sweep.py batch/configs/cytokine_rate.toml -n 3
    python3 batch/sweep.py batch/configs/collagen_vs_scar.toml --analyze
"""

import argparse
import itertools
import os
import subprocess
import sys

# Add project root to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir))

from datetime import datetime

from batch import lib

try:
    import tomllib
except ImportError:
    try:
        import tomli as tomllib
    except ImportError:
        tomllib = None


def parse_toml(path):
    """Parse TOML sweep config. Falls back to manual parsing if no tomllib."""
    if tomllib:
        with open(path, "rb") as f:
            return tomllib.load(f)
    # Minimal fallback parser for sweep configs
    return _parse_toml_fallback(path)


def _parse_toml_fallback(path):
    """Minimal TOML parser for sweep config format."""
    import ast
    import re

    config = {}
    current_section = config
    current_path = []
    array_sections = {}  # track [[section]] arrays

    with open(path) as f:
        for line in f:
            line = line.split("#")[0].strip()
            if not line:
                continue

            # Array of tables: [[sweep.params]]
            m = re.match(r"^\[\[(.+)\]\]$", line)
            if m:
                key = m.group(1)
                parts = key.split(".")
                node = config
                for p in parts[:-1]:
                    node = node.setdefault(p, {})
                arr = node.setdefault(parts[-1], [])
                new_table = {}
                arr.append(new_table)
                current_section = new_table
                continue

            # Table: [section]
            m = re.match(r"^\[(.+)\]$", line)
            if m:
                current_path = m.group(1).split(".")
                node = config
                for p in current_path:
                    node = node.setdefault(p, {})
                current_section = node
                continue

            # Key = value
            m = re.match(r"^(\w+)\s*=\s*(.+)$", line)
            if m:
                key, val = m.group(1), m.group(2).strip()
                # Parse value
                if val.startswith("["):
                    current_section[key] = ast.literal_eval(val)
                elif val.startswith('"'):
                    current_section[key] = val.strip('"')
                elif val == "true":
                    current_section[key] = True
                elif val == "false":
                    current_section[key] = False
                elif "." in val:
                    current_section[key] = float(val)
                else:
                    try:
                        current_section[key] = int(val)
                    except ValueError:
                        current_section[key] = val

    return config


def load_sweep_config(path):
    """Load and validate a sweep config."""
    cfg = parse_toml(path)
    sweep = cfg.get("sweep", {})
    outcomes = cfg.get("outcomes", {})

    name = sweep.get("name", os.path.splitext(os.path.basename(path))[0])
    preset = sweep.get("preset", None)
    skin = sweep.get("skin", None)
    runs = sweep.get("runs_per_value", 5)

    # Single param sweep
    if "param" in sweep:
        params = [{"param": sweep["param"], "values": sweep["values"]}]
    # Multi param sweep
    elif "params" in sweep:
        params = sweep["params"]
    else:
        print("ERROR: sweep config must have 'param' or [[sweep.params]]")
        sys.exit(1)

    primary = outcomes.get("primary", "wound_closure_pct")
    secondary = outcomes.get("secondary", [])
    measure = outcomes.get("measure", "final")

    return dict(
        name=name, preset=preset, skin=skin, runs_per_value=runs,
        params=params, primary=primary, secondary=secondary,
        measure=measure, source_path=path,
    )


def generate_sweep_points(params):
    """Generate all combinations of parameter values.

    Returns list of dicts: [{param_path: value, ...}, ...]
    """
    param_names = [p["param"] for p in params]
    value_lists = [p["values"] for p in params]
    points = []
    for combo in itertools.product(*value_lists):
        points.append(dict(zip(param_names, combo)))
    return points


def run_sweep(cfg, override_runs=None):
    """Execute a full parameter sweep."""
    runs_per = override_runs or cfg["runs_per_value"]
    points = generate_sweep_points(cfg["params"])
    all_outcomes = cfg["secondary"] + [cfg["primary"]]

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    result_dir = os.path.join(
        os.path.dirname(__file__), "results",
        f"{cfg['name']}_{timestamp}")
    raw_dir = os.path.join(result_dir, "raw")
    os.makedirs(raw_dir, exist_ok=True)

    # Copy config for reproducibility
    import shutil
    shutil.copy2(cfg["source_path"], os.path.join(result_dir, "config.toml"))

    print(f"=== Sweep: {cfg['name']} ===")
    print(f"  Params: {[p['param'] for p in cfg['params']]}")
    print(f"  Points: {len(points)}")
    print(f"  Runs per point: {runs_per}")
    print(f"  Total simulations: {len(points) * runs_per}")
    print(f"  Output: {result_dir}/")
    print()

    # Build once
    lib.build_if_needed()

    # Summary data
    summary_rows = []
    total_sims = len(points) * runs_per
    sim_count = 0

    for pt_idx, point in enumerate(points):
        point_label = ", ".join(f"{k.split('.')[-1]}={v}" for k, v in point.items())
        print(f"--- Point {pt_idx + 1}/{len(points)}: {point_label} ---")

        csv_paths = []
        for run_idx in range(runs_per):
            sim_count += 1
            label = f"[{sim_count}/{total_sims}]"

            # Fresh config each run
            lib.merge_config()
            if cfg["skin"]:
                lib.apply_profile(cfg["skin"])
            if cfg["preset"]:
                lib.apply_preset(cfg["preset"])

            # Apply sweep parameter overrides
            for param_path, value in point.items():
                lib.override_param(param_path, value)

            # Run
            print(f"  {label} run {run_idx + 1}/{runs_per}...", end=" ", flush=True)
            ok, elapsed = lib.run_simulation()

            if not ok:
                print(f"FAILED ({elapsed:.0f}s)")
                continue

            # Save CSV
            src = lib.get_metrics_path()
            if not os.path.isfile(src):
                print(f"FAILED (no metrics)")
                continue

            # Encode point values in filename
            val_str = "_".join(f"{v}" for v in point.values())
            dst = os.path.join(raw_dir, f"pt{pt_idx:03d}_val{val_str}_run{run_idx:03d}.csv")
            shutil.copy2(src, dst)
            csv_paths.append(dst)
            print(f"OK ({elapsed:.0f}s)")

        # Aggregate this point
        if csv_paths:
            mean_data, std_data = lib.aggregate_csvs(csv_paths)
            row = dict(point)  # param values
            row["n_runs"] = len(csv_paths)
            for col in all_outcomes:
                val = lib.extract_outcome(mean_data, col, cfg["measure"])
                std_val = lib.extract_outcome(std_data, col, cfg["measure"])
                row[col] = val
                row[f"{col}_std"] = std_val
            summary_rows.append(row)

    # Write summary CSV
    if summary_rows:
        summary_path = os.path.join(result_dir, "summary.csv")
        columns = list(summary_rows[0].keys())
        with open(summary_path, "w", newline="") as f:
            import csv
            writer = csv.DictWriter(f, fieldnames=columns)
            writer.writeheader()
            for row in summary_rows:
                writer.writerow({k: f"{v:.6g}" if isinstance(v, float) else v
                                 for k, v in row.items()})
        print(f"\nSummary: {summary_path}")

    print(f"\nSweep complete: {result_dir}/")
    return result_dir


def main():
    parser = argparse.ArgumentParser(description="Parameter sweep runner")
    parser.add_argument("config", help="Path to sweep config TOML")
    parser.add_argument("-n", "--runs", type=int, default=None,
                        help="Override runs_per_value")
    parser.add_argument("--analyze", action="store_true",
                        help="Run analyze.py after sweep completes")
    args = parser.parse_args()

    cfg = load_sweep_config(args.config)
    result_dir = run_sweep(cfg, override_runs=args.runs)

    if args.analyze:
        subprocess.run([
            "python3", os.path.join(os.path.dirname(__file__), "analyze.py"),
            result_dir])


if __name__ == "__main__":
    main()
