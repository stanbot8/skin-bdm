#!/usr/bin/env python3
"""Experiment runner: execute named experiments from TOML configs.

Each experiment TOML defines a set of named configurations with parameter
overrides, run counts, and comparison mode.  The runner handles config
merging, simulation execution, consensus aggregation, and outcome
comparison automatically.

Usage:
    python3 scripts/study/experiment_runner.py studies/diabetic-wound/experiments/wound_size.toml
    python3 scripts/study/experiment_runner.py studies/*/experiments/*.toml
    python3 scripts/study/experiment_runner.py studies/diabetic-wound/experiments/biofilm_infection.toml --runs=3
"""

import argparse
import csv
import math
import os
import shutil
import sys
import time

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, REPO)

from batch.lib import (
    merge_config,
    apply_profile,
    apply_study,
    apply_overlay,
    override_param,
    run_simulation,
    get_metrics_path,
    load_csv,
    aggregate_csvs,
    extract_outcome,
    write_csv,
)

sys.path.insert(0, os.path.join(REPO, "literature"))
from cure_criteria import evaluate_cure, print_cure_assessment

from batch.lib import get_tomllib
tomllib = get_tomllib()
# ---------------------------------------------------------------------------
# Experiment loading
# ---------------------------------------------------------------------------

def load_experiment(path):
    """Parse and validate an experiment TOML file."""
    with open(path, "rb") as f:
        data = tomllib.load(f)

    experiment = data.get("experiment")
    if not experiment:
        raise ValueError(f"{path}: missing [experiment] table")

    required = ["name", "configs"]
    for key in required:
        if key not in experiment:
            raise ValueError(f"{path}: missing experiment.{key}")

    # Defaults
    experiment.setdefault("profile", "diabetic")
    experiment.setdefault("study", "diabetic-wound")
    experiment.setdefault("runs_per_config", 5)
    experiment.setdefault("description", "")

    # Validate configs
    for i, cfg in enumerate(experiment["configs"]):
        if "label" not in cfg:
            raise ValueError(f"{path}: configs[{i}] missing 'label'")
        cfg.setdefault("overrides", {})
        cfg.setdefault("profile", None)
        cfg.setdefault("study", None)
        cfg.setdefault("treatments", [])
        cfg.setdefault("extra_overlays", [])

    return experiment


# ---------------------------------------------------------------------------
# Config preparation
# ---------------------------------------------------------------------------

def prepare_experiment_config(cfg, experiment):
    """Build bdm.toml for one experiment config entry.

    Merges core config, applies profile and study config from the experiment
    (or per-config overrides), applies treatments, then applies
    parameter overrides.
    """
    bdm_path = os.path.join(REPO, "bdm.toml")
    if os.path.exists(bdm_path):
        os.remove(bdm_path)

    # Merge base config
    merge_config()

    # Profile: per-config overrides experiment-level
    profile = cfg.get("profile") or experiment.get("profile")
    if profile:
        apply_profile(profile)

    # Study config: per-config overrides experiment-level
    study = cfg.get("study") or experiment.get("study")
    if study:
        apply_study(study)

    # Treatments (search study-scoped treatments first, then all studies)
    for tname in cfg.get("treatments", []):
        tpath = None
        if study:
            candidate = os.path.join(REPO, "studies", study, "treatments", f"{tname}.toml")
            if os.path.isfile(candidate):
                tpath = candidate
        if not tpath:
            import glob as _glob
            matches = sorted(_glob.glob(os.path.join(
                REPO, "studies", "*", "treatments", f"{tname}.toml")))
            if matches:
                tpath = matches[0]
        if tpath:
            apply_overlay(tpath)
        else:
            print(f"  WARNING: treatment '{tname}' not found", flush=True)

    # Extra overlays (arbitrary TOML files)
    for overlay in cfg.get("extra_overlays", []):
        opath = os.path.join(REPO, overlay)
        if os.path.isfile(opath):
            apply_overlay(opath)
        else:
            print(f"  WARNING: overlay '{overlay}' not found", flush=True)

    # Parameter overrides
    for param_path, value in cfg.get("overrides", {}).items():
        override_param(param_path, value)

    # Strip visualization for headless batch
    _strip_viz(bdm_path)


def _strip_viz(bdm_path):
    """Remove visualization sections and disable autoopen."""
    with open(bdm_path) as f:
        lines = f.readlines()
    out = []
    skip = False
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("[visualization") or stripped.startswith("[[visualize"):
            skip = True
            continue
        if skip and stripped.startswith("[") and not stripped.startswith("[visualization"):
            skip = False
        if skip:
            continue
        if stripped.startswith("metrics_autoopen"):
            out.append(line.replace("true", "false"))
            continue
        out.append(line)
    with open(bdm_path, "w") as f:
        f.writelines(out)


# ---------------------------------------------------------------------------
# Simulation execution with consensus
# ---------------------------------------------------------------------------

def run_consensus(label, cfg, experiment, n_runs, output_dir):
    """Run n_runs simulations for one config and aggregate.

    Returns (mean_data, std_data, csv_paths, outcomes).
    """
    csv_paths = []
    raw_dir = os.path.join(output_dir, "raw")
    os.makedirs(raw_dir, exist_ok=True)

    for i in range(n_runs):
        prepare_experiment_config(cfg, experiment)
        safe_label = label.replace(" ", "_").replace("+", "_").replace("(", "").replace(")", "")
        run_dir = os.path.join(raw_dir, f"{safe_label}_run{i:03d}")
        success, elapsed = run_simulation(output_path=run_dir)

        if not success:
            print(f"    Run {i+1}/{n_runs} failed, skipping", flush=True)
            continue

        csv_path = get_metrics_path()
        if not os.path.isfile(csv_path):
            print(f"    Run {i+1}/{n_runs} failed (no metrics)", flush=True)
            continue
        csv_paths.append(csv_path)
        print(f"    Run {i+1}/{n_runs} done ({elapsed:.0f}s)", flush=True)

    if not csv_paths:
        return {}, {}, [], {}

    # Aggregate
    mean_data, std_data = aggregate_csvs(csv_paths)

    # Extract scalar outcomes
    outcomes = {}
    if mean_data:
        outcomes["wound_closure_pct"] = extract_outcome(mean_data, "wound_closure_pct", "final")
        outcomes["peak_inflammation"] = extract_outcome(mean_data, "mean_infl_wound", "peak")
        outcomes["scar_magnitude"] = extract_outcome(mean_data, "scar_magnitude", "final")
        outcomes["peak_neutrophils"] = extract_outcome(mean_data, "n_neutrophils", "peak")
        outcomes["peak_macrophages"] = extract_outcome(mean_data, "n_macrophages", "peak")
        outcomes["peak_collagen"] = extract_outcome(mean_data, "mean_collagen_wound", "peak")
        outcomes["time_to_50pct_h"] = extract_outcome(mean_data, "wound_closure_pct", "time_to_50")
        outcomes["time_to_90pct_h"] = extract_outcome(mean_data, "wound_closure_pct", "time_to_90")
        # Convert hours to days
        for key in ["time_to_50pct_h", "time_to_90pct_h"]:
            v = outcomes.get(key)
            day_key = key.replace("_h", "_days")
            outcomes[day_key] = round(v / 24, 1) if v and not math.isnan(v) else None

    return mean_data, std_data, csv_paths, outcomes


def _detect_cure_condition(experiment):
    """Map experiment profile/study to a cure criteria condition."""
    profile = experiment.get("profile", "")
    study = experiment.get("study", "")
    if profile == "rheumatoid" or "rheumatoid" in study:
        return "rheumatoid"
    if profile == "pressure" or "pressure" in study:
        return "pressure"
    if profile == "burn" or "burn" in study:
        return "burn"
    if profile == "diabetic" or "diabetic" in study:
        return "diabetic"
    return None


# ---------------------------------------------------------------------------
# Comparison and output
# ---------------------------------------------------------------------------

def write_comparison(results, output_dir):
    """Write comparison CSV with all configs side by side."""
    if not results:
        return

    path = os.path.join(output_dir, "comparison.csv")
    outcome_keys = set()
    for r in results:
        outcome_keys.update(r.get("outcomes", {}).keys())
    outcome_keys = sorted(outcome_keys)

    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["config", "n_runs"] + outcome_keys)
        for r in results:
            row = [r["label"], r.get("n_runs", 0)]
            for key in outcome_keys:
                val = r.get("outcomes", {}).get(key, "")
                if isinstance(val, float):
                    row.append(f"{val:.6g}")
                else:
                    row.append(val if val is not None else "")
            writer.writerow(row)

    print(f"  Comparison: {path}", flush=True)
    return path


def write_summary(experiment, results, output_dir, elapsed):
    """Write human-readable summary."""
    path = os.path.join(output_dir, "summary.txt")
    with open(path, "w") as f:
        f.write(f"Experiment: {experiment['name']}\n")
        if experiment.get("description"):
            f.write(f"{experiment['description']}\n")
        f.write(f"\nProfile: {experiment.get('profile', 'default')}\n")
        f.write(f"Study: {experiment.get('study', 'default')}\n")
        f.write(f"Runs per config: {experiment.get('runs_per_config', 5)}\n")
        f.write(f"Total time: {elapsed/60:.1f} minutes\n")
        f.write(f"\n{'='*70}\n")
        f.write(f"{'Config':<30} {'Closure':>8} {'T50 (d)':>8} {'Peak Infl':>10} {'Scar':>8}\n")
        f.write(f"{'-'*70}\n")

        for r in results:
            o = r.get("outcomes", {})
            label = r["label"][:30]
            closure = o.get("wound_closure_pct", 0)
            t50 = o.get("time_to_50pct_days")
            t50_str = f"{t50:.1f}" if t50 is not None else "N/A"
            peak = o.get("peak_inflammation", 0)
            scar = o.get("scar_magnitude", 0)
            f.write(f"{label:<30} {closure:>7.1f}% {t50_str:>8} {peak:>10.4f} {scar:>8.3f}\n")

        # Comparison vs first config (baseline)
        if len(results) > 1 and results[0].get("outcomes"):
            base = results[0]["outcomes"]
            base_closure = base.get("wound_closure_pct", 0)
            f.write(f"\n{'='*70}\n")
            f.write("Comparison vs first config:\n\n")
            for r in results[1:]:
                o = r.get("outcomes", {})
                closure = o.get("wound_closure_pct", 0)
                delta = closure - base_closure
                f.write(f"  {r['label']:<30} closure: {delta:+.1f}%\n")

    print(f"  Summary: {path}", flush=True)


def print_results_table(results):
    """Print a results table to stdout."""
    if not results:
        return

    print(f"\n  {'Config':<30} {'Closure':>8} {'T50 (d)':>8} {'Peak Infl':>10} {'Scar':>8}")
    print(f"  {'-'*66}")

    for r in results:
        o = r.get("outcomes", {})
        label = r["label"][:30]
        closure = o.get("wound_closure_pct", 0)
        t50 = o.get("time_to_50pct_days")
        t50_str = f"{t50:.1f}" if t50 is not None else "N/A"
        peak = o.get("peak_inflammation", 0)
        scar = o.get("scar_magnitude", 0)
        print(f"  {label:<30} {closure:>7.1f}% {t50_str:>8} {peak:>10.4f} {scar:>8.3f}")


# ---------------------------------------------------------------------------
# Main orchestration
# ---------------------------------------------------------------------------

def run_experiment_file(experiment_path, runs_override=None):
    """Load and run a single experiment file."""
    experiment = load_experiment(experiment_path)
    name = experiment["name"]
    n_runs = runs_override or experiment["runs_per_config"]

    # Output into the study's results directory
    study_name = experiment.get("study", "diabetic-wound")
    safe_name = name.lower().replace(" ", "_").replace("(", "").replace(")", "")
    output_dir = os.path.join(REPO, "studies", study_name, "results", "experiments", safe_name)
    os.makedirs(output_dir, exist_ok=True)

    n_configs = len(experiment["configs"])
    total_runs = n_configs * n_runs

    print(f"\n{'='*60}")
    print(f"  Experiment: {name}")
    if experiment.get("description"):
        print(f"  {experiment['description']}")
    print(f"  {n_configs} configs x {n_runs} runs = {total_runs} total")
    print(f"{'='*60}\n")

    results = []
    t0 = time.time()

    for idx, cfg in enumerate(experiment["configs"]):
        label = cfg["label"]
        print(f"[{idx+1}/{n_configs}] {label}", flush=True)

        mean_data, std_data, csv_paths, outcomes = run_consensus(
            label, cfg, experiment, n_runs, output_dir)

        # Save consensus CSV
        if mean_data:
            safe_label = label.replace(" ", "_").replace("+", "_").replace("(", "").replace(")", "")
            consensus_path = os.path.join(output_dir, f"consensus_{safe_label}.csv")
            cols = sorted(mean_data.keys())
            write_csv(mean_data, consensus_path, cols)

        result = {
            "label": label,
            "n_runs": len(csv_paths),
            "outcomes": outcomes,
            "mean_data": mean_data,
        }
        results.append(result)

        closure = outcomes.get("wound_closure_pct", 0)
        print(f"  >> Closure: {closure:.1f}% ({len(csv_paths)} runs)\n", flush=True)

    elapsed = time.time() - t0

    # Write outputs
    write_comparison(results, output_dir)
    write_summary(experiment, results, output_dir, elapsed)
    print_results_table(results)

    # Cure assessment (if applicable to this condition)
    cure_cond = _detect_cure_condition(experiment)
    if cure_cond:
        print(f"\n  Cure Assessment ({cure_cond}):")
        for r in results:
            if r.get("mean_data"):
                cure_r = evaluate_cure(r["mean_data"], cure_cond)
                r["cure"] = cure_r
                print(f"    {r['label']}: {cure_r['verdict']} ({cure_r['score']:.0%})")
                for c in cure_r["criteria"]:
                    mark = "+" if c["met"] else "-"
                    print(f"      [{mark}] {c['name']}: {c['detail']}")

    mins = int(elapsed) // 60
    secs = int(elapsed) % 60
    print(f"\n  Completed in {mins}m {secs}s")
    print(f"  Output: {output_dir}/\n")

    return results


def main():
    parser = argparse.ArgumentParser(
        description="Run experiment experiments from TOML configs")
    parser.add_argument("experiments", nargs="+",
                        help="Path(s) to experiment TOML file(s)")
    parser.add_argument("--runs", type=int, default=None,
                        help="Override runs per config")
    args = parser.parse_args()

    for path in args.experiments:
        if not os.path.isfile(path):
            print(f"ERROR: experiment file not found: {path}")
            sys.exit(1)

    for path in args.experiments:
        run_experiment_file(path, runs_override=args.runs)


if __name__ == "__main__":
    main()
