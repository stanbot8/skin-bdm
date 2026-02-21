#!/usr/bin/env python3
"""Treatment study: run diabetic wound with each treatment and compare outcomes.

Usage:
    python3 scripts/study/treatment_study.py [--treatments=all|name1,name2,...] [--steps=N]
    python3 scripts/study/treatment_study.py --combos=pairs   # singles + all pairwise
    python3 scripts/study/treatment_study.py --combos=all     # all 2^N-1 subsets
"""

import argparse
import csv
import itertools
import os
import subprocess
import sys
import shutil
import time

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BINARY = os.path.join(REPO, "build", "skibidy")
OUTPUT_BASE = os.path.join(REPO, "output", "treatment_study")
TREATMENTS_DIR = os.path.join(REPO, "treatments")


def available_treatments():
    """List all treatment TOML files."""
    tdir = TREATMENTS_DIR
    return sorted(
        os.path.splitext(f)[0]
        for f in os.listdir(tdir)
        if f.endswith(".toml")
    )


def run_script(args, label):
    """Run a helper script, printing errors on failure."""
    result = subprocess.run(args, cwd=REPO, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"FAILED", flush=True)
        if result.stderr:
            for line in result.stderr.strip().splitlines():
                print(f"      {line}", flush=True)
        if result.stdout:
            for line in result.stdout.strip().splitlines():
                print(f"      {line}", flush=True)
        raise RuntimeError(f"{label} failed (exit {result.returncode})")
    print("ok", flush=True)


def prepare_config(treatment_names=None):
    """Generate bdm.toml for one or more treatments (or baseline diabetic).

    Args:
        treatment_names: None for baseline, a string for single treatment,
                        or a list of strings for combination therapy.
    """
    if isinstance(treatment_names, str):
        treatment_names = [treatment_names]

    # Step 0: clear stale config
    bdm_path = os.path.join(REPO, "bdm.toml")
    if os.path.exists(bdm_path):
        os.remove(bdm_path)

    print("    Merging config...", end=" ", flush=True)
    run_script(
        [sys.executable, os.path.join(REPO, "scripts", "config", "merge_config.py")],
        "merge_config",
    )

    print("    Applying diabetic profile...", end=" ", flush=True)
    run_script(
        [sys.executable, os.path.join(REPO, "scripts", "config", "apply_preset.py"),
         os.path.join(REPO, "profiles", "diabetic.toml"),
         os.path.join(REPO, "bdm.toml")],
        "apply diabetic profile",
    )

    print("    Applying wound preset...", end=" ", flush=True)
    run_script(
        [sys.executable, os.path.join(REPO, "scripts", "config", "apply_preset.py"),
         os.path.join(REPO, "presets", "diabetic_wound.toml"),
         os.path.join(REPO, "bdm.toml")],
        "apply wound preset",
    )

    if treatment_names:
        for tname in treatment_names:
            print(f"    Applying treatment: {tname}...", end=" ", flush=True)
            treatment_file = os.path.join(TREATMENTS_DIR, f"{tname}.toml")
            run_script(
                [sys.executable, os.path.join(REPO, "scripts", "config", "apply_preset.py"),
                 treatment_file,
                 os.path.join(REPO, "bdm.toml")],
                f"apply treatment {tname}",
            )

    # Patch config for headless batch runs:
    # - strip all visualization sections (prevents BioDynaMo GL init)
    # - disable metrics_autoopen (prevents xdg-open after each sim)
    bdm_toml = os.path.join(REPO, "bdm.toml")
    with open(bdm_toml) as f:
        lines = f.readlines()
    out = []
    skip = False
    for line in lines:
        stripped = line.strip()
        # Skip any section starting with [visualization or [[visualize
        if stripped.startswith("[visualization") or stripped.startswith("[[visualize"):
            skip = True
            continue
        # Stop skipping when we hit a non-visualization section
        if skip and stripped.startswith("[") and not stripped.startswith("[visualization"):
            skip = False
        if skip:
            continue
        # Disable metrics auto-open
        if stripped.startswith("metrics_autoopen"):
            out.append(line.replace("true", "false"))
            continue
        out.append(line)
    with open(bdm_toml, "w") as f:
        f.writelines(out)


def run_simulation(label):
    """Run the simulation and return metrics CSV path."""
    # Clean sim output only (preserve treatment_study/ results)
    sim_dir = os.path.join(REPO, "output", "skibidy")
    if os.path.exists(sim_dir):
        shutil.rmtree(sim_dir)

    print(f"    Simulating...", flush=True)
    print(f"    {'~' * 50}", flush=True)
    t0 = time.time()
    env = os.environ.copy()
    env.pop("DISPLAY", None)  # prevent BioDynaMo GL init in batch mode
    proc = subprocess.Popen(
        [BINARY], cwd=REPO, env=env,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        text=True, bufsize=1,
    )
    for line in proc.stdout:
        line = line.rstrip()
        if line:
            print(f"    {line}", flush=True)
    rc = proc.wait()
    elapsed = time.time() - t0
    mins = int(elapsed) // 60
    secs = int(elapsed) % 60
    print(f"    {'~' * 50}", flush=True)
    if rc != 0:
        print(f"    WARNING: simulation exited with code {rc}", flush=True)
    print(f"    Done ({mins}m {secs}s)", flush=True)

    metrics_path = os.path.join(REPO, "output", "skibidy", "metrics.csv")
    if not os.path.exists(metrics_path):
        print(f"  WARNING: metrics.csv not found for {label}")
        return None
    return metrics_path


def extract_outcomes(metrics_path):
    """Extract key outcome metrics from a simulation run."""
    if not metrics_path or not os.path.exists(metrics_path):
        return {}

    with open(metrics_path) as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    if not rows:
        return {}

    outcomes = {}

    # Final values
    last = rows[-1]
    outcomes["wound_closure_pct"] = float(last.get("wound_closure_pct", 0))
    outcomes["mean_infl_wound"] = float(last.get("mean_infl_wound", 0))
    outcomes["scar_magnitude"] = float(last.get("scar_magnitude", 0))
    outcomes["n_agents"] = int(float(last.get("n_total", 0)))

    # Peak inflammation
    peak_infl = 0
    peak_infl_step = 0
    for row in rows:
        infl = float(row.get("mean_infl_wound", 0))
        if infl > peak_infl:
            peak_infl = infl
            peak_infl_step = int(float(row.get("step", 0)))
    outcomes["peak_inflammation"] = peak_infl
    outcomes["peak_inflammation_step"] = peak_infl_step
    outcomes["peak_inflammation_day"] = round(peak_infl_step * 0.1 / 24, 1)

    # Time to 50% closure
    for row in rows:
        closure = float(row.get("wound_closure_pct", 0))
        if closure >= 50:
            step = int(float(row.get("step", 0)))
            outcomes["time_to_50pct_days"] = round(step * 0.1 / 24, 1)
            break
    else:
        outcomes["time_to_50pct_days"] = None

    # Time to 90% closure
    for row in rows:
        closure = float(row.get("wound_closure_pct", 0))
        if closure >= 90:
            step = int(float(row.get("step", 0)))
            outcomes["time_to_90pct_days"] = round(step * 0.1 / 24, 1)
            break
    else:
        outcomes["time_to_90pct_days"] = None

    # Peak immune cells
    peak_neut = 0
    peak_mac = 0
    for row in rows:
        n = int(float(row.get("n_neutrophils", 0)))
        m = int(float(row.get("n_macrophages", 0)))
        if n > peak_neut:
            peak_neut = n
        if m > peak_mac:
            peak_mac = m
    outcomes["peak_neutrophils"] = peak_neut
    outcomes["peak_macrophages"] = peak_mac

    # Collagen metrics
    peak_col = 0
    for row in rows:
        col = float(row.get("mean_collagen_wound", 0))
        if col > peak_col:
            peak_col = col
    outcomes["peak_collagen"] = peak_col

    # Fibroblast metrics
    peak_myofib = 0
    for row in rows:
        mf = int(float(row.get("n_myofibroblasts", 0)))
        if mf > peak_myofib:
            peak_myofib = mf
    outcomes["peak_myofibroblasts"] = peak_myofib

    return outcomes


def print_comparison(results):
    """Print a comparison table of all results."""
    if not results:
        print("No results to compare.")
        return

    print("\n" + "=" * 100)
    print("TREATMENT STUDY RESULTS")
    print("=" * 100)

    # Header
    metrics = [
        ("wound_closure_pct", "Closure %", ".1f"),
        ("time_to_50pct_days", "50% (days)", ""),
        ("time_to_90pct_days", "90% (days)", ""),
        ("peak_inflammation", "Peak Infl", ".4f"),
        ("peak_inflammation_day", "Infl Peak (d)", ".1f"),
        ("mean_infl_wound", "Final Infl", ".4f"),
        ("peak_neutrophils", "Peak Neut", "d"),
        ("peak_macrophages", "Peak Mac", "d"),
        ("peak_myofibroblasts", "Peak Myofib", "d"),
        ("peak_collagen", "Peak Collagen", ".5f"),
        ("scar_magnitude", "Scar", ".3f"),
    ]

    # Find widths
    name_width = max(len(name) for name, _ in results)
    name_width = max(name_width, 15)

    print(f"\n{'Treatment':<{name_width}}", end="")
    for _, label, _ in metrics:
        print(f"  {label:>13}", end="")
    print()
    print("-" * (name_width + len(metrics) * 15))

    baseline_outcomes = results[0][1] if results else {}

    for name, outcomes in results:
        if not outcomes:
            print(f"{name:<{name_width}}  (failed)")
            continue

        print(f"{name:<{name_width}}", end="")
        for key, _, fmt in metrics:
            val = outcomes.get(key)
            if val is None:
                print(f"  {'N/A':>13}", end="")
            elif fmt == "d":
                print(f"  {int(val):>13d}", end="")
            elif fmt:
                print(f"  {val:>13{fmt}}", end="")
            else:
                if val is None:
                    print(f"  {'N/A':>13}", end="")
                else:
                    print(f"  {val:>13.1f}", end="")
        print()

    # Improvement vs baseline
    if len(results) > 1 and baseline_outcomes:
        print()
        print(f"{'IMPROVEMENT vs baseline':<{name_width}}")
        print("-" * (name_width + len(metrics) * 15))
        baseline_closure = baseline_outcomes.get("wound_closure_pct", 0)

        for name, outcomes in results[1:]:
            if not outcomes:
                continue
            closure = outcomes.get("wound_closure_pct", 0)
            delta_closure = closure - baseline_closure

            t50_base = baseline_outcomes.get("time_to_50pct_days")
            t50_treat = outcomes.get("time_to_50pct_days")
            if t50_base and t50_treat:
                delta_t50 = f"{t50_base - t50_treat:+.1f}d"
            else:
                delta_t50 = "N/A"

            t90_base = baseline_outcomes.get("time_to_90pct_days")
            t90_treat = outcomes.get("time_to_90pct_days")
            if t90_base and t90_treat:
                delta_t90 = f"{t90_base - t90_treat:+.1f}d"
            else:
                delta_t90 = "N/A"

            peak_infl_base = baseline_outcomes.get("peak_inflammation", 0)
            peak_infl_treat = outcomes.get("peak_inflammation", 0)
            if peak_infl_base > 0:
                delta_infl = f"{(peak_infl_treat - peak_infl_base) / peak_infl_base * 100:+.0f}%"
            else:
                delta_infl = "N/A"

            print(f"{name:<{name_width}}  closure: {delta_closure:+.1f}%  "
                  f"50%: {delta_t50}  90%: {delta_t90}  "
                  f"peak_infl: {delta_infl}")

    print()


def save_results(results, output_dir):
    """Save results to CSV."""
    os.makedirs(output_dir, exist_ok=True)
    csv_path = os.path.join(output_dir, "treatment_comparison.csv")
    if not results:
        return

    # Collect all keys
    all_keys = set()
    for _, outcomes in results:
        if outcomes:
            all_keys.update(outcomes.keys())
    all_keys = sorted(all_keys)

    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["treatment"] + all_keys)
        for name, outcomes in results:
            row = [name]
            for key in all_keys:
                row.append(outcomes.get(key, "") if outcomes else "")
            writer.writerow(row)

    print(f"Results saved to {csv_path}")


def combo_label(names):
    """Human-readable label for a treatment combination."""
    return "+".join(sorted(names))


def generate_combos(treatments, mode):
    """Generate treatment combinations based on mode.

    Args:
        treatments: list of single treatment names (excludes 'combination')
        mode: 'pairs' for all 2-combos, 'triples' for 2+3-combos,
              'all' for all 2^N-1 subsets of size >= 2

    Returns:
        list of tuples, each a combination of treatment names
    """
    combos = []
    if mode == "pairs":
        combos = list(itertools.combinations(treatments, 2))
    elif mode == "triples":
        combos = list(itertools.combinations(treatments, 2))
        combos += list(itertools.combinations(treatments, 3))
    elif mode == "all":
        for r in range(2, len(treatments) + 1):
            combos += list(itertools.combinations(treatments, r))
    return combos


def main():
    parser = argparse.ArgumentParser(description="Treatment study for diabetic wounds")
    parser.add_argument("--treatments", default="all",
                        help="Comma-separated treatment names, or 'all'")
    parser.add_argument("--combos", default=None,
                        choices=["pairs", "triples", "all"],
                        help="Generate combinations: pairs (C(n,2)), "
                             "triples (2+3-combos), all (2^N-1 subsets)")
    parser.add_argument("--steps", type=int, default=None,
                        help="Override num_steps (default: use preset)")
    args = parser.parse_args()

    all_treatments = available_treatments()

    if args.treatments == "all":
        treatments = all_treatments
    else:
        treatments = [t.strip() for t in args.treatments.split(",")]
        for t in treatments:
            if t not in all_treatments:
                print(f"ERROR: unknown treatment '{t}'. Available: {all_treatments}")
                sys.exit(1)

    if not os.path.exists(BINARY):
        print("ERROR: binary not found. Run build first.")
        sys.exit(1)

    # For auto-combos, exclude the hand-tuned 'combination' preset
    combo_treatments = [t for t in treatments if t != "combination"]
    combos = generate_combos(combo_treatments, args.combos) if args.combos else []

    total = 1 + len(treatments) + len(combos)
    print(f"Treatment study: {total} simulations", flush=True)
    print(f"  1 baseline + {len(treatments)} singles", end="", flush=True)
    if combos:
        print(f" + {len(combos)} combinations ({args.combos})", end="", flush=True)
    print(flush=True)
    print(flush=True)

    results = []
    run_idx = 1

    # Baseline: diabetic wound, no treatment
    print(f"[{run_idx}/{total}] BASELINE (diabetic, no treatment)", flush=True)
    prepare_config(treatment_names=None)
    metrics = run_simulation("baseline")
    if metrics:
        os.makedirs(OUTPUT_BASE, exist_ok=True)
        shutil.copy2(metrics, os.path.join(OUTPUT_BASE, "metrics_baseline.csv"))
    outcomes = extract_outcomes(metrics)
    results.append(("baseline", outcomes))
    closure = outcomes.get("wound_closure_pct", 0)
    print(f"    >> Closure: {closure:.1f}%", flush=True)
    print(flush=True)
    run_idx += 1

    # Each single treatment
    for treatment in treatments:
        print(f"[{run_idx}/{total}] {treatment.upper()}", flush=True)
        try:
            prepare_config(treatment_names=treatment)
            metrics = run_simulation(treatment)
            if metrics:
                os.makedirs(OUTPUT_BASE, exist_ok=True)
                shutil.copy2(metrics, os.path.join(OUTPUT_BASE, f"metrics_{treatment}.csv"))
            outcomes = extract_outcomes(metrics)
            results.append((treatment, outcomes))
            closure = outcomes.get("wound_closure_pct", 0)
            base_closure = results[0][1].get("wound_closure_pct", 0) if results[0][1] else 0
            delta = closure - base_closure
            print(f"    >> Closure: {closure:.1f}% ({delta:+.1f}% vs baseline)", flush=True)
        except Exception as e:
            print(f"    ERROR: {e}", flush=True)
            results.append((treatment, None))
        print(flush=True)
        run_idx += 1

    # Combination treatments
    for combo in combos:
        label = combo_label(combo)
        print(f"[{run_idx}/{total}] COMBO: {label}", flush=True)
        try:
            prepare_config(treatment_names=list(combo))
            metrics = run_simulation(label)
            if metrics:
                os.makedirs(OUTPUT_BASE, exist_ok=True)
                safe_name = label.replace("+", "_")
                shutil.copy2(metrics, os.path.join(OUTPUT_BASE, f"metrics_{safe_name}.csv"))
            outcomes = extract_outcomes(metrics)
            results.append((label, outcomes))
            closure = outcomes.get("wound_closure_pct", 0)
            base_closure = results[0][1].get("wound_closure_pct", 0) if results[0][1] else 0
            delta = closure - base_closure
            print(f"    >> Closure: {closure:.1f}% ({delta:+.1f}% vs baseline)", flush=True)
        except Exception as e:
            print(f"    ERROR: {e}", flush=True)
            results.append((label, None))
        print(flush=True)
        run_idx += 1

    # Compare
    print_comparison(results)
    save_results(results, OUTPUT_BASE)

    # Run validation on baseline
    print("\n--- Validation (baseline diabetic) ---")
    baseline_csv = os.path.join(OUTPUT_BASE, "metrics_baseline.csv")
    if os.path.exists(baseline_csv):
        subprocess.run(
            [sys.executable, os.path.join(REPO, "literature", "validate_all.py"),
             baseline_csv],
            cwd=REPO,
        )


if __name__ == "__main__":
    main()
