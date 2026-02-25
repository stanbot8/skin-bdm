#!/usr/bin/env python3
"""Treatment study: run diabetic wound with each treatment and compare outcomes.

Usage:
    python3 scripts/study/treatment_study.py [--treatments=all|name1,name2,...] [--steps=N]
    python3 scripts/study/treatment_study.py --combos=pairs   # singles + all pairwise
    python3 scripts/study/treatment_study.py --combos=all     # all 2^N-1 subsets
    python3 scripts/study/treatment_study.py --workers=8      # parallel workers
"""

import argparse
import concurrent.futures
import csv
import itertools
import os
import shutil
import subprocess
import sys
import tempfile
import time

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BINARY = os.path.join(REPO, "build", "skibidy")
OUTPUT_BASE = os.path.join(REPO, "output", "treatment_study")
TREATMENTS_DIR = os.path.join(REPO, "treatments")
MERGE_SCRIPT = os.path.join(REPO, "scripts", "config", "merge_config.py")
APPLY_SCRIPT = os.path.join(REPO, "scripts", "config", "apply_preset.py")
DIABETIC_PROFILE = os.path.join(REPO, "profiles", "diabetic.toml")
STUDY_PRESET = os.path.join(REPO, "studies", "diabetic-wound", "preset.toml")
_num_workers = 2  # set by main() before launching jobs


def available_treatments():
    """List all treatment TOML files."""
    tdir = TREATMENTS_DIR
    return sorted(
        os.path.splitext(f)[0]
        for f in os.listdir(tdir)
        if f.endswith(".toml")
    )


def _patch_config(bdm_path):
    """Strip visualization sections and disable metrics_autoopen."""
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


def run_job(args):
    """Run one simulation in an isolated temp directory.

    Args:
        args: (label, treatment_names, safe_name)

    Returns:
        (label, outcomes, elapsed_s)
    """
    label, treatment_names, safe_name = args

    with tempfile.TemporaryDirectory(prefix="skibidy_") as work_dir:
        bdm_path = os.path.join(work_dir, "bdm.toml")

        # Generate merged config into isolated work dir
        r = subprocess.run(
            [sys.executable, MERGE_SCRIPT,
             os.path.join(REPO, "bdm.core.toml"),
             os.path.join(REPO, "modules"),
             bdm_path],
            cwd=REPO, capture_output=True, text=True)
        if r.returncode != 0:
            return label, {}, 0.0

        # Apply diabetic profile, study preset, then each treatment
        presets = [DIABETIC_PROFILE, STUDY_PRESET]
        for t in (treatment_names or []):
            presets.append(os.path.join(TREATMENTS_DIR, f"{t}.toml"))
        for preset in presets:
            r = subprocess.run(
                [sys.executable, APPLY_SCRIPT, preset, bdm_path],
                cwd=REPO, capture_output=True, text=True)
            if r.returncode != 0:
                return label, {}, 0.0

        # Patch config for headless batch (strip viz, disable autoopen)
        _patch_config(bdm_path)

        # Run binary from isolated work dir (partition cores across workers)
        t0 = time.time()
        env = os.environ.copy()
        env.pop("DISPLAY", None)
        cores = os.cpu_count() or 4
        env["OMP_NUM_THREADS"] = str(max(1, cores // _num_workers))
        subprocess.run([BINARY], cwd=work_dir, env=env, capture_output=True)
        elapsed = time.time() - t0

        # Copy metrics to output dir before tempdir is cleaned up
        metrics_src = os.path.join(work_dir, "output", "skibidy", "metrics.csv")
        saved_path = None
        if os.path.exists(metrics_src):
            os.makedirs(OUTPUT_BASE, exist_ok=True)
            saved_path = os.path.join(OUTPUT_BASE, f"metrics_{safe_name}.csv")
            shutil.copy2(metrics_src, saved_path)

        outcomes = extract_outcomes(saved_path)
        return label, outcomes, elapsed


def _find_biological_rows(rows):
    """Return rows up to (but not including) the dissolution stamp.

    WoundResolution has two effects that can appear in different metrics
    intervals: (1) stratum stamp fills remaining wound voxels, causing a
    large closure jump, and (2) agent dissolution drops n_agents to 0.
    The stamp always precedes or coincides with dissolution.

    We detect the earliest of:
      - n_agents dropping to 0 (definitive dissolution)
      - closure jumping >5% in one interval while n_agents also drops
        (stratum stamp with partial dissolution in same interval)
    And return all rows before that point.
    """
    cutoff = len(rows)
    for i in range(len(rows)):
        if int(float(rows[i].get("n_agents", 1))) == 0:
            cutoff = i
            break
    # Check if the stratum stamp happened one interval earlier
    if cutoff > 0:
        prev_closure = float(rows[cutoff - 1].get("wound_closure_pct", 0))
        # Walk back to find the stamp: closure jumped while agents were
        # still present (stratum fill happens before agent removal)
        for i in range(1, cutoff):
            prev = float(rows[i - 1].get("wound_closure_pct", 0))
            curr = float(rows[i].get("wound_closure_pct", 0))
            if curr - prev > 5.0 and curr >= 95.0:
                cutoff = i
                break
    return rows[:cutoff] if cutoff > 0 else rows[:1]


def extract_outcomes(metrics_path):
    """Extract key outcome metrics from a simulation run."""
    if not metrics_path or not os.path.exists(metrics_path):
        return {}

    with open(metrics_path) as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    if not rows:
        return {}

    # Use pre-dissolution rows for closure and time-to-threshold metrics
    bio_rows = _find_biological_rows(rows)

    outcomes = {}

    bio_last = bio_rows[-1]
    outcomes["wound_closure_pct"] = float(bio_last.get("wound_closure_pct", 0))

    # Final inflammation, scar, and agent count from actual last row
    last = rows[-1]
    outcomes["mean_infl_wound"] = float(last.get("mean_infl_wound", 0))
    outcomes["scar_magnitude"] = float(last.get("scar_magnitude", 0))
    outcomes["n_agents"] = int(float(last.get("n_agents", 0)))

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

    # Time-to-threshold uses biological rows only (no dissolution stamp)
    for row in bio_rows:
        if float(row.get("wound_closure_pct", 0)) >= 50:
            step = int(float(row.get("step", 0)))
            outcomes["time_to_50pct_days"] = round(step * 0.1 / 24, 1)
            break
    else:
        outcomes["time_to_50pct_days"] = None

    for row in bio_rows:
        if float(row.get("wound_closure_pct", 0)) >= 90:
            step = int(float(row.get("step", 0)))
            outcomes["time_to_90pct_days"] = round(step * 0.1 / 24, 1)
            break
    else:
        outcomes["time_to_90pct_days"] = None

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

    peak_col = 0
    for row in rows:
        col = float(row.get("mean_collagen_wound", 0))
        if col > peak_col:
            peak_col = col
    outcomes["peak_collagen"] = peak_col

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
            delta_t50 = f"{t50_base - t50_treat:+.1f}d" if t50_base and t50_treat else "N/A"

            t90_base = baseline_outcomes.get("time_to_90pct_days")
            t90_treat = outcomes.get("time_to_90pct_days")
            delta_t90 = f"{t90_base - t90_treat:+.1f}d" if t90_base and t90_treat else "N/A"

            peak_infl_base = baseline_outcomes.get("peak_inflammation", 0)
            peak_infl_treat = outcomes.get("peak_inflammation", 0)
            delta_infl = (f"{(peak_infl_treat - peak_infl_base) / peak_infl_base * 100:+.0f}%"
                          if peak_infl_base > 0 else "N/A")

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
                        help="Override num_steps (default: from study config)")
    parser.add_argument("--workers", type=int,
                        default=min(2, os.cpu_count() or 2),
                        help="Parallel simulations (default: min(2, cpu_count))")
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

    combo_treatments = [t for t in treatments if t != "combination"]
    combos = generate_combos(combo_treatments, args.combos) if args.combos else []

    # Build ordered job list: (label, treatment_names, safe_name)
    jobs = [("baseline", None, "baseline")]
    for t in treatments:
        jobs.append((t, [t], t))
    for combo in combos:
        lbl = combo_label(combo)
        jobs.append((lbl, list(combo), lbl.replace("+", "_")))

    total = len(jobs)
    print(f"Treatment study: {total} simulations ({args.workers} workers)", flush=True)
    print(f"  1 baseline + {len(treatments)} singles", end="", flush=True)
    if combos:
        print(f" + {len(combos)} combinations ({args.combos})", end="", flush=True)
    print(flush=True)
    print(flush=True)

    global _num_workers
    _num_workers = args.workers

    results_dict = {}
    done = 0

    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as executor:
        future_to_label = {executor.submit(run_job, job): job[0] for job in jobs}
        for future in concurrent.futures.as_completed(future_to_label):
            label, outcomes, elapsed = future.result()
            done += 1
            mins = int(elapsed) // 60
            secs = int(elapsed) % 60
            closure = outcomes.get("wound_closure_pct", 0) if outcomes else 0
            base_outcomes = results_dict.get("baseline")
            if label != "baseline" and base_outcomes is not None:
                base_closure = base_outcomes.get("wound_closure_pct", 0)
                delta = f" ({closure - base_closure:+.1f}% vs baseline)"
            else:
                delta = ""
            status = "OK" if outcomes else "FAILED"
            print(f"  [{done}/{total}] {label}: closure {closure:.1f}%{delta} "
                  f"({mins}m{secs}s) {status}", flush=True)
            results_dict[label] = outcomes

    # Reconstruct results in original job order
    results = [(label, results_dict.get(label, {})) for label, _, _ in jobs]

    print_comparison(results)
    save_results(results, OUTPUT_BASE)

    # Validation on baseline
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
