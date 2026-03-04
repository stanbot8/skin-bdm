#!/usr/bin/env python3
"""Generate treatment combinatorial experiment TOMLs.

Three-phase design for efficient exploration:

  Phase 1 (--phase screen): All 2^N-1 treatment subsets at day 0.
    511 configs for 9 treatments, ~12.8h. Answers: which combos help?

  Phase 2 (--phase refine): Read screen results, pick top N combos by T50,
    generate timing grids only for those. Efficient: explores 4-8 way combos
    without 5^8 explosion.

  Phase 3 (--phase timing): Brute-force timing grid for all combos up to
    --combo-max. Use only for small combo sizes (1-3).

Usage:
    python3 scripts/study/gen_treatment_schedule.py --phase screen
    python3 scripts/study/gen_treatment_schedule.py --phase refine --top 20
    python3 scripts/study/gen_treatment_schedule.py --phase refine --top 10 --days 0,7,14
    python3 scripts/study/gen_treatment_schedule.py --phase timing --combo-max 2

Output: studies/diabetic-wound/experiments/
"""

import argparse
import csv
import itertools
import os
import sys

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    import tomllib
except ModuleNotFoundError:
    try:
        import tomli as tomllib
    except ModuleNotFoundError:
        import pip._vendor.tomli as tomllib

ALL_TREATMENTS = [
    "npwt", "hbo", "growth_factor", "doxycycline",
    "anti_inflammatory", "msc", "moisture", "combination", "senolytic",
]


def load_treatment(name, study="diabetic-wound"):
    """Load a treatment TOML and return its parameter overrides as flat dict.

    Searches study-scoped treatments first, then shared, then all studies.
    """
    import glob as _glob
    candidates = [
        os.path.join(REPO, "studies", study, "treatments", f"{name}.toml"),
        os.path.join(REPO, "studies", "shared", "treatments", f"{name}.toml"),
    ]
    for c in candidates:
        if os.path.isfile(c):
            with open(c, "rb") as f:
                return flatten_toml(tomllib.load(f))
    for c in sorted(_glob.glob(os.path.join(
            REPO, "studies", "*", "treatments", f"{name}.toml"))):
        with open(c, "rb") as f:
            return flatten_toml(tomllib.load(f))
    print(f"ERROR: treatment not found: {name}")
    sys.exit(1)


def load_profile(name):
    """Load a profile TOML (diabetic baseline) and return flat dict."""
    path = os.path.join(REPO, "profiles", f"{name}.toml")
    if not os.path.isfile(path):
        return {}
    with open(path, "rb") as f:
        data = tomllib.load(f)
    return flatten_toml(data)


def flatten_toml(d, prefix=""):
    """Flatten nested TOML dict to dotted keys."""
    flat = {}
    for k, v in d.items():
        key = f"{prefix}{k}" if prefix else k
        if isinstance(v, dict):
            flat.update(flatten_toml(v, key + "."))
        else:
            flat[key] = v
    return flat


def interpolate(baseline, treatment, fraction):
    """Interpolate between baseline and treatment values.

    fraction=1.0 means full treatment, 0.0 means no treatment.
    Only interpolates numeric params that differ.
    """
    result = {}
    for key, treat_val in treatment.items():
        base_val = baseline.get(key)
        if base_val is None:
            if fraction > 0:
                result[key] = treat_val
            continue
        if isinstance(treat_val, (int, float)) and isinstance(base_val, (int, float)):
            interp = base_val + fraction * (treat_val - base_val)
            if abs(interp - base_val) > 1e-10:
                result[key] = round(interp, 6)
        elif isinstance(treat_val, bool) and isinstance(base_val, bool):
            result[key] = treat_val if fraction >= 0.5 else base_val
        elif isinstance(treat_val, str):
            result[key] = treat_val if fraction > 0 else base_val
    return result


def treatment_fraction(start_day, sim_days=42):
    """Fraction of treatment effect based on start day."""
    if start_day >= sim_days:
        return 0.0
    return (sim_days - start_day) / sim_days


def _merge_overrides(override_list, baseline):
    """Merge multiple override dicts. For overlapping keys, take value furthest from baseline."""
    merged = {}
    for overrides in override_list:
        for key, val in overrides.items():
            if key in merged:
                base_val = baseline.get(key, 0)
                if isinstance(val, (int, float)) and isinstance(base_val, (int, float)):
                    if abs(val - base_val) > abs(merged[key] - base_val):
                        merged[key] = val
            else:
                merged[key] = val
    return merged


def _format_override(key, val):
    """Format a single override line for TOML output."""
    if isinstance(val, str):
        return f'"{key}" = "{val}"'
    elif isinstance(val, bool):
        return f'"{key}" = {"true" if val else "false"}'
    elif isinstance(val, float):
        return f'"{key}" = {val}'
    elif isinstance(val, int):
        return f'"{key}" = {val}'
    return None


# ---------------------------------------------------------------------------
# Phase 1: Screening (all subsets at day 0)
# ---------------------------------------------------------------------------

def generate_screen_experiment(treatments, baseline, runs_per_config=3):
    """Generate one experiment TOML with all 2^N-1 treatment subsets at day 0."""
    lines = []
    lines.append(f"# Auto-generated: full combinatorial screen ({len(treatments)} treatments)")
    lines.append(f"# All subsets tested at day 0 (full treatment effect)")
    lines.append("")
    lines.append("[experiment]")
    lines.append(f'name = "Combinatorial Screen"')
    lines.append(f'description = "All {2**len(treatments)-1} treatment subsets at day 0"')
    lines.append('profile = "diabetic"')
    lines.append('study = "diabetic-wound"')
    lines.append(f"runs_per_config = {runs_per_config}")
    lines.append("")

    # Untreated baseline
    lines.append("[[experiment.configs]]")
    lines.append('label = "Untreated"')
    lines.append("")

    # All non-empty subsets, ordered by size then alphabetically
    for k in range(1, len(treatments) + 1):
        for group in itertools.combinations(treatments, k):
            group = list(group)
            tlist = ", ".join(f'"{t}"' for t in group)
            label = " + ".join(t.upper() for t in group)

            lines.append("[[experiment.configs]]")
            lines.append(f'label = "{label}"')
            lines.append(f"treatments = [{tlist}]")
            lines.append("")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Phase 2: Timing grid (cross-product of start days)
# ---------------------------------------------------------------------------

def generate_timing_experiment(treatments, days_matrix, baseline, runs_per_config=3):
    """Generate timing experiment for a specific treatment combo."""
    combo_name = "+".join(treatments)

    lines = []
    lines.append(f"# Auto-generated: timing grid ({combo_name})")
    lines.append("")
    lines.append("[experiment]")
    lines.append(f'name = "Timing: {combo_name}"')
    lines.append(f'description = "Start-day grid for {combo_name}"')
    lines.append('profile = "diabetic"')
    lines.append('study = "diabetic-wound"')
    lines.append(f"runs_per_config = {runs_per_config}")
    lines.append("")

    # Untreated baseline
    lines.append("[[experiment.configs]]")
    lines.append('label = "Untreated"')
    lines.append("")

    if len(treatments) == 1:
        # Single treatment: simple day sweep
        t = treatments[0]
        treatment = load_treatment(t)
        for day in days_matrix:
            frac = treatment_fraction(day)
            label = f"{t.upper()} Day {day} ({frac*100:.0f}%)"
            overrides = interpolate(baseline, treatment, frac)

            lines.append("[[experiment.configs]]")
            lines.append(f'label = "{label}"')
            if day == 0:
                lines.append(f'treatments = ["{t}"]')
            elif overrides:
                lines.append("[experiment.configs.overrides]")
                for key, val in sorted(overrides.items()):
                    fmt = _format_override(key, val)
                    if fmt:
                        lines.append(fmt)
            lines.append("")
    else:
        # Multi-treatment: cross-product of start days
        for day_combo in itertools.product(days_matrix, repeat=len(treatments)):
            override_list = []
            label_parts = []
            all_day0 = True
            for t, d in zip(treatments, day_combo):
                frac = treatment_fraction(d)
                treat = load_treatment(t)
                override_list.append(interpolate(baseline, treat, frac))
                label_parts.append(f"{t.upper()} D{d}")
                if d != 0:
                    all_day0 = False

            merged = _merge_overrides(override_list, baseline)
            label = " + ".join(label_parts)

            lines.append("[[experiment.configs]]")
            lines.append(f'label = "{label}"')
            if all_day0:
                tlist = ", ".join(f'"{t}"' for t in treatments)
                lines.append(f"treatments = [{tlist}]")
            elif merged:
                lines.append("[experiment.configs.overrides]")
                for key, val in sorted(merged.items()):
                    fmt = _format_override(key, val)
                    if fmt:
                        lines.append(fmt)
            lines.append("")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Phase 3: Refine (read screen results, pick top N, generate timing)
# ---------------------------------------------------------------------------

def parse_screen_results(results_dir):
    """Read comparison.csv from screen run and return ranked combos.

    Returns list of (label, treatments_list, t50_days, closure_pct, scar) sorted by T50.
    """
    path = os.path.join(results_dir, "comparison.csv")
    if not os.path.isfile(path):
        print(f"ERROR: screen results not found: {path}")
        print("Run the screen phase first:")
        print("  python3 studies/run_experiments.py studies/diabetic-wound/experiments/combo_screen.toml")
        sys.exit(1)

    results = []
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            label = row.get("config", "")
            if label == "Untreated":
                continue

            # Parse treatment names from label ("NPWT + HBO + MSC" -> ["npwt", "hbo", "msc"])
            treatments = [t.strip().lower() for t in label.split(" + ")]

            t50 = row.get("time_to_50pct_days", "")
            closure = row.get("wound_closure_pct", "")
            scar = row.get("scar_magnitude", "")

            try:
                t50 = float(t50) if t50 else 999
            except ValueError:
                t50 = 999
            try:
                closure = float(closure) if closure else 0
            except ValueError:
                closure = 0
            try:
                scar = float(scar) if scar else 999
            except ValueError:
                scar = 999

            results.append((label, treatments, t50, closure, scar))

    # Sort by T50 (lower is better), then by scar (lower is better)
    results.sort(key=lambda x: (x[2], x[4]))
    return results


def print_screen_ranking(results, top_n=20):
    """Print ranked screen results."""
    print(f"\n  {'Rank':>4}  {'Combo':<50} {'T50 (d)':>8} {'Closure':>8} {'Scar':>8} {'Size':>4}")
    print(f"  {'-'*86}")
    for i, (label, treatments, t50, closure, scar) in enumerate(results[:top_n]):
        t50_str = f"{t50:.1f}" if t50 < 999 else "N/A"
        print(f"  {i+1:>4}  {label:<50} {t50_str:>8} {closure:>7.1f}% {scar:>8.3f} {len(treatments):>4}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Generate treatment combinatorial experiment TOMLs")
    parser.add_argument("--treatments", type=str, default=",".join(ALL_TREATMENTS),
                        help="Comma-separated treatment names (default: all 9)")
    parser.add_argument("--phase", type=str, default="screen",
                        choices=["screen", "refine", "timing"],
                        help="Phase: screen (all subsets day 0), refine (top N from screen), "
                             "timing (brute-force grid)")
    parser.add_argument("--days", type=str, default="0,3,7,14,21",
                        help="Start days for timing/refine (default: 0,3,7,14,21)")
    parser.add_argument("--top", type=int, default=20,
                        help="Refine phase: number of top combos to explore (default: 20)")
    parser.add_argument("--combo-min", type=int, default=1, dest="combo_min",
                        help="Timing phase: min combo size (default: 1)")
    parser.add_argument("--combo-max", type=int, default=None, dest="combo_max",
                        help="Timing phase: max combo size (default: all)")
    parser.add_argument("--screen-dir", type=str, default=None, dest="screen_dir",
                        help="Refine phase: path to screen results dir "
                             "(default: studies/diabetic-wound/results/experiments/combinatorial_screen/)")
    parser.add_argument("--runs", type=int, default=3,
                        help="Runs per config (default: 3)")
    parser.add_argument("--dry-run", action="store_true", dest="dry_run",
                        help="Print config counts without writing files")
    parser.add_argument("--output-dir", type=str,
                        default=os.path.join(REPO, "studies", "diabetic-wound", "experiments"))
    args = parser.parse_args()

    treatments = [t.strip() for t in args.treatments.split(",")]
    days = [int(d.strip()) for d in args.days.split(",")]
    baseline = load_profile("diabetic")

    os.makedirs(args.output_dir, exist_ok=True)

    total_configs = 0
    total_files = 0

    # ---- Phase: Screen ----
    if args.phase == "screen":
        n_subsets = 2 ** len(treatments) - 1
        n_configs = n_subsets + 1
        total_configs += n_configs
        total_files += 1

        if args.dry_run:
            print(f"  [screen] combo_screen.toml: {n_configs} configs x {args.runs} runs "
                  f"= {n_configs * args.runs} sims")
        else:
            content = generate_screen_experiment(treatments, baseline, args.runs)
            path = os.path.join(args.output_dir, "combo_screen.toml")
            with open(path, "w") as f:
                f.write(content)
            print(f"  Generated: {path} ({n_configs} configs)")

    # ---- Phase: Refine ----
    elif args.phase == "refine":
        screen_dir = args.screen_dir or os.path.join(
            REPO, "studies", "diabetic-wound", "results", "experiments", "combinatorial_screen")
        ranked = parse_screen_results(screen_dir)

        if not ranked:
            print("ERROR: no results found in screen output")
            sys.exit(1)

        print(f"  Screen results: {len(ranked)} combos ranked")
        print_screen_ranking(ranked, args.top)

        # Generate timing grids for top N
        selected = ranked[:args.top]
        print(f"\n  Generating timing grids for top {len(selected)} combos:\n")

        for i, (label, treat_list, t50, closure, scar) in enumerate(selected):
            k = len(treat_list)
            n_configs = len(days) ** k + 1
            total_configs += n_configs
            total_files += 1

            name = f"refine_{'_'.join(treat_list)}.toml"

            if args.dry_run:
                print(f"  [{i+1:>2}] {name}: {n_configs} configs x {args.runs} runs "
                      f"= {n_configs * args.runs} sims (size {k})")
                continue

            content = generate_timing_experiment(treat_list, days, baseline, args.runs)
            path = os.path.join(args.output_dir, name)
            with open(path, "w") as f:
                f.write(content)
            print(f"  Generated: {path} ({n_configs} configs, size {k})")

    # ---- Phase: Timing (brute-force) ----
    elif args.phase == "timing":
        min_k = max(1, args.combo_min)
        max_k = len(treatments) if args.combo_max is None else min(args.combo_max, len(treatments))

        for k in range(min_k, max_k + 1):
            for group in itertools.combinations(treatments, k):
                group = list(group)
                n_configs = len(days) ** k + 1
                total_configs += n_configs
                total_files += 1

                name = f"timing_{'_'.join(group)}.toml"

                if args.dry_run:
                    print(f"  [timing] {name}: {n_configs} configs x {args.runs} runs "
                          f"= {n_configs * args.runs} sims")
                    continue

                content = generate_timing_experiment(group, days, baseline, args.runs)
                path = os.path.join(args.output_dir, name)
                with open(path, "w") as f:
                    f.write(content)
                print(f"  Generated: {path} ({n_configs} configs)")

    total_runs = total_configs * args.runs
    est_hours = total_runs * 30 / 3600
    print(f"\n  Total: {total_files} files, {total_configs} configs, "
          f"{total_runs} runs (~{est_hours:.1f}h at 30s/run)")

    if not args.dry_run and total_files > 0:
        print(f"\nRun with:")
        if args.phase == "screen":
            print(f"  python3 studies/run_experiments.py studies/diabetic-wound/experiments/combo_screen.toml")
        elif args.phase == "refine":
            print(f"  python3 studies/run_experiments.py studies/diabetic-wound/experiments/refine_*.toml")
        elif args.phase == "timing":
            print(f"  python3 studies/run_experiments.py studies/diabetic-wound/experiments/timing_*.toml")


if __name__ == "__main__":
    main()
