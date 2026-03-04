#!/usr/bin/env python3
"""Diabetic wound treatment study.

Runs a diabetic wound baseline plus therapeutic interventions, then generates
an Excel workbook comparing outcomes.

Usage:
    python3 studies/diabetic-wound/treatment.py                        # all treatments
    python3 studies/diabetic-wound/treatment.py --quick                # baseline + 2
    python3 studies/diabetic-wound/treatment.py --treatments hbo,npwt  # specific ones
    python3 studies/diabetic-wound/treatment.py --combos pairs         # pairwise combos
    python3 studies/diabetic-wound/treatment.py --no-excel             # skip Excel
"""

import argparse
import glob
import os
import subprocess
import sys
import time

REPO = os.path.normpath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
sys.path.insert(0, REPO)

from scripts.study.study_common import check_bdm, build_if_needed, banner, summary


def main():
    parser = argparse.ArgumentParser(description="Diabetic wound treatment study")
    parser.add_argument("--quick", action="store_true",
                        help="Run baseline + 2 treatments only")
    parser.add_argument("--treatments", default="all",
                        help="Comma-separated treatment names or 'all'")
    parser.add_argument("--combos", default=None,
                        choices=["pairs", "triples", "all"],
                        help="Combination mode")
    parser.add_argument("--workers", type=int, default=None,
                        help="Parallel workers")
    parser.add_argument("--no-excel", action="store_true",
                        help="Skip Excel workbook generation")
    args = parser.parse_args()

    treatments = "growth_factor,npwt" if args.quick else args.treatments

    check_bdm()
    build_if_needed()

    # Count treatments for banner
    if treatments == "all":
        n_treat = len(glob.glob(os.path.join(REPO, "treatments", "*.toml")))
    else:
        n_treat = len(treatments.split(","))

    detail = [f"{n_treat} treatments (+ baseline)"]
    if args.combos:
        detail.append(f"Combinations: {args.combos}")
    banner("Diabetic Treatment Study", detail)

    t0 = time.time()

    # Run treatment study
    cmd = [
        sys.executable, "-u",
        os.path.join(REPO, "scripts", "study", "treatment_study.py"),
        f"--treatments={treatments}",
    ]
    if args.combos:
        cmd.append(f"--combos={args.combos}")
    if args.workers:
        cmd.append(f"--workers={args.workers}")
    subprocess.run(cmd, cwd=REPO, check=True)

    src = os.path.join(REPO, "output", "treatment_study")

    # Generate Excel workbook
    if not args.no_excel:
        try:
            import openpyxl  # noqa: F401
            print()
            print("Generating Excel workbook...")
            subprocess.run([
                sys.executable,
                os.path.join(REPO, "scripts", "study", "generate_study.py"),
                src, "-o", src, "--name", "diabetic_study",
            ], cwd=REPO, check=False)
        except ImportError:
            print()
            print("Skipping Excel (pip install openpyxl to enable)")

    # Generate baseline plots
    baseline_csv = os.path.join(src, "metrics_baseline.csv")
    if os.path.isfile(baseline_csv):
        print()
        print("Generating baseline plots...")
        plot_dir = os.path.join(src, "plots")
        os.makedirs(plot_dir, exist_ok=True)
        subprocess.run([
            sys.executable,
            os.path.join(REPO, "scripts", "analysis", "plot_metrics.py"),
            baseline_csv, "-o", plot_dir,
        ], cwd=REPO, check=False)

    # Resilience analysis
    resilience_script = os.path.join(REPO, "scripts", "analysis", "resilience_analysis.py")
    if os.path.isfile(resilience_script):
        print()
        print("Running resilience analysis...")
        subprocess.run([sys.executable, resilience_script, src],
                       cwd=REPO, check=False)

    summary(src, t0, "Study complete")


if __name__ == "__main__":
    main()
