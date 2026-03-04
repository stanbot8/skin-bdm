#!/usr/bin/env python3
"""Adaptive surrogate-guided treatment combination search.

Efficiently explores 2^N-1 treatment combinations using adaptive consensus
(early stopping when variance is low), an additive surrogate model (predict
combos from singles), and guided selection (only run promising combos).

Usage:
    python3 studies/diabetic-wound/adaptive.py                      # full search
    python3 studies/diabetic-wound/adaptive.py --quick               # top-5, 2-3 runs
    python3 studies/diabetic-wound/adaptive.py --top-k 10            # top-10 combos
    python3 studies/diabetic-wound/adaptive.py --treatments hbo,npwt,growth_factor
"""

import argparse
import glob
import os
import shutil
import subprocess
import sys
import time

REPO = os.path.normpath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
sys.path.insert(0, REPO)

from scripts.study.study_common import check_bdm, build_if_needed, banner, summary


def main():
    parser = argparse.ArgumentParser(
        description="Adaptive surrogate-guided combination search")
    parser.add_argument("--quick", action="store_true",
                        help="Top-5, 2-3 runs per combo")
    parser.add_argument("--top-k", type=int, default=20,
                        help="Number of top combos to evaluate (default: 20)")
    parser.add_argument("--min-runs", type=int, default=3)
    parser.add_argument("--max-runs", type=int, default=10)
    parser.add_argument("--cv-threshold", type=float, default=0.05,
                        help="CV threshold for early stopping (default: 0.05)")
    parser.add_argument("--treatments", default="all",
                        help="Comma-separated treatment names or 'all'")
    parser.add_argument("--no-excel", action="store_true",
                        help="Skip Excel workbook generation")
    args = parser.parse_args()

    if args.quick:
        args.top_k = 5
        args.min_runs = 2
        args.max_runs = 3

    check_bdm()
    build_if_needed()

    banner("Adaptive Combination Search", [
        f"Top-K: {args.top_k}, Runs: {args.min_runs}-{args.max_runs}",
        f"CV threshold: {args.cv_threshold}",
    ])

    t0 = time.time()

    subprocess.run([
        sys.executable, "-u",
        os.path.join(REPO, "scripts", "study", "adaptive_study.py"),
        f"--treatments={args.treatments}",
        f"--top-k={args.top_k}",
        f"--min-runs={args.min_runs}",
        f"--max-runs={args.max_runs}",
        f"--cv-threshold={args.cv_threshold}",
    ], cwd=REPO, check=True)

    # Collect results
    src = os.path.join(REPO, "output", "adaptive_study")
    results = os.path.join(os.path.dirname(__file__), "results")
    os.makedirs(results, exist_ok=True)
    for csv_file in glob.glob(os.path.join(src, "*.csv")):
        shutil.copy2(csv_file, results)

    # Generate Excel workbook
    if not args.no_excel:
        try:
            import openpyxl  # noqa: F401
            print()
            print("Generating Excel workbook...")
            subprocess.run([
                sys.executable,
                os.path.join(REPO, "scripts", "study", "generate_study.py"),
                src, "-o", results, "--name", "adaptive_study",
            ], cwd=REPO, check=False)
        except ImportError:
            print()
            print("Skipping Excel (pip install openpyxl to enable)")

    summary(results, t0, "Adaptive study complete")


if __name__ == "__main__":
    main()
