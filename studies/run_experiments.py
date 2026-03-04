#!/usr/bin/env python3
"""Run study experiments from TOML configs.

Usage:
    python3 studies/run_experiments.py                                              # all experiments
    python3 studies/run_experiments.py studies/diabetic-wound/experiments/wound_size.toml  # one file
    python3 studies/run_experiments.py --quick                                      # all with 2 runs
    python3 studies/run_experiments.py --runs=3 studies/diabetic-wound/experiments/biofilm*
"""

import argparse
import glob
import os
import subprocess
import sys
import time

REPO = os.path.normpath(os.path.join(os.path.dirname(__file__), os.pardir))
sys.path.insert(0, REPO)

from scripts.study.study_common import check_bdm, build_if_needed, banner


def main():
    parser = argparse.ArgumentParser(description="Run study experiments")
    parser.add_argument("experiments", nargs="*",
                        help="TOML experiment files (default: all)")
    parser.add_argument("--quick", action="store_true",
                        help="Use 2 runs per experiment")
    parser.add_argument("--runs", type=int, default=None,
                        help="Runs per experiment config")
    parser.add_argument("--list", action="store_true",
                        help="List available experiments and exit")
    args = parser.parse_args()

    if args.list:
        for f in sorted(glob.glob(os.path.join(REPO, "studies", "*", "experiments", "*.toml"))):
            study = os.path.basename(os.path.dirname(os.path.dirname(f)))
            name = os.path.splitext(os.path.basename(f))[0]
            print(f"  {study:12s} {name}")
        return

    check_bdm()
    build_if_needed()

    experiments = args.experiments
    if not experiments:
        experiments = sorted(glob.glob(
            os.path.join(REPO, "studies", "*", "experiments", "*.toml")))

    if not experiments:
        print("No experiment TOMLs found.")
        return

    runs_flag = []
    if args.quick:
        runs_flag = ["--runs=2"]
    elif args.runs:
        runs_flag = [f"--runs={args.runs}"]

    runner = os.path.join(REPO, "scripts", "study", "experiment_runner.py")

    banner("Experiment Runner", [f"{len(experiments)} experiment(s) to run"])

    t0 = time.time()
    for exp in experiments:
        name = os.path.splitext(os.path.basename(exp))[0]
        print(f">>> {name}")
        subprocess.run(
            [sys.executable, "-u", runner, exp] + runs_flag,
            cwd=REPO, check=True)
        print()

    elapsed = time.time() - t0
    mins = int(elapsed) // 60
    secs = int(elapsed) % 60
    print("=" * 44)
    print(f"  All experiments completed in {mins}m {secs}s")
    print("=" * 44)


if __name__ == "__main__":
    main()
