#!/usr/bin/env python3
"""Tumor growth study.

Runs tumor and tumor+wound experiments with consensus batches, then generates
comparison plots and validation.

Usage:
    python3 studies/tumor/study.py              # 10-run consensus each
    python3 studies/tumor/study.py --quick      # 5-run consensus
    python3 studies/tumor/study.py -n 20        # 20-run consensus
"""

import argparse
import os
import sys
import time

REPO = os.path.normpath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
sys.path.insert(0, REPO)

from scripts.study.study_common import (
    check_bdm, build_if_needed, banner, run_batch, copy_results,
    plot_metrics, compare_profiles, summary,
)

EXPERIMENTS = [
    {"label": "tumor",       "skin": "normal", "study": "tumor",
     "desc": "Tumor only (normal skin)"},
    {"label": "tumor_wound", "skin": "normal", "study": "tumor-wound",
     "desc": "Tumor with wound (normal skin)"},
    {"label": "aged_tumor",  "skin": "aged",   "study": "tumor",
     "desc": "Tumor with aged skin"},
]


def main():
    parser = argparse.ArgumentParser(description="Tumor growth study")
    parser.add_argument("-n", "--num-runs", type=int, default=10)
    parser.add_argument("--quick", action="store_true", help="5-run consensus")
    args = parser.parse_args()

    n_runs = 5 if args.quick else args.num_runs

    check_bdm()
    build_if_needed()

    results = os.path.join(os.path.dirname(__file__), "results")
    os.makedirs(results, exist_ok=True)

    banner("Tumor Growth Study", [f"{n_runs}-run consensus per experiment"])

    t0 = time.time()

    for i, exp in enumerate(EXPERIMENTS, 1):
        print(f"[{i}/{len(EXPERIMENTS)}] {exp['desc']}")
        consensus_dir = run_batch(n_runs, exp["skin"], exp["study"], validate=True)
        copy_results(consensus_dir, results, exp["label"])
        print()

    plot_metrics(results)
    compare_profiles(results, [e["label"] for e in EXPERIMENTS])
    summary(results, t0, "Tumor study complete")


if __name__ == "__main__":
    main()
