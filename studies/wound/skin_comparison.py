#!/usr/bin/env python3
"""Skin type comparison study.

Runs wound healing with normal, aged, and diabetic skin profiles using
consensus batches, then compares healing trajectories.

Usage:
    python3 studies/wound/skin_comparison.py             # 10-run consensus
    python3 studies/wound/skin_comparison.py --quick      # 5-run consensus
    python3 studies/wound/skin_comparison.py -n 20        # 20-run consensus
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

PROFILES = [
    {"label": "normal",   "skin": "normal",   "study": "wound"},
    {"label": "aged",     "skin": "aged",     "study": "wound"},
    {"label": "diabetic", "skin": "diabetic", "study": "diabetic-wound"},
]


def main():
    parser = argparse.ArgumentParser(description="Skin type comparison study")
    parser.add_argument("-n", "--num-runs", type=int, default=10)
    parser.add_argument("--quick", action="store_true", help="5-run consensus")
    args = parser.parse_args()

    n_runs = 5 if args.quick else args.num_runs

    check_bdm()
    build_if_needed()

    results = os.path.join(os.path.dirname(__file__), "results")
    os.makedirs(results, exist_ok=True)

    banner("Skin Type Comparison Study", [f"{n_runs}-run consensus per profile"])

    t0 = time.time()

    for i, p in enumerate(PROFILES, 1):
        print(f"[{i}/{len(PROFILES)}] {p['skin'].capitalize()} skin + wound")
        consensus_dir = run_batch(n_runs, p["skin"], p["study"], validate=True)
        copy_results(consensus_dir, results, p["label"])
        print()

    plot_metrics(results)
    compare_profiles(results, [p["label"] for p in PROFILES])
    summary(results, t0, "Skin comparison complete")


if __name__ == "__main__":
    main()
