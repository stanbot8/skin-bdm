#!/usr/bin/env python3
"""Rheumatoid arthritis study.

Runs RA wound healing with biologic treatment comparisons, disease severity
sweeps, and cartilage protection analysis using consensus batches.

Usage:
    python3 studies/rheumatoid/study.py                # baseline + 5 treatments
    python3 studies/rheumatoid/study.py --quick         # 5-run consensus
    python3 studies/rheumatoid/study.py -n 20           # 20-run consensus
    python3 studies/rheumatoid/study.py --baseline-only  # untreated RA only
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

TREATMENTS = [
    {"label": "untreated",     "extra": []},
    {"label": "anti_tnf",      "extra": ["--treatment", "anti_tnf"]},
    {"label": "tocilizumab",   "extra": ["--treatment", "tocilizumab"]},
    {"label": "methotrexate",  "extra": ["--treatment", "methotrexate"]},
    {"label": "jak_inhibitor", "extra": ["--treatment", "jak_inhibitor"]},
    {"label": "ra_combo",      "extra": ["--treatment", "ra_combination"]},
]


def main():
    parser = argparse.ArgumentParser(description="Rheumatoid arthritis study")
    parser.add_argument("-n", "--num-runs", type=int, default=10)
    parser.add_argument("--quick", action="store_true", help="5-run consensus")
    parser.add_argument("--baseline-only", action="store_true",
                        help="Run untreated RA baseline only")
    args = parser.parse_args()

    n_runs = 5 if args.quick else args.num_runs
    configs = TREATMENTS[:1] if args.baseline_only else TREATMENTS

    check_bdm()
    build_if_needed()

    results = os.path.join(os.path.dirname(__file__), "results")
    os.makedirs(results, exist_ok=True)

    banner("Rheumatoid Arthritis Study", [
        f"{n_runs}-run consensus per config",
        f"{len(configs)} configurations",
    ])

    t0 = time.time()

    for i, cfg in enumerate(configs, 1):
        print(f"[{i}/{len(configs)}] {cfg['label']}")
        consensus_dir = run_batch(
            n_runs, "rheumatoid", "rheumatoid",
            validate=False, extra_args=cfg["extra"],
        )
        copy_results(consensus_dir, results, cfg["label"])
        print()

    plot_metrics(results)
    labels = [c["label"] for c in configs]
    if len(labels) >= 2:
        compare_profiles(results, labels)
    summary(results, t0, "RA study complete")


if __name__ == "__main__":
    main()
