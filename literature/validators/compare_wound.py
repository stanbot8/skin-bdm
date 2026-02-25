#!/usr/bin/env python3
"""Standalone wound healing validation.

Usage:
    python3 literature/validators/compare_wound.py [path/to/metrics.csv]
"""

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from lib import (load_csv, detect_condition, validate_wound,
                 plot_wound_panels, print_summary)


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    flags = [a for a in sys.argv[1:] if a.startswith("--")]
    sim_path = args[0] if args else "output/skibidy/metrics.csv"

    condition = None
    for f in flags:
        if f == "--diabetic":
            condition = "diabetic"
        elif f == "--normal":
            condition = "normal"
    if condition is None:
        condition = detect_condition()

    if not os.path.exists(sim_path):
        print(f"Error: {sim_path} not found.")
        sys.exit(1)

    sim = load_csv(sim_path)
    sim_days = [h / 24.0 for h in sim["time_h"]]

    if "wound_closure_pct" not in sim or max(sim["wound_closure_pct"]) == 0:
        print("No wound data found. Set wound_enabled = true and re-run.")
        sys.exit(1)

    r = validate_wound(sim, sim_days, condition)
    print_summary(wound=r)

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    plot_wound_panels(r, sim_days, axes)
    fig.suptitle("Wound Healing Validation", fontsize=14, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.96])

    out_dir = "output/plots"
    os.makedirs(out_dir, exist_ok=True)
    fig.savefig(os.path.join(out_dir, "wound_validation.png"), dpi=150)
    plt.close(fig)
    print(f"  Saved {out_dir}/wound_validation.png")


if __name__ == "__main__":
    main()
