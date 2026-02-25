#!/usr/bin/env python3
"""Standalone fibroblast/collagen validation.

Usage:
    python3 literature/validators/compare_fibroblast.py [path/to/metrics.csv]
"""

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from lib import load_csv, validate_fibroblast, plot_fibroblast_panels, print_summary


def main():
    sim_path = sys.argv[1] if len(sys.argv) > 1 else "output/skibidy/metrics.csv"

    if not os.path.exists(sim_path):
        print(f"Error: {sim_path} not found.")
        sys.exit(1)

    sim = load_csv(sim_path)
    sim_days = [h / 24.0 for h in sim["time_h"]]

    if "n_myofibroblasts" not in sim or max(sim["n_myofibroblasts"]) == 0:
        print("No fibroblast data found. Set fibroblast_enabled = true and re-run.")
        sys.exit(1)

    r = validate_fibroblast(sim, sim_days)
    print_summary(fibroblast=r)

    fig, axes = plt.subplots(3, 1, figsize=(8, 9), sharex=True)
    plot_fibroblast_panels(r, sim_days, axes)
    fig.suptitle("Fibroblast/Collagen Validation", fontsize=14, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.96])

    out_dir = "output/plots"
    os.makedirs(out_dir, exist_ok=True)
    fig.savefig(os.path.join(out_dir, "fibroblast_validation.png"), dpi=150)
    plt.close(fig)
    print(f"  Saved {out_dir}/fibroblast_validation.png")


if __name__ == "__main__":
    main()
