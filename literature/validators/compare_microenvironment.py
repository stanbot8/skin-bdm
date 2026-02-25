#!/usr/bin/env python3
"""Standalone microenvironment validation (TGF-b, VEGF, fibronectin, MMP).

Usage:
    python3 literature/validators/compare_microenvironment.py [path/to/metrics.csv]
"""

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from lib import (load_csv, validate_microenvironment,
                 plot_microenvironment_panels, print_summary)


def main():
    sim_path = sys.argv[1] if len(sys.argv) > 1 else "output/skibidy/metrics.csv"

    if not os.path.exists(sim_path):
        print(f"Error: {sim_path} not found.")
        sys.exit(1)

    sim = load_csv(sim_path)
    sim_days = [h / 24.0 for h in sim["time_h"]]

    if "mean_tgfb_wound" not in sim or max(sim["mean_tgfb_wound"]) == 0:
        print("No microenvironment data found.")
        sys.exit(1)

    r = validate_microenvironment(sim, sim_days)
    print_summary(microenv=r)

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    plot_microenvironment_panels(r, sim_days, axes)
    fig.suptitle("Microenvironment Validation", fontsize=14, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.96])

    out_dir = "output/plots"
    os.makedirs(out_dir, exist_ok=True)
    fig.savefig(os.path.join(out_dir, "microenvironment_validation.png"), dpi=150)
    plt.close(fig)
    print(f"  Saved {out_dir}/microenvironment_validation.png")


if __name__ == "__main__":
    main()
