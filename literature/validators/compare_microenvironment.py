#!/usr/bin/env python3
"""Standalone microenvironment validation (TGF-b, VEGF, fibronectin, MMP, pH).

Usage:
    python3 literature/validators/compare_microenvironment.py [path/to/metrics.csv]
"""

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from lib import (load_csv, plots_dir, validate_microenvironment,
                 validate_ph, plot_microenvironment_panels, plot_ph_panel,
                 print_summary)


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
    ph_r = None
    if "mean_ph_wound" in sim and max(sim["mean_ph_wound"]) > 0:
        ph_r = validate_ph(sim, sim_days)
    print_summary(microenv=r, ph=ph_r)

    # 3x2 layout: 4 microenvironment panels + pH panel + 1 empty
    fig, axes = plt.subplots(3, 2, figsize=(12, 12))
    plot_microenvironment_panels(r, sim_days, axes[:2])
    if ph_r:
        plot_ph_panel(ph_r, sim_days, axes[2, 0])
    else:
        axes[2, 0].set_visible(False)
    axes[2, 1].set_visible(False)
    fig.suptitle("Microenvironment Validation", fontsize=14, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.96])

    out_dir = plots_dir(sim_path)
    os.makedirs(out_dir, exist_ok=True)
    fig.savefig(os.path.join(out_dir, "microenvironment_validation.png"), dpi=150)
    plt.close(fig)
    print(f"  Saved {out_dir}/microenvironment_validation.png")


if __name__ == "__main__":
    main()
