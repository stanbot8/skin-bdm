#!/usr/bin/env python3
"""Standalone tumor growth validation.

Usage:
    python3 literature/validators/compare_tumor.py [path/to/metrics.csv]
"""

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from lib import (load_csv, validate_tumor, plot_tumor_panels, print_summary,
                 surface_fraction)


def main():
    sim_path = sys.argv[1] if len(sys.argv) > 1 else "output/skibidy/metrics.csv"

    if not os.path.exists(sim_path):
        print(f"Error: {sim_path} not found.")
        sys.exit(1)

    sim = load_csv(sim_path)
    sim_days = [h / 24.0 for h in sim["time_h"]]

    if "n_tumor_cells" not in sim:
        print("No tumor data found. Set tumor_enabled = true and re-run.")
        sys.exit(1)

    r = validate_tumor(sim, sim_days)
    print_summary(tumor=r)

    # Extended scale context (standalone only)
    n_init = r["n_init"]
    n_final = r["obs_final"]
    sf_init = surface_fraction(n_init)
    sf_final = surface_fraction(n_final)
    sf_established = surface_fraction(1e5)
    print(f"\n  Scale context (Gompertzian):")
    print(f"    Initial cells:              {n_init:.0f}  (surface fraction {sf_init:.0%})")
    print(f"    Final cells:                {n_final:.0f}  (surface fraction {sf_final:.0%})")
    print(f"    Established BCC (~10^5):    surface fraction {sf_established:.0%}")
    print(f"    Reference Ki-67 ({r['bcc_ki67_pct']:.1f}%) applies at established-tumor scale.")

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    plot_tumor_panels(r, sim_days, axes)
    fig.suptitle("Tumor Growth Validation", fontsize=14, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.92])

    out_dir = "output/plots"
    os.makedirs(out_dir, exist_ok=True)
    fig.savefig(os.path.join(out_dir, "tumor_validation.png"), dpi=150)
    plt.close(fig)
    print(f"  Saved {out_dir}/tumor_validation.png")


if __name__ == "__main__":
    main()
