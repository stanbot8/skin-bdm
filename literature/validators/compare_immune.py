#!/usr/bin/env python3
"""Standalone immune cell kinetics validation.

Usage:
    python3 literature/validators/compare_immune.py [path/to/metrics.csv]
"""

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from lib import (load_csv, detect_condition, validate_wound,
                 print_summary, SIM_KW, REF_KW)


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

    print(f"Immune cell kinetics validation ({len(sim_days)} sim points vs literature)")
    print(f"  Neutrophils:  RMSE = {r['neut_rmse'] * 100:.2f} %  (peak count = {r['neut_peak']:.0f})")
    print(f"  Macrophages:  RMSE = {r['mac_rmse'] * 100:.2f} %  (peak count = {r['mac_peak']:.0f})")

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)

    ax1.plot(sim_days, r["sim_neut"], **SIM_KW)
    ax1.plot(r["ref_immune"]["day"], r["ref_immune"]["neutrophils_normalized"], **REF_KW)
    ax1.set_ylabel("Neutrophils (normalized)")
    ax1.set_title("Immune Cell Kinetics: Simulation vs Literature")
    ax1.set_ylim(-0.05, 1.15)
    ax1.legend(loc="upper right")
    ax1.grid(True, alpha=0.3)
    ax1.text(0.98, 0.85, f"RMSE = {r['neut_rmse'] * 100:.1f}%",
             transform=ax1.transAxes, ha="right", va="top",
             fontsize=9, color="gray")

    ax2.plot(sim_days, r["sim_mac"], **SIM_KW)
    ax2.plot(r["ref_immune"]["day"], r["ref_immune"]["macrophages_normalized"], **REF_KW)
    ax2.set_xlabel("Time (days)")
    ax2.set_ylabel("Macrophages (normalized)")
    ax2.set_xlim(0, 30)
    ax2.set_ylim(-0.05, 1.15)
    ax2.legend(loc="upper right")
    ax2.grid(True, alpha=0.3)
    ax2.text(0.98, 0.85, f"RMSE = {r['mac_rmse'] * 100:.1f}%",
             transform=ax2.transAxes, ha="right", va="top",
             fontsize=9, color="gray")

    out_dir = "output/plots"
    os.makedirs(out_dir, exist_ok=True)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "immune_validation.png"), dpi=150)
    plt.close(fig)
    print(f"  Saved {out_dir}/immune_validation.png")


if __name__ == "__main__":
    main()
