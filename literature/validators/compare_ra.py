#!/usr/bin/env python3
"""Standalone RA cytokine and cartilage validation.

Usage:
    python3 literature/validators/compare_ra.py [path/to/metrics.csv]
"""

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from lib import load_csv, plots_dir, validate_ra, plot_ra_panels, SIM_KW, REF_KW


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    sim_path = args[0] if args else "output/skibidy/metrics.csv"

    if not os.path.exists(sim_path):
        print(f"Error: {sim_path} not found.")
        sys.exit(1)

    sim = load_csv(sim_path)
    sim_days = [h / 24.0 for h in sim["time_h"]]

    if "mean_tnf_alpha_wound" not in sim or max(sim["mean_tnf_alpha_wound"]) == 0:
        print("No RA data found. Set ra_enabled = true and re-run.")
        sys.exit(1)

    r = validate_ra(sim, sim_days)

    print(f"RA validation ({len(sim_days)} sim points vs literature)")
    print(f"  TNF-alpha:  RMSE = {r['tnf_rmse'] * 100:.2f} %  (peak = {r['tnf_peak']:.4f})")
    print(f"    Flare 0-7d:   RMSE = {r['tnf_flare_rmse'] * 100:.2f} %")
    print(f"    Chronic 7-30d: RMSE = {r['tnf_chronic_rmse'] * 100:.2f} %")
    print(f"  IL-6:       RMSE = {r['il6_rmse'] * 100:.2f} %  (peak = {r['il6_peak']:.4f})")
    print(f"  Cartilage:  RMSE = {r['cart_rmse'] * 100:.2f} %")
    if r.get("has_bone"):
        print(f"  Bone:       RMSE = {r['bone_rmse'] * 100:.2f} %")
    if r.get("has_tcell"):
        print(f"  T cells:    RMSE = {r['tcell_rmse'] * 100:.2f} %  (peak = {r['tcell_peak']:.4f})")
    if r.get("has_syn"):
        print(f"  Pannus:     RMSE = {r['syn_rmse'] * 100:.2f} %  (peak = {r['syn_peak']:.4f})")

    has_ext = r.get("has_bone") or r.get("has_tcell") or r.get("has_syn")
    if has_ext:
        fig, ax_arr = plt.subplots(3, 2, figsize=(12, 9), sharex=True)
        axes = [ax_arr[0, 0], ax_arr[0, 1], ax_arr[1, 0],
                ax_arr[1, 1], ax_arr[2, 0], ax_arr[2, 1]]
    else:
        fig, axes = plt.subplots(3, 1, figsize=(8, 9), sharex=True)
    plot_ra_panels(r, sim_days, axes)
    fig.suptitle("Rheumatoid Arthritis Validation", fontsize=14, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.96])

    out_dir = plots_dir(sim_path)
    os.makedirs(out_dir, exist_ok=True)
    fig.savefig(os.path.join(out_dir, "ra_validation.png"), dpi=150)
    plt.close(fig)
    print(f"  Saved {out_dir}/ra_validation.png")


if __name__ == "__main__":
    main()
