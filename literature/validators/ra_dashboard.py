#!/usr/bin/env python3
"""RA analysis dashboard: multi-panel overview of all RA observables.

Generates a 3x2 grid showing TNF-alpha, IL-6, cartilage, bone, T cell,
and synovial pannus dynamics from one or more metrics.csv files.
When multiple files are provided, overlays them for treatment comparison.

Usage:
    # Single run analysis
    python3 literature/validators/ra_dashboard.py path/to/metrics.csv

    # Treatment comparison (overlay multiple runs)
    python3 literature/validators/ra_dashboard.py \\
        --label "Untreated" results/baseline/metrics.csv \\
        --label "Anti-TNF"  results/anti_tnf/metrics.csv \\
        --label "Tocilizumab" results/tocilizumab/metrics.csv

    # From experiment results directory
    python3 literature/validators/ra_dashboard.py --experiment results/experiments/biologic_comparison/
"""

import glob
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from lib import (load_csv, plots_dir, validate_ra, interpolate,
                 peak_normalize, compute_rmse,
                 SIM_COLOR, REF_COLOR, REF_KW,
                 _study_ref_path)


COLORS = [
    "#D46664",  # red/salmon (default sim)
    "#4A90D9",  # blue
    "#7B9F35",  # olive
    "#9B59B6",  # purple
    "#E67E22",  # orange
    "#1ABC9C",  # teal
    "#E74C3C",  # crimson
    "#3498DB",  # sky blue
]


def load_run(csv_path):
    """Load a metrics CSV and compute RA validation."""
    sim = load_csv(csv_path)
    sim_days = [h / 24.0 for h in sim["time_h"]]
    r = validate_ra(sim, sim_days)
    return sim, sim_days, r


def load_experiment_dir(exp_dir):
    """Load all config results from an experiment output directory."""
    runs = []
    for config_dir in sorted(glob.glob(os.path.join(exp_dir, "*"))):
        if not os.path.isdir(config_dir):
            continue
        csv_path = os.path.join(config_dir, "metrics.csv")
        if not os.path.exists(csv_path):
            csv_path = os.path.join(config_dir, "consensus_metrics.csv")
        if not os.path.exists(csv_path):
            continue
        label = os.path.basename(config_dir).replace("_", " ").title()
        sim, sim_days, r = load_run(csv_path)
        runs.append((label, sim, sim_days, r))
    return runs


def _plot_panel(ax, runs, key, ref_key, ref_col, ylabel, title,
                ylim=(-0.05, 1.15), rmse_key=None, rmse_pos=(0.98, 0.85),
                rmse_va="top", legend_loc="best"):
    """Plot a single observable panel with optional multi-run overlay."""
    for i, (label, sim, sim_days, r) in enumerate(runs):
        vals = r.get(key, [])
        if not vals:
            continue
        color = COLORS[i % len(COLORS)]
        lw = 2.0 if i == 0 else 1.5
        ax.plot(sim_days, vals, color=color, linewidth=lw, label=label)

    # Reference curve (from first run)
    _, _, _, r0 = runs[0]
    if ref_key in r0 and r0[ref_key]:
        ref = r0[ref_key]
        ax.plot(ref["day"], ref[ref_col], color="#888888", linewidth=1.5,
                linestyle="--", marker="o", markersize=3, label="Literature")

    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_ylim(ylim)
    ax.set_xlim(0, 30)
    ax.legend(fontsize=7, loc=legend_loc)
    ax.grid(True, alpha=0.3)

    # RMSE annotation (first run only)
    if rmse_key and rmse_key in r0:
        ax.text(rmse_pos[0], rmse_pos[1],
                f"RMSE = {r0[rmse_key] * 100:.1f}%",
                transform=ax.transAxes, ha="right", va=rmse_va,
                fontsize=8, color="gray")


def dashboard(runs, out_path):
    """Generate the 3x2 RA dashboard PNG."""
    fig, axes = plt.subplots(3, 2, figsize=(14, 10), sharex=True)

    # Row 0: cytokines
    _plot_panel(axes[0, 0], runs, "sim_tnf", "ref_tnf", "tnf_alpha_normalized",
                "TNF-alpha (normalized)", "TNF-alpha Kinetics",
                rmse_key="tnf_rmse")
    _plot_panel(axes[0, 1], runs, "sim_il6", "ref_il6", "il6_normalized",
                "IL-6 (normalized)", "IL-6 Kinetics",
                rmse_key="il6_rmse")

    # Row 1: structural integrity
    _plot_panel(axes[1, 0], runs, "sim_cart", "ref_cart", "cartilage_integrity",
                "Cartilage integrity", "Cartilage Erosion",
                rmse_key="cart_rmse", rmse_pos=(0.98, 0.15), rmse_va="bottom",
                legend_loc="lower left")
    _plot_panel(axes[1, 1], runs, "sim_bone", "ref_bone", "bone_integrity",
                "Bone integrity", "Subchondral Bone Erosion",
                rmse_key="bone_rmse", rmse_pos=(0.98, 0.15), rmse_va="bottom",
                legend_loc="lower left")

    # Row 2: cellular/tissue
    _plot_panel(axes[2, 0], runs, "sim_tcell", "ref_tcell", "tcell_normalized",
                "T cell density (normalized)", "T Cell Infiltration",
                rmse_key="tcell_rmse")
    _plot_panel(axes[2, 1], runs, "sim_syn", "ref_syn", "synovial_normalized",
                "Pannus density (normalized)", "Synovial Pannus Growth",
                rmse_key="syn_rmse")

    for ax in axes[-1]:
        ax.set_xlabel("Time (days)")

    n_runs = len(runs)
    title = "RA Dashboard" if n_runs == 1 else f"RA Treatment Comparison ({n_runs} configs)"
    fig.suptitle(title, fontsize=14, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"  Saved {out_path}")


def print_table(runs):
    """Print RMSE summary table."""
    header = f"{'Config':<30} {'TNF':>7} {'IL-6':>7} {'Cart':>7} {'Bone':>7} {'Tcell':>7} {'Pannus':>7}"
    print(header)
    print("-" * len(header))
    for label, sim, sim_days, r in runs:
        tnf = f"{r['tnf_rmse']*100:.1f}%" if r.get("tnf_rmse") else "  N/A"
        il6 = f"{r['il6_rmse']*100:.1f}%" if r.get("il6_rmse") else "  N/A"
        cart = f"{r['cart_rmse']*100:.1f}%" if r.get("cart_rmse") else "  N/A"
        bone = f"{r['bone_rmse']*100:.1f}%" if r.get("has_bone") else "  N/A"
        tcell = f"{r['tcell_rmse']*100:.1f}%" if r.get("has_tcell") else "  N/A"
        syn = f"{r['syn_rmse']*100:.1f}%" if r.get("has_syn") else "  N/A"
        print(f"{label:<30} {tnf:>7} {il6:>7} {cart:>7} {bone:>7} {tcell:>7} {syn:>7}")


def main():
    args = sys.argv[1:]
    runs = []
    exp_dir = None

    i = 0
    while i < len(args):
        if args[i] == "--label" and i + 2 < len(args):
            label = args[i + 1]
            csv_path = args[i + 2]
            if not os.path.exists(csv_path):
                print(f"Error: {csv_path} not found.")
                sys.exit(1)
            sim, sim_days, r = load_run(csv_path)
            runs.append((label, sim, sim_days, r))
            i += 3
        elif args[i] == "--experiment" and i + 1 < len(args):
            exp_dir = args[i + 1]
            i += 2
        elif not args[i].startswith("--"):
            csv_path = args[i]
            if not os.path.exists(csv_path):
                print(f"Error: {csv_path} not found.")
                sys.exit(1)
            sim, sim_days, r = load_run(csv_path)
            label = os.path.basename(os.path.dirname(csv_path))
            runs.append((label, sim, sim_days, r))
            i += 1
        else:
            i += 1

    if exp_dir:
        if not os.path.isdir(exp_dir):
            print(f"Error: {exp_dir} not found.")
            sys.exit(1)
        runs.extend(load_experiment_dir(exp_dir))

    if not runs:
        # Default: try output/skibidy/metrics.csv
        default = "output/skibidy/metrics.csv"
        if os.path.exists(default):
            sim, sim_days, r = load_run(default)
            runs.append(("Simulation", sim, sim_days, r))
        else:
            print("Usage: ra_dashboard.py [--label NAME] path/to/metrics.csv ...")
            print("       ra_dashboard.py --experiment path/to/experiment_results/")
            sys.exit(1)

    # Check first run has RA data
    _, sim0, _, r0 = runs[0]
    if "mean_tnf_alpha_wound" not in sim0 or max(sim0["mean_tnf_alpha_wound"]) == 0:
        print("No RA data found. Set ra_enabled = true and re-run.")
        sys.exit(1)

    print_table(runs)

    # Output path: next to first CSV
    first_csv = None
    for a in sys.argv[1:]:
        if not a.startswith("--") and os.path.exists(a):
            first_csv = a
            break
    if exp_dir:
        out_dir = exp_dir
    elif first_csv:
        out_dir = plots_dir(first_csv)
    else:
        out_dir = plots_dir("output/skibidy/metrics.csv")
    os.makedirs(out_dir, exist_ok=True)
    dashboard(runs, os.path.join(out_dir, "ra_dashboard.png"))


if __name__ == "__main__":
    main()
