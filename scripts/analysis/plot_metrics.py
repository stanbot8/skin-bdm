#!/usr/bin/env python3
"""Generate publication-quality figures from metrics.csv.

Adaptive: only plots modules that have non-zero data.

Always produced:
  - cell_populations.png -- stacked area of cell types over time

Wound module (wound_closure_pct > 0):
  - wound_closure.png  -- wound closure % over time
  - field_recovery.png -- O2 and calcium recovery in wound cylinder

Immune module (n_neutrophils > 0):
  - immune_kinetics.png -- neutrophil and macrophage counts

Fibroblast module (n_myofibroblasts > 0):
  - fibroblast_kinetics.png -- myofibroblast count + collagen mean

Tumor module (n_tumor_cells > 0):
  - tumor_growth.png -- tumor cell count (log scale)

Usage:
    python3 scripts/analysis/plot_metrics.py [path/to/metrics.csv]

Defaults to output/skibidy/metrics.csv if no argument given.
"""

import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")  # headless rendering
import matplotlib.pyplot as plt


def load_csv(path):
    with open(path) as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    data = {}
    for key in rows[0]:
        try:
            data[key] = [float(r[key]) for r in rows]
        except ValueError:
            data[key] = [r[key] for r in rows]
    return data


def plot_wound_closure(data, out_dir):
    fig, ax = plt.subplots(figsize=(8, 4))
    days = [h / 24.0 for h in data["time_h"]]
    ax.plot(days, data["wound_closure_pct"], color="#D46664", linewidth=2)
    ax.axhline(90, color="gray", linestyle="--", linewidth=0.8, label="90% threshold")
    ax.set_xlabel("Time (days)")
    ax.set_ylabel("Wound closure (%)")
    ax.set_title("Wound Closure Over Time")
    ax.set_ylim(-2, 100)
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "wound_closure.png"), dpi=150)
    plt.close(fig)
    print(f"  wound_closure.png")


def plot_cell_populations(data, out_dir):
    fig, ax = plt.subplots(figsize=(8, 4))
    days = [h / 24.0 for h in data["time_h"]]

    colors = {
        "n_basal": "#D46664",
        "n_spinous": "#E89679",
        "n_granular": "#F0C987",
        "n_cornified": "#F5E6D3",
    }
    labels = {
        "n_basal": "Basal",
        "n_spinous": "Spinous",
        "n_granular": "Granular",
        "n_cornified": "Cornified",
    }
    keys = ["n_basal", "n_spinous", "n_granular", "n_cornified"]
    stacks = [data[k] for k in keys]

    ax.stackplot(
        days, *stacks,
        labels=[labels[k] for k in keys],
        colors=[colors[k] for k in keys],
        alpha=0.85,
    )
    ax.set_xlabel("Time (days)")
    ax.set_ylabel("Cell count")
    ax.set_title("Cell Populations by Stratum")
    ax.legend(loc="upper left")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "cell_populations.png"), dpi=150)
    plt.close(fig)
    print(f"  cell_populations.png")


def plot_field_recovery(data, out_dir):
    fig, ax1 = plt.subplots(figsize=(8, 4))
    days = [h / 24.0 for h in data["time_h"]]

    color_o2 = "#4A90D9"
    color_ca = "#B5534B"

    ax1.plot(days, data["mean_o2_wound"], color=color_o2, linewidth=2, label="O2")
    ax1.set_xlabel("Time (days)")
    ax1.set_ylabel("Mean O2 (wound)", color=color_o2)
    ax1.tick_params(axis="y", labelcolor=color_o2)

    ax2 = ax1.twinx()
    ax2.plot(days, data["mean_ca_wound"], color=color_ca, linewidth=2, label="Ca2+")
    ax2.set_ylabel("Mean Ca2+ (wound)", color=color_ca)
    ax2.tick_params(axis="y", labelcolor=color_ca)

    ax1.set_title("Wound Field Recovery: O2 and Calcium")
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc="center right")
    ax1.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "field_recovery.png"), dpi=150)
    plt.close(fig)
    print(f"  field_recovery.png")


def plot_immune_kinetics(data, out_dir):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    days = [h / 24.0 for h in data["time_h"]]

    ax1.plot(days, data["n_neutrophils"], color="#D46664", linewidth=2)
    ax1.set_ylabel("Neutrophils")
    ax1.set_title("Immune Cell Kinetics")
    ax1.grid(True, alpha=0.3)

    ax2.plot(days, data["n_macrophages"], color="#4A90D9", linewidth=2)
    ax2.set_xlabel("Time (days)")
    ax2.set_ylabel("Macrophages")
    ax2.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "immune_kinetics.png"), dpi=150)
    plt.close(fig)
    print(f"  immune_kinetics.png")


def plot_fibroblast_kinetics(data, out_dir):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    days = [h / 24.0 for h in data["time_h"]]

    ax1.plot(days, data["n_myofibroblasts"], color="#D46664", linewidth=2)
    ax1.set_ylabel("Myofibroblasts")
    ax1.set_title("Fibroblast / Collagen Kinetics")
    ax1.grid(True, alpha=0.3)

    ax2.plot(days, data["mean_collagen_wound"], color="#7B9F35", linewidth=2)
    ax2.set_xlabel("Time (days)")
    ax2.set_ylabel("Mean collagen (wound)")
    ax2.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "fibroblast_kinetics.png"), dpi=150)
    plt.close(fig)
    print(f"  fibroblast_kinetics.png")


def plot_ecm_dynamics(data, out_dir):
    fig, axes = plt.subplots(3, 1, figsize=(8, 8), sharex=True)
    days = [h / 24.0 for h in data["time_h"]]

    axes[0].plot(days, data["mean_mmp_wound"], color="#D46664", linewidth=2)
    axes[0].set_ylabel("Mean MMP (wound)")
    axes[0].set_title("ECM Remodeling Dynamics")
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(days, data["mean_fibronectin_wound"], color="#4A90D9", linewidth=2)
    axes[1].set_ylabel("Mean fibronectin (wound)")
    axes[1].grid(True, alpha=0.3)

    axes[2].plot(days, data["mean_collagen_wound"], color="#7B9F35", linewidth=2)
    axes[2].set_xlabel("Time (days)")
    axes[2].set_ylabel("Mean collagen (wound)")
    axes[2].grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "ecm_dynamics.png"), dpi=150)
    plt.close(fig)
    print(f"  ecm_dynamics.png")


def plot_tumor_growth(data, out_dir):
    fig, ax = plt.subplots(figsize=(8, 4))
    days = [h / 24.0 for h in data["time_h"]]

    sim_agents = data["n_tumor_cells"]
    if "tumor_field_cells" in data:
        sim_total = [a + f for a, f in zip(sim_agents, data["tumor_field_cells"])]
    else:
        sim_total = sim_agents

    ax.plot(days, sim_total, color="#D46664", linewidth=2, label="Total (agents + field)")
    ax.plot(days, sim_agents, color="#4A90D9", linewidth=1.5,
            linestyle="--", label="Active agents")
    ax.set_xlabel("Time (days)")
    ax.set_ylabel("Tumor cells")
    ax.set_title("Tumor Growth")
    ax.set_yscale("log")
    ax.set_ylim(bottom=1)
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "tumor_growth.png"), dpi=150)
    plt.close(fig)
    print(f"  tumor_growth.png")


def main():
    csv_path = "output/skibidy/metrics.csv"
    out_dir = None
    args = sys.argv[1:]
    i = 0
    while i < len(args):
        if args[i] in ("-o", "--output") and i + 1 < len(args):
            out_dir = args[i + 1]
            i += 2
        else:
            csv_path = args[i]
            i += 1

    if not os.path.exists(csv_path):
        print(f"Error: {csv_path} not found. Run the simulation first.")
        sys.exit(1)

    if out_dir is None:
        out_dir = "output/plots"
    os.makedirs(out_dir, exist_ok=True)
    data = load_csv(csv_path)

    # Detect active modules
    has_wound = ("wound_closure_pct" in data and max(data["wound_closure_pct"]) > 0)
    has_immune = ("n_neutrophils" in data and max(data["n_neutrophils"]) > 0)
    has_fibroblast = ("n_myofibroblasts" in data and max(data["n_myofibroblasts"]) > 0)
    has_tumor = ("n_tumor_cells" in data and max(data.get("n_tumor_cells", [0])) > 0)
    has_ecm = ("mean_mmp_wound" in data and max(data.get("mean_mmp_wound", [0])) > 0)

    print(f"Plotting {csv_path} -> {out_dir}/")
    plot_cell_populations(data, out_dir)
    if has_wound:
        plot_wound_closure(data, out_dir)
        plot_field_recovery(data, out_dir)
    if has_immune:
        plot_immune_kinetics(data, out_dir)
    if has_fibroblast:
        plot_fibroblast_kinetics(data, out_dir)
    if has_ecm:
        plot_ecm_dynamics(data, out_dir)
    if has_tumor:
        plot_tumor_growth(data, out_dir)
    print("Done.")


if __name__ == "__main__":
    main()
