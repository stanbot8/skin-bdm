#!/usr/bin/env python3
"""Overlay multiple simulation runs on shared axes.

Accepts 2+ metrics.csv files and generates one figure per column group.
Only plots columns that have non-zero data in at least one input file.

Usage:
    python3 scripts/analysis/compare_runs.py run1/metrics.csv run2/metrics.csv
    python3 scripts/analysis/compare_runs.py --labels "Normal" "Aged" run1.csv run2.csv

Output: output/plots/compare_*.png
"""

import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load_csv(path):
    with open(path) as f:
        reader = csv.DictReader(
            (row for row in f if not row.startswith("#"))
        )
        rows = list(reader)
    data = {}
    for key in rows[0]:
        try:
            data[key] = [float(r[key]) for r in rows]
        except ValueError:
            data[key] = [r[key] for r in rows]
    return data


def has_data(datasets, key):
    """Check if any dataset has non-zero data for a key."""
    for data in datasets:
        if key in data and max(data[key]) > 0:
            return True
    return False


COLORS = ["#D46664", "#4A90D9", "#7B9F35", "#B5534B", "#8E6FBF",
          "#E89679", "#3D7A6B", "#C4943B"]

GROUPS = {
    "wound_closure": {
        "title": "Wound Closure",
        "columns": [("wound_closure_pct", "Closure (%)")],
    },
    "cell_populations": {
        "title": "Cell Populations",
        "columns": [
            ("n_basal", "Basal"),
            ("n_spinous", "Spinous"),
            ("n_granular", "Granular"),
            ("n_cornified", "Cornified"),
        ],
    },
    "field_means": {
        "title": "Wound Field Recovery",
        "columns": [
            ("mean_o2_wound", "O2"),
            ("mean_ca_wound", "Ca2+"),
            ("mean_perfusion_wound", "Perfusion"),
        ],
    },
    "immune": {
        "title": "Immune Kinetics",
        "columns": [
            ("n_neutrophils", "Neutrophils"),
            ("n_macrophages", "Macrophages"),
        ],
    },
    "fibroblast": {
        "title": "Fibroblast / Collagen",
        "columns": [
            ("n_myofibroblasts", "Myofibroblasts"),
            ("mean_collagen_wound", "Collagen"),
        ],
    },
    "tumor": {
        "title": "Tumor Growth",
        "columns": [("n_tumor_cells", "Tumor cells")],
    },
}


def main():
    args = sys.argv[1:]
    labels = []
    paths = []

    # Parse --labels
    if "--labels" in args:
        idx = args.index("--labels")
        args.pop(idx)
        while idx < len(args) and not args[idx].endswith(".csv"):
            labels.append(args.pop(idx))

    paths = args
    if len(paths) < 2:
        print("Usage: python3 scripts/analysis/compare_runs.py [--labels L1 L2 ...] CSV1 CSV2 [CSV3 ...]")
        sys.exit(1)

    for p in paths:
        if not os.path.exists(p):
            print(f"Error: {p} not found.")
            sys.exit(1)

    datasets = [load_csv(p) for p in paths]

    if not labels:
        labels = [os.path.basename(os.path.dirname(p)) or os.path.basename(p)
                  for p in paths]

    out_dir = "output/plots"
    os.makedirs(out_dir, exist_ok=True)

    for group_name, group in GROUPS.items():
        # Check if any column in this group has data
        active_cols = [(col, label) for col, label in group["columns"]
                       if has_data(datasets, col)]
        if not active_cols:
            continue

        nrows = len(active_cols)
        fig, axes = plt.subplots(nrows, 1, figsize=(8, 3 * nrows), sharex=True)
        if nrows == 1:
            axes = [axes]

        for row, (col, ylabel) in enumerate(active_cols):
            ax = axes[row]
            for i, (data, lbl) in enumerate(zip(datasets, labels)):
                if col not in data:
                    continue
                days = [h / 24.0 for h in data["time_h"]]
                color = COLORS[i % len(COLORS)]
                ax.plot(days, data[col], color=color, linewidth=2,
                        label=lbl, alpha=0.85)
            ax.set_ylabel(ylabel)
            ax.legend(fontsize=8)
            ax.grid(True, alpha=0.3)

        axes[0].set_title(group["title"])
        axes[-1].set_xlabel("Time (days)")
        fig.tight_layout()
        out_path = os.path.join(out_dir, f"compare_{group_name}.png")
        fig.savefig(out_path, dpi=150)
        plt.close(fig)
        print(f"  {out_path}")


if __name__ == "__main__":
    main()
