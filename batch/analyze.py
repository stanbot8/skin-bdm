#!/usr/bin/env python3
"""Post-hoc analysis of sweep or batch results.

Reads a results directory (from sweep.py or batch.py) and generates
plots and summary statistics.

Usage:
    python3 batch/analyze.py batch/results/cytokine_rate_20260225_143000/
    python3 batch/analyze.py batch/results/collagen_vs_scar_20260225_150000/
    python3 batch/analyze.py batch/results/consensus_20260225_160000/ --validate
"""

import argparse
import csv
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir))

from batch import lib


def load_summary(result_dir):
    """Load summary.csv from a sweep result directory."""
    path = os.path.join(result_dir, "summary.csv")
    if not os.path.isfile(path):
        return None
    with open(path) as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    # Convert numeric fields
    for row in rows:
        for k, v in row.items():
            try:
                row[k] = float(v)
            except (ValueError, TypeError):
                pass
    return rows


def detect_sweep_params(rows):
    """Detect which columns are sweep parameters (not outcomes)."""
    if not rows:
        return []
    skip = {"n_runs"}
    params = []
    for k in rows[0]:
        if k.endswith("_std") or k in skip:
            continue
        # Check if it varies across rows
        vals = set()
        for r in rows:
            vals.add(r[k])
        if len(vals) > 1 and not k.startswith("mean_") and not k.startswith("n_"):
            # Heuristic: sweep params are things like skin.immune.cytokine_rate
            if "." in str(list(vals)[0]) or isinstance(list(vals)[0], float):
                params.append(k)
    return params


def plot_single_param(rows, param, outcomes, result_dir):
    """Plot outcome vs single parameter with error bars."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available, skipping plots.")
        return

    x = [r[param] for r in rows]
    param_short = param.split(".")[-1] if "." in param else param

    n_outcomes = len(outcomes)
    fig, axes = plt.subplots(1, n_outcomes, figsize=(5 * n_outcomes, 4),
                             squeeze=False)
    axes = axes[0]

    for i, outcome in enumerate(outcomes):
        ax = axes[i]
        y = [r.get(outcome, 0) for r in rows]
        yerr = [r.get(f"{outcome}_std", 0) for r in rows]
        ax.errorbar(x, y, yerr=yerr, marker="o", capsize=4,
                    color="#D46664", linewidth=2, markersize=6)
        ax.set_xlabel(param_short)
        ax.set_ylabel(outcome.replace("_", " "))
        ax.set_title(outcome)
        ax.grid(True, alpha=0.3)

    fig.suptitle(f"Sweep: {param_short}", fontsize=12, fontweight="bold")
    fig.tight_layout()

    plot_path = os.path.join(result_dir, "sweep.png")
    fig.savefig(plot_path, dpi=150)
    plt.close(fig)
    print(f"Plot: {plot_path}")


def plot_two_param(rows, params, outcome, result_dir):
    """Plot heatmap for two-parameter sweep."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        print("matplotlib/numpy not available, skipping heatmap.")
        return

    p1, p2 = params[0], params[1]
    p1_short = p1.split(".")[-1]
    p2_short = p2.split(".")[-1]

    # Extract unique values
    x_vals = sorted(set(r[p1] for r in rows))
    y_vals = sorted(set(r[p2] for r in rows))

    # Build grid
    grid = np.full((len(y_vals), len(x_vals)), np.nan)
    for r in rows:
        xi = x_vals.index(r[p1])
        yi = y_vals.index(r[p2])
        grid[yi, xi] = r.get(outcome, float("nan"))

    fig, ax = plt.subplots(figsize=(6, 5))
    im = ax.imshow(grid, aspect="auto", origin="lower",
                   extent=[min(x_vals), max(x_vals), min(y_vals), max(y_vals)],
                   cmap="RdYlGn")
    ax.set_xlabel(p1_short)
    ax.set_ylabel(p2_short)
    ax.set_title(f"{outcome} ({p1_short} vs {p2_short})")
    fig.colorbar(im, ax=ax, label=outcome.replace("_", " "))

    # Annotate cells
    for r in rows:
        xi = r[p1]
        yi = r[p2]
        val = r.get(outcome, float("nan"))
        if not (val != val):  # not NaN
            ax.text(xi, yi, f"{val:.1f}", ha="center", va="center",
                    fontsize=7, color="black")

    fig.tight_layout()
    plot_path = os.path.join(result_dir, "heatmap.png")
    fig.savefig(plot_path, dpi=150)
    plt.close(fig)
    print(f"Heatmap: {plot_path}")


def print_summary_table(rows, params, outcomes):
    """Print a text summary table."""
    param_shorts = [p.split(".")[-1] for p in params]

    # Header
    header = param_shorts + ["n"] + outcomes
    widths = [max(len(h), 10) for h in header]

    print("\n" + "=" * sum(w + 2 for w in widths))
    print("  ".join(h.ljust(w) for h, w in zip(header, widths)))
    print("-" * sum(w + 2 for w in widths))

    for r in rows:
        parts = []
        for p in params:
            parts.append(f"{r[p]:.4g}".ljust(widths[len(parts)]))
        parts.append(f"{int(r.get('n_runs', 0))}".ljust(widths[len(parts)]))
        for o in outcomes:
            val = r.get(o, float("nan"))
            std = r.get(f"{o}_std", 0)
            text = f"{val:.3f}" if std == 0 else f"{val:.3f}+/-{std:.3f}"
            parts.append(text.ljust(widths[len(parts)]))
        print("  ".join(parts))

    print("=" * sum(w + 2 for w in widths))


def analyze_consensus(result_dir):
    """Analyze a batch consensus result (no sweep params)."""
    metrics_path = os.path.join(result_dir, "metrics.csv")
    if not os.path.isfile(metrics_path):
        print(f"No metrics.csv in {result_dir}")
        return

    data = lib.load_csv(metrics_path)
    print(f"\nConsensus metrics from: {result_dir}")
    print(f"  Steps: {len(data.get('step', []))}")
    print(f"  Final wound closure: {data.get('wound_closure_pct', [0])[-1]:.1f}%")

    # Count raw CSVs
    raw_dir = os.path.join(result_dir, "raw")
    if os.path.isdir(raw_dir):
        n_raw = len([f for f in os.listdir(raw_dir) if f.endswith(".csv")])
        print(f"  Runs: {n_raw}")


def main():
    parser = argparse.ArgumentParser(description="Analyze sweep/batch results")
    parser.add_argument("result_dir", help="Path to results directory")
    parser.add_argument("--validate", action="store_true",
                        help="Run validation on consensus metrics")
    args = parser.parse_args()

    result_dir = args.result_dir

    # Check for summary.csv (sweep result)
    rows = load_summary(result_dir)
    if rows:
        params = detect_sweep_params(rows)
        # Detect outcome columns
        skip = set(params) | {"n_runs"}
        outcomes = [k for k in rows[0]
                    if not k.endswith("_std") and k not in skip
                    and isinstance(rows[0][k], float)]

        print(f"Sweep results: {len(rows)} points, {len(params)} params")
        print(f"  Params: {params}")
        print(f"  Outcomes: {outcomes}")

        print_summary_table(rows, params, outcomes)

        if len(params) == 1:
            plot_single_param(rows, params[0], outcomes, result_dir)
        elif len(params) == 2 and outcomes:
            plot_two_param(rows, params, outcomes[0], result_dir)
            if len(outcomes) > 1:
                plot_single_param(rows, params[0], outcomes, result_dir)
    else:
        analyze_consensus(result_dir)

    if args.validate:
        metrics_path = os.path.join(result_dir, "metrics.csv")
        if os.path.isfile(metrics_path):
            lib.run_validation(metrics_path)
        else:
            print("No metrics.csv for validation.")


if __name__ == "__main__":
    main()
