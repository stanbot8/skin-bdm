"""Shared utilities for study runner scripts.

Provides environment checks, build orchestration, banner printing,
timing, and result summary for all study scripts.
"""

import glob
import os
import subprocess
import sys
import time

REPO = os.path.normpath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
sys.path.insert(0, REPO)

from batch.lib import build_if_needed


def check_bdm():
    """Verify BioDynaMo environment is sourced."""
    if not os.environ.get("BDMSYS"):
        print("ERROR: BioDynaMo not sourced. Run: source <path>/bin/thisbdm.sh")
        sys.exit(1)


def banner(title, lines=None):
    """Print a study banner with optional detail lines."""
    print()
    print("=" * 44)
    print(f"  {title}")
    for line in (lines or []):
        print(f"  {line}")
    print("=" * 44)
    print()


def run_batch(n_runs, skin, study, validate=False, extra_args=None):
    """Run batch.py and return the study results directory path."""
    cmd = [
        sys.executable, "-u", os.path.join(REPO, "batch", "batch.py"),
        "-n", str(n_runs),
        "--skin", skin,
        "--study", study,
    ]
    if validate:
        cmd.append("--validate")
    if extra_args:
        cmd.extend(extra_args)
    subprocess.run(cmd, cwd=REPO, check=True)

    # batch.py writes to studies/{study}/results/
    return os.path.join(REPO, "studies", study, "results")


def copy_results(consensus_dir, results_dir, label):
    """Copy consensus and metrics CSVs from a batch run to results dir."""
    if not consensus_dir:
        return
    os.makedirs(results_dir, exist_ok=True)
    for src_name, dst_name in [
        ("consensus.csv", f"consensus_{label}.csv"),
        ("metrics.csv", f"metrics_{label}.csv"),
    ]:
        src = os.path.join(consensus_dir, src_name)
        if os.path.isfile(src):
            import shutil
            shutil.copy2(src, os.path.join(results_dir, dst_name))


def plot_metrics(results_dir, prefix=None):
    """Generate plots for a metrics CSV."""
    plot_script = os.path.join(REPO, "scripts", "analysis", "plot_metrics.py")
    csv_path = os.path.join(results_dir, f"metrics_{prefix}.csv") if prefix else None

    if prefix and csv_path and os.path.isfile(csv_path):
        plot_dir = os.path.join(results_dir, "plots")
        os.makedirs(plot_dir, exist_ok=True)
        print(f"Plotting: {prefix}")
        subprocess.run(
            [sys.executable, plot_script, csv_path, "-o", plot_dir,
             "--prefix", f"{prefix}_"],
            cwd=REPO, check=False)
    else:
        # Plot all metrics_*.csv files
        for csv_file in glob.glob(os.path.join(results_dir, "metrics_*.csv")):
            name = os.path.basename(csv_file).replace("metrics_", "").replace(".csv", "")
            plot_dir = os.path.join(results_dir, "plots")
            os.makedirs(plot_dir, exist_ok=True)
            print(f"Plotting: {name}")
            subprocess.run(
                [sys.executable, plot_script, csv_file, "-o", plot_dir,
                 "--prefix", f"{name}_"],
                cwd=REPO, check=False)


def compare_profiles(results_dir, labels):
    """Run cross-profile comparison on available metrics CSVs."""
    compare_script = os.path.join(REPO, "scripts", "analysis", "compare_sims.py")
    csvs = []
    found_labels = []
    for label in labels:
        f = os.path.join(results_dir, f"metrics_{label}.csv")
        if os.path.isfile(f):
            csvs.append(f)
            found_labels.append(label)
    if len(csvs) >= 2:
        print()
        print("Comparing profiles...")
        subprocess.run(
            [sys.executable, compare_script] + csvs +
            ["--labels", ",".join(found_labels)],
            cwd=REPO, check=False)


def summary(results_dir, t0, title="Study complete"):
    """Print result summary with timing."""
    elapsed = time.time() - t0
    mins = int(elapsed) // 60
    secs = int(elapsed) % 60

    n_csv = len(glob.glob(os.path.join(results_dir, "*.csv")))
    n_plots = len(glob.glob(os.path.join(results_dir, "plots", "*.png")))
    n_xlsx = len(glob.glob(os.path.join(results_dir, "*.xlsx")))

    print()
    print("=" * 44)
    print(f"  {title} ({mins}m {secs}s)")
    print("=" * 44)
    print(f"  Results in: {os.path.relpath(results_dir, REPO)}/")
    print(f"  CSVs:  {n_csv}")
    if n_xlsx:
        print(f"  Excel: {n_xlsx}")
    if n_plots:
        print(f"  Plots: {n_plots}")
    print()
    for f in sorted(os.listdir(results_dir)):
        print(f"    {f}")
    print()
