#!/usr/bin/env python3
"""Multi-run consensus: same config, multiple seeds.

Runs N simulations, computes mean/std across runs, writes
validation-compatible consensus metrics, and optionally validates.

Usage:
    python3 batch/batch.py                         # 20 runs, default config
    python3 batch/batch.py -n 5                    # 5 runs
    python3 batch/batch.py -n 10 --preset wound    # 10 runs, wound preset
    python3 batch/batch.py --skin aged --preset wound --validate
"""

import argparse
import os
import shutil
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir))

from datetime import datetime

from batch import lib


def main():
    parser = argparse.ArgumentParser(
        description="Multi-run consensus (same config, different seeds)")
    parser.add_argument("-n", "--num-runs", type=int, default=20,
                        help="Number of runs (default: 20)")
    parser.add_argument("--preset", type=str, default=None)
    parser.add_argument("--skin", type=str, default=None)
    parser.add_argument("--validate", action="store_true",
                        help="Run validation on consensus")
    parser.add_argument("--no-validate", action="store_true",
                        help="Skip validation (default: skip)")
    args = parser.parse_args()
    do_validate = args.validate and not args.no_validate

    skin_label = args.skin or "default"
    preset_label = args.preset or "default"
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    result_dir = os.path.join(
        os.path.dirname(__file__), "results",
        f"consensus_{skin_label}_{preset_label}_{timestamp}")
    raw_dir = os.path.join(result_dir, "raw")
    os.makedirs(raw_dir, exist_ok=True)

    print(f"=== Batch: {args.num_runs} runs ===")
    print(f"  Skin: {skin_label}, Preset: {preset_label}")
    print(f"  Output: {result_dir}/")
    print()

    # Setup config once
    lib.merge_config()
    if args.skin:
        lib.apply_profile(args.skin)
    if args.preset:
        lib.apply_preset(args.preset)
    lib.build_if_needed()

    # Run simulations
    csv_paths = []
    t_start = time.time()
    for i in range(args.num_runs):
        # Re-merge config each run (fresh bdm.toml)
        lib.merge_config()
        if args.skin:
            lib.apply_profile(args.skin)
        if args.preset:
            lib.apply_preset(args.preset)

        label = f"[{i + 1}/{args.num_runs}]"
        print(f"  {label}", end=" ", flush=True)
        ok, elapsed = lib.run_simulation()

        if not ok:
            print(f"FAILED ({elapsed:.0f}s)")
            continue

        src = lib.get_metrics_path()
        if not os.path.isfile(src):
            print(f"FAILED (no metrics)")
            continue

        dst = os.path.join(raw_dir, f"run_{i:03d}.csv")
        shutil.copy2(src, dst)
        csv_paths.append(dst)
        print(f"OK ({elapsed:.0f}s)")

    t_total = time.time() - t_start
    print(f"\n{len(csv_paths)}/{args.num_runs} completed "
          f"({t_total:.0f}s total, {t_total / max(1, args.num_runs):.0f}s avg)")

    if not csv_paths:
        print("No successful runs.")
        sys.exit(1)

    # Compute consensus
    mean_data, std_data = lib.aggregate_csvs(csv_paths)
    columns = list(mean_data.keys())

    # Write validation-compatible metrics.csv (means only)
    lib.write_csv(mean_data, os.path.join(result_dir, "metrics.csv"), columns)

    # Write full consensus with std
    full_path = os.path.join(result_dir, "consensus.csv")
    full_data = {}
    for col in columns:
        full_data[col] = mean_data[col]
        full_data[f"{col}_std"] = std_data[col]
    full_cols = []
    for col in columns:
        full_cols.append(col)
        full_cols.append(f"{col}_std")
    lib.write_csv(full_data, full_path, full_cols)

    print(f"Consensus: {full_path}")
    print(f"Metrics:   {os.path.join(result_dir, 'metrics.csv')}")

    # Validate
    if do_validate:
        print("\n=== Validation ===")
        lib.run_validation(os.path.join(result_dir, "metrics.csv"))

    print(f"\nResults: {result_dir}/")


if __name__ == "__main__":
    main()
