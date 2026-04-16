#!/usr/bin/env python3
"""Multi-run consensus: same config, multiple seeds.

Runs N simulations, computes mean/std across runs, writes
validation-compatible consensus metrics, and optionally validates.

Usage:
    python3 batch/batch.py                         # 20 runs, default config
    python3 batch/batch.py -n 5                    # 5 runs
    python3 batch/batch.py -n 10 --study wound     # 10 runs, wound study
    python3 batch/batch.py --skin aged --study wound --validate
    python3 batch/batch.py --skin diabetic --study diabetic-wound --validate
"""

import argparse
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir))

from batch import lib


def main():
    parser = argparse.ArgumentParser(
        description="Multi-run consensus (same config, different seeds)")
    parser.add_argument("-n", "--num-runs", type=int, default=20,
                        help="Number of runs (default: 20)")
    parser.add_argument("--study", type=str, default=None)
    parser.add_argument("--skin", type=str, default=None)
    parser.add_argument("--site", type=str, default=None,
                        help="Body site (e.g. foot_plantar, scalp)")
    parser.add_argument("--treatment", type=str, default=None,
                        help="Treatment overlay (e.g. anti_tnf, tocilizumab)")
    parser.add_argument("--validate", action="store_true",
                        help="Run validation on consensus")
    parser.add_argument("--no-validate", action="store_true",
                        help="Skip validation (default: skip)")
    parser.add_argument("--seed", type=int, default=None,
                        help="Base random seed (run i uses seed+i for reproducibility)")
    args = parser.parse_args()
    do_validate = args.validate and not args.no_validate

    study_label = args.study or "default"
    skin_label = args.skin or "default"
    study_dir = os.path.join(os.path.dirname(__file__), os.pardir,
                             "studies", study_label, "results")
    result_dir = study_dir
    raw_dir = os.path.join(result_dir, "raw")
    os.makedirs(raw_dir, exist_ok=True)

    print(f"=== Batch: {args.num_runs} runs ===")
    print(f"  Skin: {skin_label}, Study: {study_label}")
    if args.seed is not None:
        print(f"  Base seed: {args.seed}")
    print(f"  Output: {result_dir}/")
    print()

    # Setup config once
    lib.setup_run(args.skin, args.site, args.study, args.treatment)
    lib.build_if_needed()

    # Run simulations
    csv_paths = []
    t_start = time.time()
    for i in range(args.num_runs):
        lib.setup_run(args.skin, args.site, args.study, args.treatment)

        # Set reproducible seed: base_seed + run_index
        run_seed = (args.seed + i) if args.seed is not None else None
        if run_seed is not None:
            lib.override_param("simulation.random_seed", run_seed)

        label = f"[{i + 1}/{args.num_runs}]"
        print(f"  {label}", end=" ", flush=True)
        run_dir = os.path.join(raw_dir, f"run_{i:03d}")
        ok, elapsed = lib.run_simulation(output_path=run_dir)

        if not ok:
            print(f"FAILED ({elapsed:.0f}s)")
            continue

        csv_path = lib.get_metrics_path()
        if not os.path.isfile(csv_path):
            print(f"FAILED (no metrics)")
            continue

        csv_paths.append(csv_path)
        print(f"OK ({elapsed:.0f}s)")

    t_total = time.time() - t_start
    print(f"\n{len(csv_paths)}/{args.num_runs} completed "
          f"({t_total:.0f}s total, {t_total / max(1, args.num_runs):.0f}s avg)")

    # Write seed manifest for reproducibility
    if args.seed is not None:
        manifest_path = os.path.join(result_dir, "seed_manifest.csv")
        with open(manifest_path, "w", newline="") as f:
            import csv as _csv
            w = _csv.writer(f)
            w.writerow(["run", "seed", "status"])
            for i in range(args.num_runs):
                run_csv = os.path.join(raw_dir, f"run_{i:03d}",
                                       "skibidy", "metrics.csv")
                status = "ok" if os.path.isfile(run_csv) else "failed"
                w.writerow([i, args.seed + i, status])
        print(f"Seeds:     {manifest_path}")

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
