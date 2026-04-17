"""Run N simulation replicates with distinct seeds and save each CSV.

Used by the wiggum loop to get multi-replicate signal for mechanistic
validation. Runs sequentially; full-config sim takes ~4 min each in WSL.

Usage:
    python3 scripts/wiggum_replicate.py [--n N] [--skin S] [--study ST] [--out-dir DIR]
"""
import argparse
import os
import shutil
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from batch import lib


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--n", type=int, default=3)
    p.add_argument("--skin", default="normal")
    p.add_argument("--study", default="wound")
    p.add_argument("--seed-start", type=int, default=42)
    p.add_argument("--out-dir", default="/tmp/replicate_runs")
    p.add_argument("--build-dir", default="build")
    args = p.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    for i in range(args.n):
        seed = args.seed_start + i
        print(f"[replicate {i + 1}/{args.n}] seed={seed}", flush=True)
        lib.setup_run(skin=args.skin, study=args.study)
        lib.override_param("simulation.random_seed", seed)
        shutil.copy("bdm.toml", os.path.join(args.build_dir, "bdm.toml"))
        out = os.path.join(args.build_dir, "output")
        if os.path.exists(out):
            shutil.rmtree(out)
        r = subprocess.run(
            ["./skibidy"],
            cwd=args.build_dir,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            timeout=400,
        )
        if r.returncode != 0:
            print(f"  FAIL: skibidy exit={r.returncode}", flush=True)
            continue
        src = os.path.join(out, "skibidy", "metrics.csv")
        dst = os.path.join(args.out_dir, f"metrics_seed{seed}.csv")
        shutil.copy(src, dst)
        print(f"  saved {dst}", flush=True)

    print("DONE")


if __name__ == "__main__":
    main()
