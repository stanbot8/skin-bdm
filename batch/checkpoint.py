"""Checkpoint-based batch execution for fork-style parallelism.

Usage:
    # Create a checkpoint at step 700 (day 7) from untreated baseline:
    python3 batch/checkpoint.py save --step 700

    # Run a simulation forking from a checkpoint with treatment overrides:
    python3 batch/checkpoint.py fork --checkpoint checkpoints/step_700 \
        --treatments npwt,msc

    # Full fork workflow: create checkpoint then run N configs from it:
    python3 batch/checkpoint.py workflow --step 700 --treatments npwt msc combination

Saves ~24% of sim time for day-7 forks, ~48% for day-14 forks.
"""

import argparse
import os
import shutil
import sys
import time

ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), os.pardir))
sys.path.insert(0, ROOT)

from batch.lib import (
    merge_config, apply_profile, apply_study, override_param,
    get_metrics_path, build_if_needed,
)


def _run_binary(env_extra=None):
    """Run the skibidy binary with optional extra env vars."""
    import subprocess
    env = os.environ.copy()
    if env_extra:
        env.update(env_extra)
    t0 = time.time()
    subprocess.run(
        [os.path.join(ROOT, "build", "skibidy")],
        cwd=ROOT, capture_output=True, text=True, env=env)
    elapsed = time.time() - t0
    return os.path.isfile(get_metrics_path()), elapsed


def prepare_config(profile="diabetic", study="diabetic-wound"):
    """Prepare a fresh bdm.toml with profile and study."""
    merge_config()
    if profile:
        apply_profile(profile)
    if study:
        apply_study(study)
    override_param("visualization.export", False)
    override_param("skin.headless", True)


def save_checkpoint(step, ckpt_dir, profile="diabetic", study="diabetic-wound"):
    """Run untreated baseline to given step and save checkpoint."""
    prepare_config(profile, study)

    ckpt_path = os.path.abspath(ckpt_dir)
    os.makedirs(ckpt_path, exist_ok=True)

    print(f"  Saving checkpoint at step {step} to {ckpt_path}")
    success, elapsed = _run_binary({
        "SKIBIDY_CKPT_SAVE_DIR": ckpt_path,
        "SKIBIDY_CKPT_STEP": str(step),
    })
    print(f"  {'OK' if success else 'FAILED'} ({elapsed:.0f}s)")
    return success


def fork_from_checkpoint(ckpt_dir, treatments=None, overrides=None,
                         profile="diabetic", study="diabetic-wound"):
    """Run simulation forking from a checkpoint with new params.

    Returns (success, elapsed).
    """
    prepare_config(profile, study)

    # Apply treatments
    if treatments:
        for tname in treatments:
            tpath = os.path.join(ROOT, "treatments", f"{tname}.toml")
            if os.path.isfile(tpath):
                from batch.lib import apply_overlay
                apply_overlay(tpath)

    # Apply overrides
    if overrides:
        for param_path, value in overrides.items():
            override_param(param_path, value)

    override_param("visualization.export", False)
    override_param("skin.headless", True)

    # Clean output
    output_dir = os.path.join(ROOT, "output")
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir, ignore_errors=True)

    return _run_binary({
        "SKIBIDY_CKPT_LOAD_DIR": os.path.abspath(ckpt_dir),
    })


def main():
    parser = argparse.ArgumentParser(description="Checkpoint-based batch execution")
    sub = parser.add_subparsers(dest="command")

    # Save subcommand
    sp = sub.add_parser("save", help="Save a checkpoint at a given step")
    sp.add_argument("--step", type=int, required=True,
                    help="Step number to checkpoint at")
    sp.add_argument("--dir", type=str, default=None,
                    help="Checkpoint directory (default: checkpoints/step_N)")
    sp.add_argument("--profile", type=str, default="diabetic")
    sp.add_argument("--study", type=str, default="diabetic-wound")

    # Fork subcommand
    fp = sub.add_parser("fork", help="Run from a checkpoint with treatments")
    fp.add_argument("--checkpoint", type=str, required=True,
                    help="Checkpoint directory")
    fp.add_argument("--treatments", type=str, nargs="+", default=[],
                    help="Treatment names to apply")
    fp.add_argument("--profile", type=str, default="diabetic")
    fp.add_argument("--study", type=str, default="diabetic-wound")

    args = parser.parse_args()

    if args.command == "save":
        ckpt_dir = args.dir or os.path.join(ROOT, "checkpoints", f"step_{args.step}")
        build_if_needed()
        save_checkpoint(args.step, ckpt_dir, args.profile, args.study)

    elif args.command == "fork":
        build_if_needed()
        success, elapsed = fork_from_checkpoint(
            args.checkpoint, treatments=args.treatments,
            profile=args.profile, study=args.study)
        print(f"  Fork: {'OK' if success else 'FAILED'} ({elapsed:.0f}s)")

    else:
        parser.print_help()


if __name__ == "__main__":
    main()
