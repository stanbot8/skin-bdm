#!/usr/bin/env python3
"""Validation dashboard for skibidy.

Single entry point: loads metrics once, computes once, prints once, plots once.
Generates validation_dashboard.png plus per-module PNGs.

Usage:
    python3 literature/validate_all.py [path/to/metrics.csv] [--normal|--diabetic|--burn|--pressure|--surgical|--rheumatoid]

Condition auto-detection: reads bdm.toml for study-specific sections.
Override with explicit flags.
"""

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(__file__))
from lib import (load_csv, plots_dir, detect_condition, detect_modules,
                 validate_wound, validate_fibroblast, validate_tumor,
                 validate_microenvironment, validate_ph, validate_ra,
                 plot_wound_panels, plot_fibroblast_panels, plot_tumor_panels,
                 plot_microenvironment_panels, plot_ph_panel, plot_ra_panels,
                 print_summary)
from check_sources import run_checks as check_sources


def main():
    # --- Source integrity check (no sim data needed) ---
    src_errors, src_warnings = check_sources()
    for w in src_warnings:
        print(f"  WARN: {w}")
    if src_errors:
        for e in src_errors:
            print(f"  ERROR: {e}")
        print(f"  SOURCES check FAILED: {len(src_errors)} error(s)")
        sys.exit(1)
    n = len(src_warnings)
    print(f"  SOURCES check passed ({n} warning{'s' if n != 1 else ''})")

    # --- Parse CLI args ---
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    flags = [a for a in sys.argv[1:] if a.startswith("--")]
    sim_path = args[0] if args else "output/skibidy/metrics.csv"

    # Condition: CLI flag > auto-detect from bdm.toml
    condition = None
    for f in flags:
        if f == "--diabetic":
            condition = "diabetic"
        elif f == "--burn":
            condition = "burn"
        elif f == "--pressure":
            condition = "pressure"
        elif f == "--surgical":
            condition = "surgical"
        elif f == "--rheumatoid":
            condition = "rheumatoid"
        elif f == "--normal":
            condition = "normal"
    if condition is None:
        condition = detect_condition()
    print(f"  Condition: {condition}")

    if not os.path.exists(sim_path):
        print(f"Error: {sim_path} not found. Run the simulation first.")
        sys.exit(1)

    sim = load_csv(sim_path)
    sim_days = [h / 24.0 for h in sim["time_h"]]
    has_wound, has_fibroblast, has_tumor, has_microenv, has_ph, has_ra = detect_modules(sim)

    if not has_wound and not has_tumor and not has_ra:
        print("No wound, tumor, or RA data found in metrics. Nothing to validate.")
        sys.exit(0)

    # --- Compute once ---
    wound_r = validate_wound(sim, sim_days, condition) if has_wound else None
    fibro_r = validate_fibroblast(sim, sim_days) if has_fibroblast else None
    micro_r = validate_microenvironment(sim, sim_days, condition) if has_microenv else None
    ph_r = validate_ph(sim, sim_days) if has_ph else None
    tumor_r = validate_tumor(sim, sim_days) if has_tumor else None
    ra_r = validate_ra(sim, sim_days) if has_ra else None

    # --- Print once ---
    print_summary(wound=wound_r, fibroblast=fibro_r, tumor=tumor_r,
                  microenv=micro_r, ph=ph_r, ra=ra_r)

    out_dir = plots_dir(sim_path)
    os.makedirs(out_dir, exist_ok=True)

    # --- Per-module PNGs ---
    if wound_r:
        fig, axes = plt.subplots(2, 2, figsize=(12, 8))
        plot_wound_panels(wound_r, sim_days, axes)
        wound_tag = f" [{condition}]" if condition != "normal" else ""
        fig.suptitle(f"Wound Healing Validation{wound_tag}", fontsize=14, fontweight="bold")
        fig.tight_layout(rect=[0, 0, 1, 0.96])
        fig.savefig(os.path.join(out_dir, "wound_validation.png"), dpi=150)
        plt.close(fig)

    if fibro_r:
        fig, axes = plt.subplots(3, 1, figsize=(8, 9), sharex=True)
        plot_fibroblast_panels(fibro_r, sim_days, axes)
        fig.suptitle("Fibroblast/Collagen Validation", fontsize=14, fontweight="bold")
        fig.tight_layout(rect=[0, 0, 1, 0.96])
        fig.savefig(os.path.join(out_dir, "fibroblast_validation.png"), dpi=150)
        plt.close(fig)

    if micro_r:
        mr = 3 if ph_r else 2
        fig, axes = plt.subplots(mr, 2, figsize=(12, 4 * mr))
        plot_microenvironment_panels(micro_r, sim_days, axes[:2])
        if ph_r:
            plot_ph_panel(ph_r, sim_days, axes[2, 0])
            axes[2, 1].set_visible(False)
        fig.suptitle("Microenvironment Validation", fontsize=14, fontweight="bold")
        fig.tight_layout(rect=[0, 0, 1, 0.96])
        fig.savefig(os.path.join(out_dir, "microenvironment_validation.png"), dpi=150)
        plt.close(fig)

    if tumor_r:
        fig, axes = plt.subplots(1, 2, figsize=(12, 4))
        plot_tumor_panels(tumor_r, sim_days, axes)
        fig.suptitle("Tumor Growth Validation", fontsize=14, fontweight="bold")
        fig.tight_layout(rect=[0, 0, 1, 0.92])
        fig.savefig(os.path.join(out_dir, "tumor_validation.png"), dpi=150)
        plt.close(fig)

    if ra_r:
        has_ext = ra_r.get("has_bone") or ra_r.get("has_tcell") or ra_r.get("has_syn")
        ra_rows = 3 if has_ext else 3
        ra_cols = 2 if has_ext else 1
        if has_ext:
            fig, ax_arr = plt.subplots(3, 2, figsize=(12, 9), sharex=True)
            ra_axes = [ax_arr[0, 0], ax_arr[0, 1], ax_arr[1, 0],
                       ax_arr[1, 1], ax_arr[2, 0], ax_arr[2, 1]]
        else:
            fig, ra_axes = plt.subplots(3, 1, figsize=(8, 9), sharex=True)
        plot_ra_panels(ra_r, sim_days, ra_axes)
        fig.suptitle("Rheumatoid Arthritis Validation", fontsize=14, fontweight="bold")
        fig.tight_layout(rect=[0, 0, 1, 0.96])
        fig.savefig(os.path.join(out_dir, "ra_validation.png"), dpi=150)
        plt.close(fig)

    # --- Dashboard PNG (adaptive rows) ---
    nrows = 0
    if has_wound:
        nrows += 2
    if has_fibroblast:
        nrows += 2  # fibroblast now has 3 panels across 2 rows
    if has_microenv:
        nrows += 2
    if has_ph:
        nrows += 1
    if has_tumor:
        nrows += 1
    if has_ra:
        ra_ext = (ra_r and (ra_r.get("has_bone") or ra_r.get("has_tcell")
                            or ra_r.get("has_syn")))
        nrows += 3 if ra_ext else 2  # 6 panels (3x2) or 3 panels (2 rows)

    if nrows > 0:
        fig, axes = plt.subplots(nrows, 2, figsize=(12, 4 * nrows))
        if nrows == 1:
            axes = axes.reshape(1, 2)

        row = 0
        if wound_r:
            plot_wound_panels(wound_r, sim_days, axes[row:row+2])
            row += 2
        if fibro_r:
            # 3 panels: spread across row[0..1] and row[2] (left)
            fibro_axes = [axes[row, 0], axes[row, 1], axes[row+1, 0]]
            plot_fibroblast_panels(fibro_r, sim_days, fibro_axes)
            axes[row+1, 1].set_visible(False)
            row += 2
        if micro_r:
            plot_microenvironment_panels(micro_r, sim_days, axes[row:row+2])
            row += 2
        if ph_r:
            plot_ph_panel(ph_r, sim_days, axes[row, 0])
            axes[row, 1].set_visible(False)
            row += 1
        if tumor_r:
            plot_tumor_panels(tumor_r, sim_days, axes[row])
            row += 1
        if ra_r:
            if ra_ext:
                ra_axes = [axes[row, 0], axes[row, 1], axes[row+1, 0],
                           axes[row+1, 1], axes[row+2, 0], axes[row+2, 1]]
                plot_ra_panels(ra_r, sim_days, ra_axes)
                row += 3
            else:
                ra_axes = [axes[row, 0], axes[row, 1], axes[row+1, 0]]
                plot_ra_panels(ra_r, sim_days, ra_axes)
                axes[row+1, 1].set_visible(False)
                row += 2

        max_day = max(sim_days) if sim_days else 30
        if has_tumor:
            for c in range(2):
                axes[-1, c].set_xlim(0, max(30, max_day))

        cond_tag = f" [{condition}]" if condition != "normal" else ""
        fig.suptitle(f"skibidy Validation Dashboard{cond_tag}", fontsize=14, fontweight="bold")
        fig.tight_layout(rect=[0, 0, 1, 0.96])
        path = os.path.join(out_dir, "validation_dashboard.png")
        fig.savefig(path, dpi=150)
        plt.close(fig)
        print(f"  Saved {path}")


if __name__ == "__main__":
    main()
