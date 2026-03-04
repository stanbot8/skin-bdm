#!/usr/bin/env python3
"""Standalone single-module validation.

Usage:
    python3 literature/validators/compare.py wound [path/to/metrics.csv] [--diabetic|--normal|...]
    python3 literature/validators/compare.py fibroblast [path/to/metrics.csv]
    python3 literature/validators/compare.py microenvironment [path/to/metrics.csv]
    python3 literature/validators/compare.py tumor [path/to/metrics.csv]
    python3 literature/validators/compare.py immune [path/to/metrics.csv] [--diabetic|--normal]
    python3 literature/validators/compare.py ra [path/to/metrics.csv]
"""

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from lib import (load_csv, plots_dir, detect_condition, surface_fraction,
                 validate_wound, validate_fibroblast, validate_tumor,
                 validate_microenvironment, validate_ph, validate_ra,
                 plot_wound_panels, plot_fibroblast_panels, plot_tumor_panels,
                 plot_microenvironment_panels, plot_ph_panel, plot_ra_panels,
                 print_summary, SIM_KW, REF_KW)


def _parse_args():
    positional = [a for a in sys.argv[1:] if not a.startswith("--")]
    flags = [a for a in sys.argv[1:] if a.startswith("--")]
    if not positional:
        print(__doc__)
        sys.exit(1)
    module = positional[0]
    sim_path = positional[1] if len(positional) > 1 else "output/skibidy/metrics.csv"
    condition = None
    for f in flags:
        key = f.lstrip("-")
        if key in ("diabetic", "normal", "burn", "pressure", "surgical", "rheumatoid"):
            condition = key
    if condition is None:
        condition = detect_condition()
    return module, sim_path, condition


def _load_sim(sim_path):
    if not os.path.exists(sim_path):
        print(f"Error: {sim_path} not found.")
        sys.exit(1)
    sim = load_csv(sim_path)
    sim_days = [h / 24.0 for h in sim["time_h"]]
    return sim, sim_days


def _save(fig, out_dir, filename):
    os.makedirs(out_dir, exist_ok=True)
    fig.savefig(os.path.join(out_dir, filename), dpi=150)
    plt.close(fig)
    print(f"  Saved {out_dir}/{filename}")


def run_wound(sim, sim_days, condition, out_dir):
    if "wound_closure_pct" not in sim or max(sim["wound_closure_pct"]) == 0:
        print("No wound data found. Set wound_enabled = true and re-run.")
        sys.exit(1)
    r = validate_wound(sim, sim_days, condition)
    print_summary(wound=r)
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    plot_wound_panels(r, sim_days, axes)
    fig.suptitle("Wound Healing Validation", fontsize=14, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    _save(fig, out_dir, "wound_validation.png")


def run_fibroblast(sim, sim_days, _condition, out_dir):
    if "n_myofibroblasts" not in sim or max(sim["n_myofibroblasts"]) == 0:
        print("No fibroblast data found. Set fibroblast_enabled = true and re-run.")
        sys.exit(1)
    r = validate_fibroblast(sim, sim_days)
    print_summary(fibroblast=r)
    fig, axes = plt.subplots(3, 1, figsize=(8, 9), sharex=True)
    plot_fibroblast_panels(r, sim_days, axes)
    fig.suptitle("Fibroblast/Collagen Validation", fontsize=14, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    _save(fig, out_dir, "fibroblast_validation.png")


def run_microenvironment(sim, sim_days, _condition, out_dir):
    if "mean_tgfb_wound" not in sim or max(sim["mean_tgfb_wound"]) == 0:
        print("No microenvironment data found.")
        sys.exit(1)
    r = validate_microenvironment(sim, sim_days)
    ph_r = None
    if "mean_ph_wound" in sim and max(sim["mean_ph_wound"]) > 0:
        ph_r = validate_ph(sim, sim_days)
    print_summary(microenv=r, ph=ph_r)
    fig, axes = plt.subplots(3, 2, figsize=(12, 12))
    plot_microenvironment_panels(r, sim_days, axes[:2])
    if ph_r:
        plot_ph_panel(ph_r, sim_days, axes[2, 0])
    else:
        axes[2, 0].set_visible(False)
    axes[2, 1].set_visible(False)
    fig.suptitle("Microenvironment Validation", fontsize=14, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    _save(fig, out_dir, "microenvironment_validation.png")


def run_tumor(sim, sim_days, _condition, out_dir):
    if "n_tumor_cells" not in sim:
        print("No tumor data found. Set tumor_enabled = true and re-run.")
        sys.exit(1)
    r = validate_tumor(sim, sim_days)
    print_summary(tumor=r)
    n_init = r["n_init"]
    n_final = r["obs_final"]
    print(f"\n  Scale context (Gompertzian):")
    print(f"    Initial cells:              {n_init:.0f}  (surface fraction {surface_fraction(n_init):.0%})")
    print(f"    Final cells:                {n_final:.0f}  (surface fraction {surface_fraction(n_final):.0%})")
    print(f"    Established BCC (~10^5):    surface fraction {surface_fraction(1e5):.0%}")
    print(f"    Reference Ki-67 ({r['bcc_ki67_pct']:.1f}%) applies at established-tumor scale.")
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    plot_tumor_panels(r, sim_days, axes)
    fig.suptitle("Tumor Growth Validation", fontsize=14, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    _save(fig, out_dir, "tumor_validation.png")


def run_immune(sim, sim_days, condition, out_dir):
    if "wound_closure_pct" not in sim or max(sim["wound_closure_pct"]) == 0:
        print("No wound data found. Set wound_enabled = true and re-run.")
        sys.exit(1)
    r = validate_wound(sim, sim_days, condition)
    print(f"Immune cell kinetics validation ({len(sim_days)} sim points vs literature)")
    print(f"  Neutrophils:  RMSE = {r['neut_rmse'] * 100:.2f} %  (peak count = {r['neut_peak']:.0f})")
    print(f"  Macrophages:  RMSE = {r['mac_rmse'] * 100:.2f} %  (peak count = {r['mac_peak']:.0f})")
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    ax1.plot(sim_days, r["sim_neut"], **SIM_KW)
    ax1.plot(r["ref_immune"]["day"], r["ref_immune"]["neutrophils_normalized"], **REF_KW)
    ax1.set_ylabel("Neutrophils (normalized)")
    ax1.set_title("Immune Cell Kinetics: Simulation vs Literature")
    ax1.set_ylim(-0.05, 1.15)
    ax1.legend(loc="upper right")
    ax1.grid(True, alpha=0.3)
    ax1.text(0.98, 0.85, f"RMSE = {r['neut_rmse'] * 100:.1f}%",
             transform=ax1.transAxes, ha="right", va="top", fontsize=9, color="gray")
    ax2.plot(sim_days, r["sim_mac"], **SIM_KW)
    ax2.plot(r["ref_immune"]["day"], r["ref_immune"]["macrophages_normalized"], **REF_KW)
    ax2.set_xlabel("Time (days)")
    ax2.set_ylabel("Macrophages (normalized)")
    ax2.set_xlim(0, 30)
    ax2.set_ylim(-0.05, 1.15)
    ax2.legend(loc="upper right")
    ax2.grid(True, alpha=0.3)
    ax2.text(0.98, 0.85, f"RMSE = {r['mac_rmse'] * 100:.1f}%",
             transform=ax2.transAxes, ha="right", va="top", fontsize=9, color="gray")
    fig.tight_layout()
    _save(fig, out_dir, "immune_validation.png")


def run_ra(sim, sim_days, _condition, out_dir):
    if "mean_tnf_alpha_wound" not in sim or max(sim["mean_tnf_alpha_wound"]) == 0:
        print("No RA data found. Set ra_enabled = true and re-run.")
        sys.exit(1)
    r = validate_ra(sim, sim_days)
    print(f"RA validation ({len(sim_days)} sim points vs literature)")
    print(f"  TNF-alpha:  RMSE = {r['tnf_rmse'] * 100:.2f} %  (peak = {r['tnf_peak']:.4f})")
    print(f"    Flare 0-7d:   RMSE = {r['tnf_flare_rmse'] * 100:.2f} %")
    print(f"    Chronic 7-30d: RMSE = {r['tnf_chronic_rmse'] * 100:.2f} %")
    print(f"  IL-6:       RMSE = {r['il6_rmse'] * 100:.2f} %  (peak = {r['il6_peak']:.4f})")
    print(f"  Cartilage:  RMSE = {r['cart_rmse'] * 100:.2f} %")
    if r.get("has_bone"):
        print(f"  Bone:       RMSE = {r['bone_rmse'] * 100:.2f} %")
    if r.get("has_tcell"):
        print(f"  T cells:    RMSE = {r['tcell_rmse'] * 100:.2f} %  (peak = {r['tcell_peak']:.4f})")
    if r.get("has_syn"):
        print(f"  Pannus:     RMSE = {r['syn_rmse'] * 100:.2f} %  (peak = {r['syn_peak']:.4f})")
    has_ext = r.get("has_bone") or r.get("has_tcell") or r.get("has_syn")
    if has_ext:
        fig, ax_arr = plt.subplots(3, 2, figsize=(12, 9), sharex=True)
        axes = [ax_arr[0, 0], ax_arr[0, 1], ax_arr[1, 0],
                ax_arr[1, 1], ax_arr[2, 0], ax_arr[2, 1]]
    else:
        fig, axes = plt.subplots(3, 1, figsize=(8, 9), sharex=True)
    plot_ra_panels(r, sim_days, axes)
    fig.suptitle("Rheumatoid Arthritis Validation", fontsize=14, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    _save(fig, out_dir, "ra_validation.png")


MODULES = {
    "wound": run_wound,
    "fibroblast": run_fibroblast,
    "microenvironment": run_microenvironment,
    "microenv": run_microenvironment,
    "tumor": run_tumor,
    "immune": run_immune,
    "ra": run_ra,
}


def main():
    module, sim_path, condition = _parse_args()
    if module not in MODULES:
        print(f"Unknown module: {module}")
        print(f"Available: {', '.join(sorted(set(MODULES.values().__class__(MODULES.keys()))))}")
        sys.exit(1)
    sim, sim_days = _load_sim(sim_path)
    out_dir = plots_dir(sim_path)
    MODULES[module](sim, sim_days, condition, out_dir)


if __name__ == "__main__":
    main()
