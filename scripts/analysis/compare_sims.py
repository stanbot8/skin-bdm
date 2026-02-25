#!/usr/bin/env python3
"""Compare normal vs diabetic wound simulation metrics.

Reads two metrics CSV files and prints a side-by-side comparison of key
healing indicators at selected time points.

Usage:
    python3 scripts/analysis/compare_sims.py output/metrics_normal.csv output/metrics_diabetic.csv
"""

import csv
import sys


def load_metrics(path):
    """Load metrics CSV into list of dicts."""
    with open(path) as f:
        reader = csv.DictReader(f)
        return list(reader)


def find_row(rows, step):
    """Find row closest to given step."""
    best = None
    for r in rows:
        s = int(r["step"])
        if best is None or abs(s - step) < abs(int(best["step"]) - step):
            best = r
    return best


def closure_step(rows, threshold=99.0):
    """Find first step where wound_closure_pct >= threshold."""
    for r in rows:
        if float(r["wound_closure_pct"]) >= threshold:
            return int(r["step"])
    return None


def peak_value(rows, col):
    """Find peak value and step for a column."""
    best_val = 0
    best_step = 0
    for r in rows:
        v = float(r[col])
        if v > best_val:
            best_val = v
            best_step = int(r["step"])
    return best_val, best_step


def main():
    if len(sys.argv) < 3:
        print("Usage: python3 scripts/analysis/compare_sims.py <normal.csv> <diabetic.csv>")
        sys.exit(1)

    normal = load_metrics(sys.argv[1])
    diabetic = load_metrics(sys.argv[2])

    print("=" * 72)
    print("WOUND HEALING COMPARISON: Normal vs Diabetic")
    print("=" * 72)

    # Closure kinetics
    nc = closure_step(normal)
    dc = closure_step(diabetic)
    print(f"\n--- Wound Closure ---")
    print(f"  Normal:   {'step ' + str(nc) + f' ({nc*0.1:.0f}h = {nc*0.1/24:.1f}d)' if nc else 'NOT CLOSED'}")
    print(f"  Diabetic: {'step ' + str(dc) + f' ({dc*0.1:.0f}h = {dc*0.1/24:.1f}d)' if dc else 'NOT CLOSED'}")
    if nc and dc:
        print(f"  Delay:    {(dc - nc) * 0.1:.0f}h ({(dc - nc) * 0.1 / 24:.1f} days)")
    elif nc and not dc:
        last_d = diabetic[-1]
        print(f"  Diabetic wound FAILED TO CLOSE (last: {float(last_d['wound_closure_pct']):.1f}% at step {last_d['step']})")

    # Key timepoints comparison
    timepoints = [500, 1000, 2000, 3000, 5000, 7000]
    cols = [
        ("wound_closure_pct", "Closure %", ".1f"),
        ("mean_infl_wound", "Inflammation", ".4f"),
        ("n_neutrophils", "Neutrophils", ".0f"),
        ("n_macrophages", "Macrophages", ".0f"),
        ("n_fibroblasts", "Fibroblasts", ".0f"),
        ("n_myofibroblasts", "Myofibroblasts", ".0f"),
        ("mean_tgfb_wound", "TGF-beta", ".4f"),
        ("mean_collagen_wound", "Collagen", ".4f"),
        ("mean_mmp_wound", "MMP", ".4f"),
        ("mean_fibronectin_wound", "Fibronectin", ".4f"),
        ("mean_perfusion_wound", "Perfusion", ".4f"),
        ("scar_magnitude", "Scar", ".4f"),
    ]

    print(f"\n--- Time Series Comparison ---")
    for col_name, label, fmt in cols:
        print(f"\n  {label}:")
        print(f"    {'Step':>6}  {'Time':>6}  {'Normal':>10}  {'Diabetic':>10}  {'Ratio':>8}")
        for step in timepoints:
            nr = find_row(normal, step)
            dr = find_row(diabetic, step)
            if nr and dr and col_name in nr:
                nv = float(nr[col_name])
                dv = float(dr[col_name])
                ratio = dv / nv if nv > 1e-10 else float('inf') if dv > 1e-10 else 1.0
                t_h = step * 0.1
                print(f"    {step:>6}  {t_h:>5.0f}h  {nv:>10{fmt}}  {dv:>10{fmt}}  {ratio:>7.2f}x")

    # Peak analysis
    print(f"\n--- Peak Values ---")
    peak_cols = [
        ("mean_infl_wound", "Inflammation"),
        ("mean_mmp_wound", "MMP"),
        ("mean_tgfb_wound", "TGF-beta"),
        ("mean_collagen_wound", "Collagen"),
        ("mean_fibronectin_wound", "Fibronectin"),
        ("scar_magnitude", "Scar"),
    ]
    print(f"  {'Metric':<16} {'Normal':>20}  {'Diabetic':>20}")
    for col_name, label in peak_cols:
        nv, ns = peak_value(normal, col_name)
        dv, ds = peak_value(diabetic, col_name)
        print(f"  {label:<16} {nv:>8.4f} @ step {ns:<5}  {dv:>8.4f} @ step {ds:<5}")

    # Biological insights
    print(f"\n--- Biological Insights ---")
    ni_peak, _ = peak_value(normal, "mean_infl_wound")
    di_peak, _ = peak_value(diabetic, "mean_infl_wound")
    nm_peak, _ = peak_value(normal, "mean_mmp_wound")
    dm_peak, _ = peak_value(diabetic, "mean_mmp_wound")
    nc_peak, _ = peak_value(normal, "mean_collagen_wound")
    dc_peak, _ = peak_value(diabetic, "mean_collagen_wound")

    if di_peak > ni_peak * 1.5:
        print(f"  [!] Diabetic inflammation {di_peak/ni_peak:.1f}x higher -- chronic inflammation phenotype")
    if dm_peak > nm_peak * 1.5:
        print(f"  [!] Diabetic MMP {dm_peak/nm_peak:.1f}x higher -- ECM degradation imbalance")
    if dc_peak < nc_peak * 0.7:
        print(f"  [!] Diabetic collagen {dc_peak/nc_peak:.1f}x lower -- impaired matrix deposition")
    if not dc:
        print(f"  [!] Diabetic wound failed to close -- chronic wound phenotype")
    elif dc and nc and dc > nc * 1.3:
        print(f"  [!] Diabetic healing delayed by {(dc-nc)*0.1/24:.1f} days")

    print()


if __name__ == "__main__":
    main()
