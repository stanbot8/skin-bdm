#!/usr/bin/env python3
"""Resilience metrics for wound healing treatment trajectories.

Computes perturbation-response metrics from simulation timeseries, drawing on
established frameworks from ecological resilience theory (Holling 1973),
pharmacodynamics (PK/PD), and control theory:

  - Maximum deviation (Holling): peak displacement from homeostasis
  - Resistance (Van Meerbeek et al. 2021): inverse of max deviation
  - Return time / elasticity: time to recover from max deviation
  - AUC burden (PK/PD): area under curve for cumulative exposure
  - Overshoot & settling time (control theory): transient response metrics

Usage:
    python3 scripts/analysis/resilience_analysis.py output/treatment_study/
    python3 scripts/analysis/resilience_analysis.py output/treatment_study/ --baseline=metrics_baseline.csv
"""

import argparse
import csv
import glob
import os
import sys

# Simulation time step in hours
DT_HOURS = 0.1
HOURS_PER_DAY = 24.0

# Metrics to analyze and their reference (homeostatic) values.
# "direction" indicates whether healing moves the value up or down.
#   up   = healthy state is HIGH  (closure, collagen, perfusion)
#   down = healthy state is LOW   (inflammation, MMP, neutrophils)
TRACKED_METRICS = {
    "wound_closure_pct": {
        "label": "Wound closure",
        "unit": "%",
        "reference": 100.0,     # fully closed
        "direction": "up",
    },
    "mean_infl_wound": {
        "label": "Inflammation",
        "unit": "a.u.",
        "reference": 0.0,       # no inflammation = homeostasis
        "direction": "down",
    },
    "mean_collagen_wound": {
        "label": "Collagen",
        "unit": "a.u.",
        "reference": None,      # no fixed reference; use peak
        "direction": "up",
    },
    "mean_perfusion_wound": {
        "label": "Perfusion",
        "unit": "a.u.",
        "reference": 1.0,       # basal perfusion
        "direction": "up",
    },
    "scar_magnitude": {
        "label": "Scar",
        "unit": "a.u.",
        "reference": 0.0,       # no scar = ideal
        "direction": "down",
    },
    "n_neutrophils": {
        "label": "Neutrophils",
        "unit": "cells",
        "reference": 0.0,       # resolved = no neutrophils
        "direction": "down",
    },
    "mean_mmp_wound": {
        "label": "MMP",
        "unit": "a.u.",
        "reference": 0.0,
        "direction": "down",
    },
    "mean_fibronectin_wound": {
        "label": "Fibronectin",
        "unit": "a.u.",
        "reference": None,
        "direction": "up",
    },
}


def load_timeseries(csv_path):
    """Load a metrics CSV into a list of dicts with float values."""
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        rows = []
        for row in reader:
            parsed = {}
            for k, v in row.items():
                try:
                    parsed[k] = float(v)
                except (ValueError, TypeError):
                    parsed[k] = 0.0
            rows.append(parsed)
    return rows


def step_to_day(step):
    return step * DT_HOURS / HOURS_PER_DAY


def compute_resilience(rows, metric_key, metric_info):
    """Compute resilience metrics for a single timeseries.

    Returns a dict of resilience metrics, or None if the metric is absent.
    """
    if not rows or metric_key not in rows[0]:
        return None

    values = [r[metric_key] for r in rows]
    steps = [r.get("step", i) for i, r in enumerate(rows)]
    n = len(values)
    if n < 2:
        return None

    direction = metric_info["direction"]
    ref = metric_info["reference"]

    # --- Maximum deviation (ecological resilience) ---
    # For "down" metrics (inflammation): max deviation = peak value above ref
    # For "up" metrics (closure): max deviation = deepest trough below ref
    if direction == "down":
        max_dev_val = max(values)
        max_dev_idx = values.index(max_dev_val)
        max_deviation = max_dev_val - (ref if ref is not None else 0.0)
    else:
        if ref is not None:
            # deviation = how far below reference
            min_val = min(values)
            max_dev_idx = values.index(min_val)
            max_deviation = ref - min_val
        else:
            # no reference; use range
            max_dev_val = max(values)
            min_val = min(values)
            max_dev_idx = values.index(min_val)
            max_deviation = max_dev_val - min_val

    max_dev_day = step_to_day(steps[max_dev_idx])

    # --- Resistance (inverse of max deviation; Van Meerbeek et al. 2021) ---
    # Clamp to avoid nonsensical values when deviation is near-zero
    resistance = 1.0 / max(max_deviation, 0.01)

    # --- Peak response (C_max / nadir in PK/PD) ---
    if direction == "down":
        peak_response = max(values)
        peak_idx = values.index(peak_response)
    else:
        peak_response = max(values)
        peak_idx = values.index(peak_response)

    t_max = step_to_day(steps[peak_idx])  # T_max (time to peak)

    # --- Return time (time from max deviation to within 10% of final value) ---
    final_val = values[-1]
    tolerance = abs(final_val) * 0.10 + 1e-10  # 10% tolerance band

    return_time = None
    if max_dev_idx < n - 1:
        for i in range(max_dev_idx, n):
            if abs(values[i] - final_val) <= tolerance:
                return_time = step_to_day(steps[i]) - max_dev_day
                break

    # --- Elasticity (rate of recovery from max deviation) ---
    if return_time and return_time > 0:
        elasticity = max_deviation / return_time  # units per day
    else:
        elasticity = None

    # --- Settling time (control theory: time to stay within 5% of final) ---
    settle_tol = abs(final_val) * 0.05 + 1e-10
    settling_time = None
    for i in range(n - 1, -1, -1):
        if abs(values[i] - final_val) > settle_tol:
            if i < n - 1:
                settling_time = step_to_day(steps[i + 1])
            break

    # --- Overshoot (control theory) ---
    # For "down" metrics: overshoot = peak above steady-state
    # For "up" metrics: undershoot = deepest trough below final
    # Use max_deviation as denominator to avoid division-near-zero artifacts
    if direction == "down":
        overshoot = max(values) - final_val
        denom = max(max_deviation, 0.01)
        overshoot_pct = (overshoot / denom) * 100
    else:
        if ref is not None:
            overshoot = final_val - min(values)
            denom = max(max_deviation, 0.01)
            overshoot_pct = (overshoot / denom) * 100
        else:
            overshoot = max(values) - final_val
            denom = max(max(values), 0.01)
            overshoot_pct = (overshoot / denom) * 100

    # --- AUC burden (PK/PD: area under the curve) ---
    # For "down" metrics: AUC above reference = cumulative inflammatory burden
    # For "up" metrics: AUC below reference = cumulative healing deficit
    auc = 0.0
    dt_days = DT_HOURS / HOURS_PER_DAY
    step_interval = steps[1] - steps[0] if n > 1 else 1
    dt_sample = step_interval * dt_days

    for i in range(1, n):
        if direction == "down":
            val_ref = ref if ref is not None else 0.0
            # trapezoidal integration of excess above reference
            y0 = max(values[i - 1] - val_ref, 0)
            y1 = max(values[i] - val_ref, 0)
        else:
            if ref is not None:
                # deficit below reference
                y0 = max(ref - values[i - 1], 0)
                y1 = max(ref - values[i], 0)
            else:
                y0 = values[i - 1]
                y1 = values[i]
        auc += 0.5 * (y0 + y1) * dt_sample

    # --- Steady-state error (control theory) ---
    if ref is not None:
        if direction == "down":
            ss_error = final_val - ref  # positive = still elevated
        else:
            ss_error = ref - final_val  # positive = still below target
    else:
        ss_error = None

    # --- Rise time (control theory: time from 10% to 90% of final response) ---
    rise_time = None
    if direction == "up" and final_val > values[0]:
        target_range = final_val - values[0]
        t10 = None
        t90 = None
        for i, v in enumerate(values):
            progress = (v - values[0]) / max(target_range, 1e-10)
            if t10 is None and progress >= 0.10:
                t10 = step_to_day(steps[i])
            if t90 is None and progress >= 0.90:
                t90 = step_to_day(steps[i])
                break
        if t10 is not None and t90 is not None:
            rise_time = t90 - t10

    return {
        "max_deviation": max_deviation,
        "max_deviation_day": max_dev_day,
        "resistance": resistance,
        "peak_response": peak_response,
        "t_max_days": t_max,
        "return_time_days": return_time,
        "elasticity": elasticity,
        "settling_time_days": settling_time,
        "overshoot": overshoot,
        "overshoot_pct": overshoot_pct,
        "auc_burden": auc,
        "steady_state_error": ss_error,
        "rise_time_days": rise_time,
        "final_value": final_val,
        "initial_value": values[0],
    }


def analyze_treatment(csv_path):
    """Compute resilience metrics for all tracked metrics in a run."""
    rows = load_timeseries(csv_path)
    if not rows:
        return {}

    results = {}
    for key, info in TRACKED_METRICS.items():
        r = compute_resilience(rows, key, info)
        if r:
            results[key] = r
    return results


def fmt(val, precision=3):
    if val is None:
        return "N/A"
    if isinstance(val, float):
        return f"{val:.{precision}f}"
    return str(val)


def print_report(all_results, baseline_name="baseline"):
    """Print a comparative resilience report."""
    treatments = list(all_results.keys())
    if not treatments:
        print("No results.")
        return

    baseline = all_results.get(baseline_name, {})

    # Per-metric comparison tables
    for metric_key, metric_info in TRACKED_METRICS.items():
        label = metric_info["label"]
        direction = metric_info["direction"]

        # Check if any treatment has this metric
        has_metric = any(metric_key in all_results[t] for t in treatments)
        if not has_metric:
            continue

        print(f"\n{'='*90}")
        print(f"  {label} ({metric_key})")
        print(f"  Direction: {direction} is healthy | Reference: {metric_info['reference']}")
        print(f"{'='*90}")

        header = (f"{'Treatment':<20} {'Final':>8} {'MaxDev':>8} {'Resist':>8} "
                  f"{'Peak':>8} {'Tmax(d)':>8} {'RetTime':>8} "
                  f"{'Settle':>8} {'AUC':>10} {'OS%':>8}")
        print(header)
        print("-" * len(header))

        for t in treatments:
            r = all_results[t].get(metric_key)
            if not r:
                print(f"{t:<20} (no data)")
                continue

            print(f"{t:<20} "
                  f"{fmt(r['final_value'], 2):>8} "
                  f"{fmt(r['max_deviation'], 3):>8} "
                  f"{fmt(r['resistance'], 1):>8} "
                  f"{fmt(r['peak_response'], 3):>8} "
                  f"{fmt(r['t_max_days'], 1):>8} "
                  f"{fmt(r['return_time_days'], 1):>8} "
                  f"{fmt(r['settling_time_days'], 1):>8} "
                  f"{fmt(r['auc_burden'], 3):>10} "
                  f"{fmt(r['overshoot_pct'], 0):>8}")

    # --- Composite resilience ranking ---
    print(f"\n{'='*90}")
    print("  COMPOSITE RESILIENCE RANKING")
    print(f"{'='*90}")
    print("  Lower rank = better healing trajectory (less deviation, faster recovery, lower burden)\n")

    scores = {}
    rank_metrics = [
        ("wound_closure_pct", "auc_burden", False),       # lower AUC deficit = better
        ("mean_infl_wound", "auc_burden", False),          # lower inflammatory burden = better
        ("mean_infl_wound", "return_time_days", False),    # faster resolution = better
        ("mean_infl_wound", "max_deviation", False),       # lower peak = better
        ("scar_magnitude", "final_value", False),          # lower scar = better
    ]

    for t in treatments:
        rank_sum = 0
        n_ranked = 0
        for metric_key, res_key, higher_is_better in rank_metrics:
            r = all_results[t].get(metric_key)
            if not r or r.get(res_key) is None:
                continue
            val = r[res_key]
            # Collect (treatment, value) for ranking
            all_vals = []
            for t2 in treatments:
                r2 = all_results[t2].get(metric_key)
                if r2 and r2.get(res_key) is not None:
                    all_vals.append(r2[res_key])
            all_vals.sort(reverse=higher_is_better)
            rank = all_vals.index(val) + 1
            rank_sum += rank
            n_ranked += 1

        if n_ranked > 0:
            scores[t] = rank_sum / n_ranked

    ranked = sorted(scores.items(), key=lambda x: x[1])
    print(f"  {'Rank':<6} {'Treatment':<20} {'Mean Rank Score':>15}")
    print(f"  {'-'*45}")
    for i, (t, score) in enumerate(ranked, 1):
        marker = " <-- baseline" if t == baseline_name else ""
        print(f"  {i:<6} {t:<20} {score:>15.2f}{marker}")

    # --- Clinical insights ---
    print(f"\n{'='*90}")
    print("  TRAJECTORY INSIGHTS")
    print(f"{'='*90}")

    if baseline:
        bl_infl = baseline.get("mean_infl_wound", {})
        bl_closure = baseline.get("wound_closure_pct", {})

        for t in treatments:
            if t == baseline_name:
                continue
            r_infl = all_results[t].get("mean_infl_wound", {})
            r_closure = all_results[t].get("wound_closure_pct", {})

            insights = []

            # Inflammatory trajectory
            if r_infl and bl_infl:
                bl_peak = bl_infl.get("peak_response", 0)
                t_peak = r_infl.get("peak_response", 0)
                if bl_peak > 0:
                    reduction = (1 - t_peak / bl_peak) * 100
                    if reduction > 20:
                        insights.append(f"peak inflammation reduced {reduction:.0f}%")
                    elif reduction < -10:
                        insights.append(f"WARNING: peak inflammation increased {-reduction:.0f}%")

                bl_rt = bl_infl.get("return_time_days")
                t_rt = r_infl.get("return_time_days")
                if bl_rt and t_rt:
                    if t_rt < bl_rt * 0.7:
                        insights.append(f"inflammation resolves {(1-t_rt/bl_rt)*100:.0f}% faster")
                    elif t_rt > bl_rt * 1.3:
                        insights.append(f"inflammation resolution {(t_rt/bl_rt-1)*100:.0f}% slower")

                bl_auc = bl_infl.get("auc_burden", 0)
                t_auc = r_infl.get("auc_burden", 0)
                if bl_auc > 0:
                    auc_reduction = (1 - t_auc / bl_auc) * 100
                    if auc_reduction > 20:
                        insights.append(f"cumulative inflammatory burden reduced {auc_reduction:.0f}%")

            # Closure trajectory
            if r_closure and bl_closure:
                bl_rise = bl_closure.get("rise_time_days")
                t_rise = r_closure.get("rise_time_days")
                if bl_rise and t_rise and bl_rise > 0:
                    if t_rise < bl_rise * 0.8:
                        insights.append(f"closure rise time {(1-t_rise/bl_rise)*100:.0f}% faster")

            # Scar
            r_scar = all_results[t].get("scar_magnitude", {})
            bl_scar = baseline.get("scar_magnitude", {})
            if r_scar and bl_scar:
                bl_final = bl_scar.get("final_value", 0)
                t_final = r_scar.get("final_value", 0)
                if bl_final > 0:
                    scar_reduction = (1 - t_final / bl_final) * 100
                    if scar_reduction > 20:
                        insights.append(f"scarring reduced {scar_reduction:.0f}%")

            if insights:
                print(f"\n  {t}:")
                for ins in insights:
                    print(f"    - {ins}")

    print()


def save_csv(all_results, output_path):
    """Save all resilience metrics to a flat CSV."""
    rows = []
    for treatment, metrics in all_results.items():
        for metric_key, r in metrics.items():
            row = {"treatment": treatment, "metric": metric_key}
            row.update(r)
            rows.append(row)

    if not rows:
        return

    keys = list(rows[0].keys())
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Resilience metrics saved to {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Resilience analysis of wound healing treatment trajectories")
    parser.add_argument("study_dir",
                        help="Directory containing metrics_*.csv files")
    parser.add_argument("--baseline", default="metrics_baseline.csv",
                        help="Baseline metrics CSV filename")
    args = parser.parse_args()

    study_dir = args.study_dir
    pattern = os.path.join(study_dir, "metrics_*.csv")
    csv_files = sorted(glob.glob(pattern))

    if not csv_files:
        print(f"No metrics_*.csv files found in {study_dir}")
        sys.exit(1)

    all_results = {}
    baseline_name = None

    for csv_path in csv_files:
        fname = os.path.basename(csv_path)
        # Extract treatment name from metrics_<name>.csv
        name = fname.replace("metrics_", "").replace(".csv", "")
        if fname == args.baseline:
            baseline_name = name

        print(f"Analyzing {name}...", flush=True)
        all_results[name] = analyze_treatment(csv_path)

    if baseline_name is None and "baseline" in all_results:
        baseline_name = "baseline"

    print_report(all_results, baseline_name=baseline_name or "baseline")
    save_csv(all_results, os.path.join(study_dir, "resilience_metrics.csv"))


if __name__ == "__main__":
    main()
