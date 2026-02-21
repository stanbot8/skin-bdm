#!/usr/bin/env python3
"""Adaptive surrogate-guided treatment combination search.

Instead of brute-forcing all 2^N-1 combinations (255 configs x 10 runs each),
this uses three strategies to cut runtime by ~30x:

1. Adaptive consensus: start with 3 runs, stop early if CV < threshold
2. Additive surrogate: predict combo outcomes from single-treatment effects
3. Guided selection: only run the top-K most promising combos

Phases:
  1. Baseline (diabetic untreated) with adaptive consensus
  2. Singles (each treatment alone) with adaptive consensus
  3. Surrogate fit: additive model predicts all combos from singles
  4. Selection: rank combos by composite score, pick top-K
  5. Run selected combos with adaptive consensus
  6. Synergy report: compare observed vs predicted, flag interactions

Usage:
    python3 scripts/study/adaptive_study.py
    python3 scripts/study/adaptive_study.py --top-k 10 --min-runs 3 --max-runs 5
    python3 scripts/study/adaptive_study.py --treatments=hbo,npwt,growth_factor
"""

import argparse
import csv
import itertools
import math
import os
import shutil
import sys
import time

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, REPO)

from scripts.study.treatment_study import (
    available_treatments,
    prepare_config,
    run_simulation,
    extract_outcomes,
    combo_label,
)

OUTPUT_BASE = os.path.join(REPO, "output", "adaptive_study")

# Outcome keys used for scoring
SCORE_KEYS = ["time_to_50pct_days", "scar_magnitude", "peak_inflammation"]


# ---------------------------------------------------------------------------
# Adaptive consensus
# ---------------------------------------------------------------------------

def adaptive_consensus(label, treatment_names, min_runs, max_runs, cv_threshold):
    """Run simulations until CV stabilizes or max_runs reached.

    Returns (mean_outcomes, std_outcomes, n_runs, csv_paths).
    """
    outcomes_list = []
    csv_paths = []
    raw_dir = os.path.join(OUTPUT_BASE, "raw")
    os.makedirs(raw_dir, exist_ok=True)

    for i in range(max_runs):
        prepare_config(treatment_names)
        metrics_path = run_simulation(f"{label} run {i+1}")

        if not metrics_path:
            print(f"    Run {i+1} failed, skipping", flush=True)
            continue

        outcomes = extract_outcomes(metrics_path)
        if not outcomes:
            continue

        safe_label = label.replace("+", "_").replace(" ", "_")
        dst = os.path.join(raw_dir, f"{safe_label}_run{i:03d}.csv")
        shutil.copy2(metrics_path, dst)
        csv_paths.append(dst)
        outcomes_list.append(outcomes)

        # Check early stopping after min_runs
        if len(outcomes_list) >= min_runs:
            cv = _compute_cv(outcomes_list)
            print(f"    CV after {len(outcomes_list)} runs: {cv:.3f}", flush=True)
            if cv < cv_threshold:
                print(f"    Converged (CV {cv:.3f} < {cv_threshold})", flush=True)
                break

    if not outcomes_list:
        return {}, {}, 0, []

    mean_out, std_out = _aggregate_outcomes(outcomes_list)
    return mean_out, std_out, len(outcomes_list), csv_paths


def _compute_cv(outcomes_list):
    """Compute mean coefficient of variation across score keys."""
    cvs = []
    for key in SCORE_KEYS:
        vals = [o.get(key) for o in outcomes_list if o.get(key) is not None]
        if len(vals) < 2:
            continue
        mean = sum(vals) / len(vals)
        if abs(mean) < 1e-12:
            continue
        std = math.sqrt(sum((v - mean) ** 2 for v in vals) / len(vals))
        cvs.append(std / abs(mean))
    return sum(cvs) / len(cvs) if cvs else 1.0


def _aggregate_outcomes(outcomes_list):
    """Compute mean and std of outcome dicts."""
    keys = outcomes_list[0].keys()
    mean_out = {}
    std_out = {}
    for key in keys:
        vals = [o[key] for o in outcomes_list if o.get(key) is not None]
        if not vals or not isinstance(vals[0], (int, float)):
            mean_out[key] = outcomes_list[-1].get(key)
            std_out[key] = 0
            continue
        m = sum(vals) / len(vals)
        v = sum((x - m) ** 2 for x in vals) / max(1, len(vals))
        mean_out[key] = m
        std_out[key] = math.sqrt(v)
    return mean_out, std_out


# ---------------------------------------------------------------------------
# Surrogate model
# ---------------------------------------------------------------------------

def fit_surrogate(baseline, singles):
    """Fit additive surrogate: effect(X) = outcome(X) - baseline.

    Args:
        baseline: dict of baseline outcomes
        singles: dict mapping treatment name -> outcomes dict

    Returns:
        effects: dict mapping treatment name -> dict of deltas
    """
    effects = {}
    for name, outcomes in singles.items():
        delta = {}
        for key in SCORE_KEYS:
            b = baseline.get(key)
            t = outcomes.get(key)
            if b is not None and t is not None:
                delta[key] = t - b
            else:
                delta[key] = 0
        effects[name] = delta
    return effects


def predict_combo(baseline, effects, combo):
    """Predict combo outcome: baseline + sum of individual effects."""
    predicted = {}
    for key in SCORE_KEYS:
        b = baseline.get(key, 0) or 0
        total_effect = sum(effects.get(t, {}).get(key, 0) for t in combo)
        predicted[key] = b + total_effect
    return predicted


def composite_score(outcomes, baseline):
    """Score a treatment outcome (lower is better).

    Weights: 40% time-to-50%, 30% scar, 30% peak inflammation.
    Each component is normalized relative to baseline.
    """
    score = 0
    b_t50 = baseline.get("time_to_50pct_days") or 30
    t_t50 = outcomes.get("time_to_50pct_days") or 30
    score += 0.4 * (t_t50 / max(b_t50, 0.1))

    b_scar = baseline.get("scar_magnitude") or 1
    t_scar = outcomes.get("scar_magnitude") or 1
    score += 0.3 * (t_scar / max(b_scar, 1e-6))

    b_infl = baseline.get("peak_inflammation") or 1
    t_infl = outcomes.get("peak_inflammation") or 1
    score += 0.3 * (t_infl / max(b_infl, 1e-6))

    return score


def select_combos(baseline, effects, all_combos, top_k):
    """Rank all combos by predicted composite score, return top-K."""
    scored = []
    for combo in all_combos:
        predicted = predict_combo(baseline, effects, combo)
        sc = composite_score(predicted, baseline)
        scored.append((combo, predicted, sc))
    scored.sort(key=lambda x: x[2])
    return scored[:top_k]


def compute_synergy(observed, predicted, baseline):
    """Compute synergy scores normalized to baseline scale.

    synergy = (predicted - observed) / baseline
    Positive = observed is better than additive prediction (synergistic).
    Negative = observed is worse than additive prediction (antagonistic).
    For time_to_50pct, scar, and inflammation, lower is better.
    """
    synergy = {}
    for key in SCORE_KEYS:
        obs = observed.get(key)
        pred = predicted.get(key)
        base = baseline.get(key)
        if obs is not None and pred is not None and base is not None and abs(base) > 1e-12:
            synergy[key] = (pred - obs) / abs(base)
        else:
            synergy[key] = 0
    synergy["mean"] = sum(synergy[k] for k in SCORE_KEYS) / len(SCORE_KEYS)
    return synergy


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def write_results_csv(results, path):
    """Write all results to a CSV."""
    if not results:
        return

    fieldnames = ["config", "n_runs", "type"]
    # Collect all outcome keys
    outcome_keys = set()
    for r in results:
        outcome_keys.update(r.get("outcomes", {}).keys())
        outcome_keys.update(k for k in r if k.startswith("synergy_"))
        outcome_keys.update(k for k in r if k.startswith("predicted_"))
    outcome_keys = sorted(outcome_keys)
    fieldnames += outcome_keys

    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for r in results:
            row = {
                "config": r["config"],
                "n_runs": r.get("n_runs", ""),
                "type": r.get("type", ""),
            }
            for key in outcome_keys:
                if key in r:
                    row[key] = r[key]
                elif key in r.get("outcomes", {}):
                    val = r["outcomes"][key]
                    row[key] = f"{val:.6g}" if isinstance(val, float) else val
            writer.writerow(row)

    print(f"Results: {path}", flush=True)


def write_predictions_csv(all_scored, path):
    """Write surrogate predictions for all combos."""
    fieldnames = ["combo", "predicted_score"] + [f"predicted_{k}" for k in SCORE_KEYS]
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for combo, predicted, score in all_scored:
            row = {
                "combo": combo_label(combo),
                "predicted_score": f"{score:.4f}",
            }
            for key in SCORE_KEYS:
                row[f"predicted_{key}"] = f"{predicted[key]:.6g}" if predicted.get(key) is not None else ""
            writer.writerow(row)


def print_synergy_report(combo_results):
    """Print synergy analysis."""
    if not combo_results:
        return

    print("\n" + "=" * 80)
    print("SYNERGY ANALYSIS")
    print("=" * 80)
    print(f"{'Combo':<30} {'Synergy':>8}  {'Closure':>8}  {'Scar':>8}  {'Infl':>8}  {'Verdict'}")
    print("-" * 80)

    for r in sorted(combo_results, key=lambda x: -x.get("synergy_mean", 0)):
        label = r["config"]
        syn = r.get("synergy_mean", 0)
        if syn > 0.10:
            verdict = "SYNERGISTIC"
        elif syn < -0.10:
            verdict = "ANTAGONISTIC"
        else:
            verdict = "additive"

        closure = r["outcomes"].get("wound_closure_pct", 0)
        scar = r["outcomes"].get("scar_magnitude", 0)
        infl = r["outcomes"].get("peak_inflammation", 0)
        print(f"{label:<30} {syn:>+8.1%}  {closure:>7.1f}%  {scar:>8.3f}  {infl:>8.4f}  {verdict}")

    print()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Adaptive surrogate-guided treatment combination search")
    parser.add_argument("--treatments", default="all",
                        help="Comma-separated treatments or 'all' (default: all)")
    parser.add_argument("--min-runs", type=int, default=3,
                        help="Minimum runs per config (default: 3)")
    parser.add_argument("--max-runs", type=int, default=10,
                        help="Maximum runs per config (default: 10)")
    parser.add_argument("--cv-threshold", type=float, default=0.05,
                        help="CV threshold for early stopping (default: 0.05)")
    parser.add_argument("--top-k", type=int, default=20,
                        help="Number of top combos to run (default: 20)")
    args = parser.parse_args()

    all_treatments = available_treatments()
    # Exclude 'combination' preset (hand-tuned, not a single treatment)
    all_treatments = [t for t in all_treatments if t != "combination"]

    if args.treatments == "all":
        treatments = all_treatments
    else:
        treatments = [t.strip() for t in args.treatments.split(",")]
        treatments = [t for t in treatments if t != "combination"]
        for t in treatments:
            if t not in all_treatments:
                print(f"ERROR: unknown treatment '{t}'. Available: {all_treatments}")
                sys.exit(1)

    # Generate all possible combos (size >= 2)
    all_combos = []
    for r in range(2, len(treatments) + 1):
        all_combos += list(itertools.combinations(treatments, r))

    n_total_combos = len(all_combos)
    brute_force_runs = (1 + len(treatments) + n_total_combos) * args.max_runs
    adaptive_est = (1 + len(treatments)) * args.min_runs + args.top_k * args.min_runs

    print("=" * 60)
    print("  Adaptive Treatment Combination Search")
    print("=" * 60)
    print(f"  Treatments:      {len(treatments)} ({', '.join(treatments)})")
    print(f"  Total combos:    {n_total_combos}")
    print(f"  Top-K to run:    {args.top_k}")
    print(f"  Runs/config:     {args.min_runs}-{args.max_runs} (adaptive)")
    print(f"  CV threshold:    {args.cv_threshold}")
    print(f"  Brute force:     ~{brute_force_runs} runs")
    print(f"  Adaptive est:    ~{adaptive_est} runs")
    print("=" * 60)
    print()

    os.makedirs(OUTPUT_BASE, exist_ok=True)
    all_results = []
    t_global = time.time()

    # ------------------------------------------------------------------
    # Phase 1: Baseline
    # ------------------------------------------------------------------
    print("PHASE 1: Baseline (diabetic, no treatment)")
    print("-" * 40)
    baseline_mean, baseline_std, n_base, _ = adaptive_consensus(
        "baseline", None, args.min_runs, args.max_runs, args.cv_threshold)

    if not baseline_mean:
        print("ERROR: baseline failed completely")
        sys.exit(1)

    all_results.append({
        "config": "baseline",
        "type": "baseline",
        "n_runs": n_base,
        "outcomes": baseline_mean,
    })
    closure = baseline_mean.get("wound_closure_pct", 0)
    print(f"  Baseline closure: {closure:.1f}% ({n_base} runs)\n")

    # ------------------------------------------------------------------
    # Phase 2: Singles
    # ------------------------------------------------------------------
    print("PHASE 2: Single treatments")
    print("-" * 40)
    singles = {}
    for i, treatment in enumerate(treatments):
        print(f"\n[{i+1}/{len(treatments)}] {treatment.upper()}")
        mean_out, std_out, n_runs, _ = adaptive_consensus(
            treatment, treatment, args.min_runs, args.max_runs, args.cv_threshold)

        if not mean_out:
            print(f"  {treatment} failed")
            continue

        singles[treatment] = mean_out
        all_results.append({
            "config": treatment,
            "type": "single",
            "n_runs": n_runs,
            "outcomes": mean_out,
        })

        closure = mean_out.get("wound_closure_pct", 0)
        base_closure = baseline_mean.get("wound_closure_pct", 0)
        delta = closure - base_closure
        print(f"  {treatment}: closure {closure:.1f}% ({delta:+.1f}% vs baseline)")

    if not singles:
        print("ERROR: no singles succeeded")
        sys.exit(1)

    # ------------------------------------------------------------------
    # Phase 3: Surrogate model
    # ------------------------------------------------------------------
    print(f"\n\nPHASE 3: Surrogate model")
    print("-" * 40)
    effects = fit_surrogate(baseline_mean, singles)

    for name, delta in effects.items():
        parts = []
        for key in SCORE_KEYS:
            d = delta.get(key, 0)
            parts.append(f"{key}={d:+.3g}")
        print(f"  {name}: {', '.join(parts)}")

    # Predict all combos
    all_scored = []
    for combo in all_combos:
        predicted = predict_combo(baseline_mean, effects, combo)
        sc = composite_score(predicted, baseline_mean)
        all_scored.append((combo, predicted, sc))
    all_scored.sort(key=lambda x: x[2])

    # Write predictions
    pred_path = os.path.join(OUTPUT_BASE, "surrogate_predictions.csv")
    write_predictions_csv(all_scored, pred_path)
    print(f"\n  Predictions for {len(all_scored)} combos: {pred_path}")

    # ------------------------------------------------------------------
    # Phase 4: Selection
    # ------------------------------------------------------------------
    print(f"\n\nPHASE 4: Select top-{args.top_k} combos")
    print("-" * 40)
    top_k = min(args.top_k, len(all_scored))
    selected = all_scored[:top_k]

    print(f"  Top {top_k} predicted combos:")
    for i, (combo, predicted, score) in enumerate(selected):
        print(f"    {i+1}. {combo_label(combo)} (score: {score:.3f})")

    # ------------------------------------------------------------------
    # Phase 5: Run selected combos
    # ------------------------------------------------------------------
    print(f"\n\nPHASE 5: Run {top_k} selected combos")
    print("-" * 40)
    combo_results = []

    for i, (combo, predicted, pred_score) in enumerate(selected):
        label = combo_label(combo)
        print(f"\n[{i+1}/{top_k}] COMBO: {label}")

        mean_out, std_out, n_runs, _ = adaptive_consensus(
            label, list(combo), args.min_runs, args.max_runs, args.cv_threshold)

        if not mean_out:
            print(f"  {label} failed")
            continue

        # Compute synergy
        synergy = compute_synergy(mean_out, predicted, baseline_mean)

        result = {
            "config": label,
            "type": "combo",
            "n_runs": n_runs,
            "outcomes": mean_out,
            "synergy_mean": synergy["mean"],
        }
        for key in SCORE_KEYS:
            result[f"synergy_{key}"] = synergy[key]
            result[f"predicted_{key}"] = predicted[key]

        combo_results.append(result)
        all_results.append(result)

        closure = mean_out.get("wound_closure_pct", 0)
        syn_str = f"synergy: {synergy['mean']:+.1%}"
        print(f"  {label}: closure {closure:.1f}% ({syn_str})")

    # ------------------------------------------------------------------
    # Phase 6: Report
    # ------------------------------------------------------------------
    print(f"\n\nPHASE 6: Results")
    print("-" * 40)

    # Synergy report
    print_synergy_report(combo_results)

    # Best configs overall
    scored_results = []
    for r in all_results:
        if r.get("outcomes"):
            sc = composite_score(r["outcomes"], baseline_mean)
            scored_results.append((r["config"], sc, r["outcomes"]))
    scored_results.sort(key=lambda x: x[1])

    print("\nTOP 10 CONFIGS (by composite score, lower = better):")
    print(f"{'Rank':<6} {'Config':<30} {'Score':>7} {'Closure':>8} {'Scar':>8} {'T50':>8}")
    print("-" * 70)
    for i, (name, sc, out) in enumerate(scored_results[:10]):
        closure = out.get("wound_closure_pct", 0)
        scar = out.get("scar_magnitude", 0)
        t50 = out.get("time_to_50pct_days")
        t50_str = f"{t50:.1f}d" if t50 is not None else "N/A"
        print(f"{i+1:<6} {name:<30} {sc:>7.3f} {closure:>7.1f}% {scar:>8.3f} {t50_str:>8}")

    # Write CSV
    results_path = os.path.join(OUTPUT_BASE, "adaptive_study.csv")
    write_results_csv(all_results, results_path)

    # Timing
    elapsed = time.time() - t_global
    total_runs = sum(r.get("n_runs", 0) for r in all_results)
    print(f"\nTotal: {total_runs} runs in {elapsed/60:.1f} minutes")
    print(f"Output: {OUTPUT_BASE}/")


if __name__ == "__main__":
    main()
