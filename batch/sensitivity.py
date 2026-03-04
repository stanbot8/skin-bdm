#!/usr/bin/env python3
"""Morris method sensitivity analysis.

Computes elementary effects (EE) for each parameter to rank their influence
on simulation outputs. The Morris method is a global screening technique
that identifies which parameters have: (a) negligible effects, (b) linear
and additive effects, or (c) nonlinear or interaction effects.

For each parameter, reports:
  mu*   = mean of |EE|  (overall importance)
  sigma = std of EE     (nonlinearity/interactions)

Parameters with high mu* and low sigma have linear, additive effects.
Parameters with high mu* and high sigma are involved in interactions
or have nonlinear effects.

Usage:
    python3 batch/sensitivity.py batch/configs/sensitivity_wound.toml
    python3 batch/sensitivity.py batch/configs/sensitivity_wound.toml -r 10
    python3 batch/sensitivity.py batch/configs/sensitivity_wound.toml --levels 6
"""

import argparse
import csv
import math
import os
import shutil
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir))

from datetime import datetime

from batch import lib

try:
    import tomllib
except ImportError:
    try:
        import tomli as tomllib
    except ImportError:
        tomllib = None


def parse_toml(path):
    if tomllib:
        with open(path, "rb") as f:
            return tomllib.load(f)
    from batch.sweep import _parse_toml_fallback
    return _parse_toml_fallback(path)


def load_config(path):
    """Load sensitivity analysis config."""
    cfg = parse_toml(path)
    sa = cfg.get("sensitivity", {})
    outcomes = cfg.get("outcomes", {})
    params = cfg.get("params", cfg.get("sensitivity", {}).get("params", []))
    if isinstance(params, dict):
        params = params.get("params", [])

    # Handle [[params]] array of tables
    if "params" in cfg and isinstance(cfg["params"], list):
        params = cfg["params"]
    elif "sensitivity" in cfg and "params" not in cfg:
        params = sa.get("params", [])

    return dict(
        name=sa.get("name", os.path.splitext(os.path.basename(path))[0]),
        study=sa.get("study"),
        skin=sa.get("skin"),
        treatment=sa.get("treatment"),
        trajectories=sa.get("trajectories", 10),
        levels=sa.get("levels", 4),
        replicates=sa.get("replicates", 1),
        params=params,
        primary=outcomes.get("primary", "wound_closure_pct"),
        secondary=outcomes.get("secondary", []),
        measure=outcomes.get("measure", "final"),
        source_path=path,
    )


# ---------------------------------------------------------------------------
# Morris trajectory generation
# ---------------------------------------------------------------------------

def generate_morris_trajectories(params, r, p):
    """Generate r Morris trajectories for k parameters with p levels.

    Each trajectory has k+1 points (base + one perturbation per param).
    Returns list of trajectories, each a list of k+1 dicts mapping
    param name to its normalized level (0 to p-1).

    Uses optimized random sampling with distance maximization.
    """
    import random
    k = len(params)
    delta = p // 2  # standard Morris delta = p/(2*(p-1)) but we use grid levels

    trajectories = []
    for _ in range(r * 10):  # generate candidates, keep best r
        # Random base point on the grid
        base = [random.randint(0, p - 1) for _ in range(k)]
        # Random parameter order
        order = list(range(k))
        random.shuffle(order)

        traj = [list(base)]
        current = list(base)
        for idx in order:
            # Perturb parameter idx by +delta or -delta
            if current[idx] + delta <= p - 1:
                current[idx] += delta
            elif current[idx] - delta >= 0:
                current[idx] -= delta
            else:
                current[idx] = min(current[idx] + delta, p - 1)
            traj.append(list(current))

        trajectories.append((traj, order))

    # Select r trajectories maximizing minimum pairwise distance
    if len(trajectories) <= r:
        selected = trajectories
    else:
        selected = _select_spread_trajectories(trajectories, r)

    return selected


def _select_spread_trajectories(candidates, r):
    """Greedy selection of r trajectories maximizing spread."""
    import random
    # Start with a random trajectory
    selected = [random.choice(candidates)]
    remaining = [t for t in candidates if t is not selected[0]]

    while len(selected) < r and remaining:
        best_dist = -1
        best_idx = 0
        for i, cand in enumerate(remaining):
            min_dist = min(_traj_distance(cand, s) for s in selected)
            if min_dist > best_dist:
                best_dist = min_dist
                best_idx = i
        selected.append(remaining.pop(best_idx))

    return selected


def _traj_distance(t1, t2):
    """Sum of squared differences between all point pairs."""
    total = 0
    for p1 in t1[0]:
        for p2 in t2[0]:
            total += sum((a - b) ** 2 for a, b in zip(p1, p2))
    return total


def levels_to_values(param_cfg, level, p):
    """Convert a grid level (0..p-1) to actual parameter value."""
    lo = param_cfg["min"]
    hi = param_cfg["max"]
    if param_cfg.get("log", False):
        # Log-uniform spacing
        log_lo = math.log(lo) if lo > 0 else math.log(1e-10)
        log_hi = math.log(hi)
        log_val = log_lo + (log_hi - log_lo) * level / max(p - 1, 1)
        return math.exp(log_val)
    else:
        return lo + (hi - lo) * level / max(p - 1, 1)


# ---------------------------------------------------------------------------
# Elementary effect computation
# ---------------------------------------------------------------------------

def compute_elementary_effects(trajectories, results, params, p):
    """Compute elementary effects from trajectory results.

    results[traj_idx][point_idx] = {outcome: value}

    Returns dict: param_name -> {outcome -> [list of EE values]}
    """
    k = len(params)
    ee = {pc["param"]: {out: [] for out in results[0][0]} for pc in params}

    for traj_idx, (traj, order) in enumerate(trajectories):
        for step, param_idx in enumerate(order):
            pc = params[param_idx]
            pname = pc["param"]

            # Points before and after perturbation
            pt_before = traj[step]
            pt_after = traj[step + 1]

            # Level change
            delta_level = pt_after[param_idx] - pt_before[param_idx]
            if delta_level == 0:
                continue

            # Value change (for normalization)
            val_before = levels_to_values(pc, pt_before[param_idx], p)
            val_after = levels_to_values(pc, pt_after[param_idx], p)
            delta_val = val_after - val_before

            # Output change
            r_before = results[traj_idx][step]
            r_after = results[traj_idx][step + 1]

            for outcome in r_before:
                y_before = r_before[outcome]
                y_after = r_after[outcome]
                # Skip if either point has NaN (failed simulation)
                if math.isnan(y_before) or math.isnan(y_after):
                    continue
                dy = y_after - y_before
                # Normalized elementary effect: (dy/y) / (dx/x)
                # Use absolute normalization to avoid division by zero
                y_ref = y_before
                x_ref = val_before
                if abs(x_ref) > 1e-10 and abs(y_ref) > 1e-10:
                    ee_val = (dy / y_ref) / (delta_val / x_ref)
                elif abs(delta_val) > 1e-10:
                    # Fallback: absolute EE
                    ee_val = dy / delta_val
                else:
                    ee_val = 0.0
                ee[pname][outcome].append(ee_val)

    return ee


def summarize_effects(ee):
    """Compute mu* (mean |EE|) and sigma (std EE) for each param/outcome.

    Returns dict: param_name -> {outcome -> {mu_star, sigma, mu, n}}
    """
    summary = {}
    for pname, outcomes in ee.items():
        summary[pname] = {}
        for outcome, values in outcomes.items():
            if not values:
                summary[pname][outcome] = {
                    "mu_star": 0, "sigma": 0, "mu": 0, "n": 0}
                continue
            n = len(values)
            mu = sum(values) / n
            mu_star = sum(abs(v) for v in values) / n
            variance = sum((v - mu) ** 2 for v in values) / max(n - 1, 1)
            sigma = math.sqrt(variance)
            summary[pname][outcome] = {
                "mu_star": mu_star, "sigma": sigma, "mu": mu, "n": n}
    return summary


# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------

def run_sensitivity(cfg, override_r=None, override_levels=None):
    """Execute Morris sensitivity analysis."""
    r = override_r or cfg["trajectories"]
    p = override_levels or cfg["levels"]
    params = cfg["params"]
    k = len(params)
    replicates = cfg["replicates"]
    all_outcomes = [cfg["primary"]] + cfg["secondary"]
    measure = cfg["measure"]

    total_sims = r * (k + 1) * replicates
    print(f"=== Morris Sensitivity Analysis: {cfg['name']} ===")
    print(f"  Parameters: {k}")
    print(f"  Trajectories: {r}")
    print(f"  Levels: {p}")
    print(f"  Replicates per point: {replicates}")
    print(f"  Total simulations: {total_sims}")
    est_time = total_sims * 36  # ~36s per sim
    print(f"  Estimated time: {est_time // 3600}h {(est_time % 3600) // 60}m")
    print()

    # Generate trajectories
    trajectories = generate_morris_trajectories(params, r, p)
    print(f"  Generated {len(trajectories)} trajectories")

    # Output directory
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    result_dir = os.path.join(
        os.path.dirname(__file__), "results",
        f"sensitivity_{cfg['name']}_{timestamp}")
    raw_dir = os.path.join(result_dir, "raw")
    os.makedirs(raw_dir, exist_ok=True)
    shutil.copy2(cfg["source_path"], os.path.join(result_dir, "config.toml"))

    # Build once
    lib.build_if_needed()

    # Run all trajectory points
    results = []  # results[traj][point] = {outcome: value}
    sim_count = 0
    t_start = time.time()

    for traj_idx, (traj, order) in enumerate(trajectories):
        traj_results = []
        print(f"\n--- Trajectory {traj_idx + 1}/{len(trajectories)} ---")

        for pt_idx, point in enumerate(traj):
            sim_count += 1
            # Convert levels to actual parameter values
            param_values = {}
            for i, pc in enumerate(params):
                param_values[pc["param"]] = levels_to_values(pc, point[i], p)

            label = f"[{sim_count}/{total_sims}]"
            short = ", ".join(
                f"{pc['param'].split('.')[-1]}={param_values[pc['param']]:.3g}"
                for pc in params[:3])
            if len(params) > 3:
                short += f" +{len(params) - 3}"
            print(f"  {label} {short}", end=" ", flush=True)

            # Run with replicates and average
            point_outcomes = {out: [] for out in all_outcomes}
            all_ok = True

            for rep in range(replicates):
                lib.merge_config()
                if cfg["skin"]:
                    lib.apply_profile(cfg["skin"])
                if cfg["study"]:
                    lib.apply_study(cfg["study"])
                if cfg.get("treatment"):
                    lib.apply_treatment(cfg["treatment"], cfg["study"])

                for param_path, value in param_values.items():
                    lib.override_param(param_path, value)

                run_dir = os.path.join(
                    raw_dir, f"t{traj_idx:03d}_p{pt_idx:02d}_r{rep:02d}")
                ok, elapsed = lib.run_simulation(output_path=run_dir)

                if not ok:
                    print(f"FAIL", end=" ")
                    all_ok = False
                    continue

                csv_path = lib.get_metrics_path()
                if not os.path.isfile(csv_path):
                    all_ok = False
                    continue

                data = lib.load_csv(csv_path)
                for out in all_outcomes:
                    val = lib.extract_outcome(data, out, measure)
                    if not math.isnan(val):
                        point_outcomes[out].append(val)

            # Average replicates
            avg = {}
            for out in all_outcomes:
                vals = point_outcomes[out]
                avg[out] = sum(vals) / len(vals) if vals else float("nan")

            traj_results.append(avg)
            if all_ok:
                print(f"OK ({time.time() - t_start:.0f}s total)")
            else:
                print(f"partial")

        results.append(traj_results)

    elapsed_total = time.time() - t_start

    # Compute elementary effects
    ee = compute_elementary_effects(trajectories, results, params, p)
    summary = summarize_effects(ee)

    # Print results table
    print(f"\n{'=' * 72}")
    print(f"  Morris Sensitivity Analysis Results ({cfg['name']})")
    print(f"  {sim_count} simulations in {elapsed_total:.0f}s")
    print(f"{'=' * 72}")

    for outcome in all_outcomes:
        print(f"\n  Observable: {outcome} ({measure})")
        print(f"  {'Parameter':<40} {'mu*':>8} {'sigma':>8} {'mu':>8} {'n':>4}")
        print(f"  {'-' * 68}")

        # Sort by mu* descending
        ranked = sorted(
            summary.items(),
            key=lambda x: x[1].get(outcome, {}).get("mu_star", 0),
            reverse=True)

        for pname, outcomes_data in ranked:
            s = outcomes_data.get(outcome, {})
            mu_star = s.get("mu_star", 0)
            sigma = s.get("sigma", 0)
            mu = s.get("mu", 0)
            n = s.get("n", 0)
            flag = " ***" if mu_star > 0.5 else " *" if mu_star > 0.1 else ""
            short_name = pname.split(".")[-1]
            print(f"  {short_name:<40} {mu_star:>8.4f} {sigma:>8.4f} {mu:>+8.4f} {n:>4}{flag}")

    # Write summary CSV
    summary_path = os.path.join(result_dir, "sensitivity.csv")
    rows = []
    for pname, outcomes_data in summary.items():
        for outcome, s in outcomes_data.items():
            rows.append({
                "param": pname,
                "outcome": outcome,
                "mu_star": s["mu_star"],
                "sigma": s["sigma"],
                "mu": s["mu"],
                "n": s["n"],
            })
    if rows:
        with open(summary_path, "w", newline="") as f:
            writer = csv.DictWriter(
                f, fieldnames=["param", "outcome", "mu_star", "sigma", "mu", "n"])
            writer.writeheader()
            for row in rows:
                writer.writerow({
                    k: f"{v:.6g}" if isinstance(v, float) else v
                    for k, v in row.items()})
        print(f"\nSaved: {summary_path}")

    # Write raw elementary effects
    ee_path = os.path.join(result_dir, "elementary_effects.csv")
    with open(ee_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["param", "outcome", "trajectory", "ee"])
        for pname, outcomes_data in ee.items():
            for outcome, values in outcomes_data.items():
                for i, val in enumerate(values):
                    writer.writerow([pname, outcome, i, f"{val:.6g}"])
    print(f"Saved: {ee_path}")

    print(f"\nResults: {result_dir}/")
    return result_dir


def main():
    parser = argparse.ArgumentParser(
        description="Morris method sensitivity analysis")
    parser.add_argument("config", help="Path to sensitivity config TOML")
    parser.add_argument("-r", "--trajectories", type=int, default=None,
                        help="Override number of trajectories")
    parser.add_argument("--levels", type=int, default=None,
                        help="Override number of grid levels (default: 4)")
    parser.add_argument("--replicates", type=int, default=None,
                        help="Override replicates per point")
    args = parser.parse_args()

    cfg = load_config(args.config)
    if args.replicates:
        cfg["replicates"] = args.replicates
    run_sensitivity(cfg, override_r=args.trajectories,
                    override_levels=args.levels)


if __name__ == "__main__":
    main()
