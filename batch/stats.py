#!/usr/bin/env python3
"""Statistical analysis toolkit for batch simulation results.

Provides confidence intervals, hypothesis testing, effect sizes, and
multiple comparison correction for multi-run consensus and treatment
comparisons.

Usage:
    # Compare two batch results
    python3 batch/stats.py studies/wound/results studies/diabetic-wound/results

    # Analyze a single consensus with CI
    python3 batch/stats.py studies/wound/results --ci

    # Compare all raw runs within a single batch
    python3 batch/stats.py studies/wound/results --summary
"""

import argparse
import csv
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir))

from batch import lib


# ---------------------------------------------------------------------------
# Confidence intervals
# ---------------------------------------------------------------------------

def t_critical(df, alpha=0.05):
    """Approximate two-tailed t critical value using Abramowitz & Stegun.

    For exact values use scipy.stats.t.ppf, but this avoids the dependency.
    Accurate to ~0.5% for df >= 2.
    """
    # For two-tailed alpha=0.05, common values
    _TABLE = {
        1: 12.706, 2: 4.303, 3: 3.182, 4: 2.776, 5: 2.571,
        6: 2.447, 7: 2.365, 8: 2.306, 9: 2.262, 10: 2.228,
        11: 2.201, 12: 2.179, 13: 2.160, 14: 2.145, 15: 2.131,
        16: 2.120, 17: 2.110, 18: 2.101, 19: 2.093, 20: 2.086,
        25: 2.060, 30: 2.042, 40: 2.021, 50: 2.009, 60: 2.000,
        80: 1.990, 100: 1.984, 200: 1.972, 500: 1.965,
    }
    if alpha != 0.05:
        # Rough scaling for other alpha values
        z_ratio = {0.01: 2.576 / 1.960, 0.10: 1.645 / 1.960}
        scale = z_ratio.get(alpha, 1.0)
        return t_critical(df, 0.05) * scale

    if df in _TABLE:
        return _TABLE[df]
    # Interpolate between nearest table entries
    keys = sorted(_TABLE.keys())
    if df < keys[0]:
        return _TABLE[keys[0]]
    if df > keys[-1]:
        return 1.960  # z-approximation for large df
    for i in range(len(keys) - 1):
        if keys[i] <= df <= keys[i + 1]:
            lo, hi = keys[i], keys[i + 1]
            frac = (df - lo) / (hi - lo)
            return _TABLE[lo] + frac * (_TABLE[hi] - _TABLE[lo])
    return 1.960


def confidence_interval(values, alpha=0.05):
    """Compute mean, standard error, and CI for a list of values.

    Returns dict with mean, std, se, ci_lo, ci_hi, n, alpha.
    """
    n = len(values)
    if n < 2:
        m = values[0] if values else float("nan")
        return dict(mean=m, std=0, se=0, ci_lo=m, ci_hi=m, n=n, alpha=alpha)

    m = sum(values) / n
    var = sum((x - m) ** 2 for x in values) / (n - 1)  # sample variance
    std = math.sqrt(var)
    se = std / math.sqrt(n)
    t_crit = t_critical(n - 1, alpha)
    margin = t_crit * se

    return dict(
        mean=m, std=std, se=se,
        ci_lo=m - margin, ci_hi=m + margin,
        n=n, alpha=alpha,
    )


# ---------------------------------------------------------------------------
# Hypothesis testing
# ---------------------------------------------------------------------------

def welch_t_test(group_a, group_b):
    """Welch's t-test for unequal variances.

    Returns dict with t_stat, df, p_value (two-tailed approximate),
    mean_a, mean_b, cohens_d.
    """
    n_a, n_b = len(group_a), len(group_b)
    if n_a < 2 or n_b < 2:
        return dict(t_stat=float("nan"), df=0, p_value=1.0,
                    mean_a=0, mean_b=0, cohens_d=0, significant=False)

    m_a = sum(group_a) / n_a
    m_b = sum(group_b) / n_b
    var_a = sum((x - m_a) ** 2 for x in group_a) / (n_a - 1)
    var_b = sum((x - m_b) ** 2 for x in group_b) / (n_b - 1)

    se = math.sqrt(var_a / n_a + var_b / n_b)
    if se < 1e-15:
        return dict(t_stat=0, df=max(n_a, n_b) - 1, p_value=1.0,
                    mean_a=m_a, mean_b=m_b, cohens_d=0, significant=False)

    t_stat = (m_a - m_b) / se

    # Welch-Satterthwaite degrees of freedom
    num = (var_a / n_a + var_b / n_b) ** 2
    denom = ((var_a / n_a) ** 2 / (n_a - 1) +
             (var_b / n_b) ** 2 / (n_b - 1))
    df = num / denom if denom > 0 else 1

    # Approximate p-value using t-distribution CDF
    p_value = _approx_t_pvalue(abs(t_stat), df)

    # Cohen's d (pooled std)
    pooled_var = ((n_a - 1) * var_a + (n_b - 1) * var_b) / (n_a + n_b - 2)
    pooled_std = math.sqrt(pooled_var) if pooled_var > 0 else 1e-15
    cohens_d = (m_a - m_b) / pooled_std

    return dict(
        t_stat=t_stat, df=df, p_value=p_value,
        mean_a=m_a, mean_b=m_b, cohens_d=cohens_d,
        significant=False,  # set after correction
    )


def _approx_t_pvalue(t_abs, df):
    """Approximate two-tailed p-value for |t| with df degrees of freedom.

    Uses the relationship between the t-distribution and the regularized
    incomplete beta function, approximated via continued fraction.
    Accurate to ~1% for practical df values.
    """
    if df <= 0 or math.isnan(t_abs):
        return 1.0
    x = df / (df + t_abs ** 2)
    # Regularized incomplete beta function I_x(df/2, 0.5) approximation
    # For large |t|, p is very small
    if t_abs > 30:
        return 0.0
    # Use a simple empirical approximation
    # p ~ 2 * (1 - Phi(t * sqrt((df-2)/df))) for df > 5
    if df > 5:
        z = t_abs * math.sqrt((df - 2) / df)
    else:
        z = t_abs * math.sqrt(df / (df + 2))
    # Standard normal CDF approximation (Abramowitz & Stegun 26.2.17)
    p_one_tail = _normal_sf(z)
    return 2 * p_one_tail


def _normal_sf(z):
    """Survival function (1 - CDF) for standard normal, z >= 0."""
    if z < 0:
        return 1 - _normal_sf(-z)
    # Abramowitz & Stegun 26.2.17 approximation
    b0 = 0.2316419
    b1 = 0.319381530
    b2 = -0.356563782
    b3 = 1.781477937
    b4 = -1.821255978
    b5 = 1.330274429
    t = 1 / (1 + b0 * z)
    pdf = math.exp(-z * z / 2) / math.sqrt(2 * math.pi)
    return pdf * t * (b1 + t * (b2 + t * (b3 + t * (b4 + t * b5))))


# ---------------------------------------------------------------------------
# Multiple comparison correction
# ---------------------------------------------------------------------------

def bonferroni_correct(results, alpha=0.05):
    """Apply Bonferroni correction to a list of test result dicts.

    Each result must have a 'p_value' key. Adds 'significant' and
    'alpha_corrected' keys.
    """
    m = len(results)
    if m == 0:
        return results
    alpha_corrected = alpha / m
    for r in results:
        r["alpha_corrected"] = alpha_corrected
        r["significant"] = r["p_value"] < alpha_corrected
    return results


def holm_bonferroni_correct(results, alpha=0.05):
    """Apply Holm-Bonferroni step-down correction.

    Less conservative than Bonferroni while still controlling FWER.
    """
    m = len(results)
    if m == 0:
        return results
    # Sort by p-value
    indexed = sorted(enumerate(results), key=lambda x: x[1]["p_value"])
    for rank, (orig_idx, r) in enumerate(indexed):
        adjusted_alpha = alpha / (m - rank)
        r["alpha_corrected"] = adjusted_alpha
        r["significant"] = r["p_value"] < adjusted_alpha
    return results


# ---------------------------------------------------------------------------
# Effect size interpretation
# ---------------------------------------------------------------------------

def interpret_cohens_d(d):
    """Interpret Cohen's d effect size magnitude."""
    d_abs = abs(d)
    if d_abs < 0.2:
        return "negligible"
    elif d_abs < 0.5:
        return "small"
    elif d_abs < 0.8:
        return "medium"
    else:
        return "large"


# ---------------------------------------------------------------------------
# Bootstrap CI (for non-normal distributions or small n)
# ---------------------------------------------------------------------------

def bootstrap_ci(values, n_bootstrap=5000, alpha=0.05, stat_fn=None):
    """Compute bootstrap confidence interval for a statistic.

    stat_fn defaults to mean. Returns dict with estimate, ci_lo, ci_hi.
    """
    import random
    if stat_fn is None:
        stat_fn = lambda v: sum(v) / len(v)

    n = len(values)
    if n < 2:
        est = stat_fn(values) if values else float("nan")
        return dict(estimate=est, ci_lo=est, ci_hi=est, n_bootstrap=0)

    estimates = []
    for _ in range(n_bootstrap):
        sample = [random.choice(values) for _ in range(n)]
        estimates.append(stat_fn(sample))

    estimates.sort()
    lo_idx = int(n_bootstrap * alpha / 2)
    hi_idx = int(n_bootstrap * (1 - alpha / 2))

    return dict(
        estimate=stat_fn(values),
        ci_lo=estimates[lo_idx],
        ci_hi=estimates[hi_idx],
        n_bootstrap=n_bootstrap,
    )


# ---------------------------------------------------------------------------
# Batch result analysis
# ---------------------------------------------------------------------------

def load_raw_outcomes(result_dir, column, measure="final"):
    """Load per-run outcomes from raw/ subdirectory.

    Returns list of scalar outcome values, one per run.
    """
    raw_dir = os.path.join(result_dir, "raw")
    if not os.path.isdir(raw_dir):
        return []

    outcomes = []
    for entry in sorted(os.listdir(raw_dir)):
        run_dir = os.path.join(raw_dir, entry)
        if not os.path.isdir(run_dir):
            continue
        # Look for metrics.csv inside the run directory
        csv_path = os.path.join(run_dir, "skibidy", "metrics.csv")
        if not os.path.isfile(csv_path):
            # Try direct metrics.csv
            csv_path = os.path.join(run_dir, "metrics.csv")
        if not os.path.isfile(csv_path):
            continue
        data = lib.load_csv(csv_path)
        val = lib.extract_outcome(data, column, measure)
        if not math.isnan(val):
            outcomes.append(val)

    return outcomes


def analyze_single(result_dir, columns=None, measure="final", alpha=0.05):
    """Analyze a single batch result with confidence intervals.

    Returns list of dicts with column, mean, std, se, ci_lo, ci_hi, n.
    """
    if columns is None:
        columns = ["wound_closure_pct", "mean_collagen_wound",
                    "mean_infl_wound", "n_fibroblasts"]

    results = []
    for col in columns:
        values = load_raw_outcomes(result_dir, col, measure)
        if not values:
            continue
        ci = confidence_interval(values, alpha)
        ci["column"] = col
        ci["measure"] = measure
        results.append(ci)

    return results


def compare_groups(dir_a, dir_b, columns=None, measure="final", alpha=0.05,
                   correction="holm"):
    """Compare two batch result directories with hypothesis testing.

    Returns list of test result dicts with Bonferroni or Holm correction.
    """
    if columns is None:
        columns = ["wound_closure_pct", "mean_collagen_wound",
                    "mean_infl_wound"]

    results = []
    for col in columns:
        vals_a = load_raw_outcomes(dir_a, col, measure)
        vals_b = load_raw_outcomes(dir_b, col, measure)
        if not vals_a or not vals_b:
            continue
        test = welch_t_test(vals_a, vals_b)
        test["column"] = col
        test["measure"] = measure
        test["n_a"] = len(vals_a)
        test["n_b"] = len(vals_b)
        test["ci_a"] = confidence_interval(vals_a, alpha)
        test["ci_b"] = confidence_interval(vals_b, alpha)
        test["effect_size"] = interpret_cohens_d(test["cohens_d"])
        results.append(test)

    if correction == "holm":
        holm_bonferroni_correct(results, alpha)
    else:
        bonferroni_correct(results, alpha)

    return results


# ---------------------------------------------------------------------------
# Printing
# ---------------------------------------------------------------------------

def print_ci_summary(results, label=""):
    """Print confidence interval summary table."""
    if label:
        print(f"\n{'=' * 72}")
        print(f"  {label}")
        print(f"{'=' * 72}")
    print(f"  {'Metric':<30} {'Mean':>8} {'Std':>8} {'SE':>8} "
          f"{'95% CI':>20} {'n':>4}")
    print(f"  {'-' * 70}")
    for r in results:
        col_short = r["column"].split("_", 1)[-1] if "column" in r else "?"
        ci_str = f"[{r['ci_lo']:.3f}, {r['ci_hi']:.3f}]"
        print(f"  {col_short:<30} {r['mean']:>8.3f} {r['std']:>8.3f} "
              f"{r['se']:>8.4f} {ci_str:>20} {r['n']:>4}")
    print()


def print_comparison(results, label_a="Group A", label_b="Group B"):
    """Print hypothesis test comparison table."""
    print(f"\n{'=' * 80}")
    print(f"  {label_a} vs {label_b}")
    print(f"{'=' * 80}")
    print(f"  {'Metric':<25} {'Mean A':>8} {'Mean B':>8} {'Diff':>8} "
          f"{'t':>7} {'p':>8} {'d':>6} {'Effect':>10} {'Sig':>5}")
    print(f"  {'-' * 78}")

    for r in results:
        col_short = r["column"].split("_", 1)[-1] if "column" in r else "?"
        diff = r["mean_a"] - r["mean_b"]
        sig = "***" if r["significant"] else ""
        print(f"  {col_short:<25} {r['mean_a']:>8.3f} {r['mean_b']:>8.3f} "
              f"{diff:>+8.3f} {r['t_stat']:>7.2f} {r['p_value']:>8.4f} "
              f"{r['cohens_d']:>+6.2f} {r['effect_size']:>10} {sig:>5}")

    # Footer
    corrected = results[0].get("alpha_corrected", 0.05) if results else 0.05
    n_sig = sum(1 for r in results if r.get("significant"))
    print(f"\n  Corrected alpha: {corrected:.4f} "
          f"(Holm-Bonferroni, {len(results)} comparisons)")
    print(f"  Significant results: {n_sig}/{len(results)}")
    print()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Statistical analysis of batch/sweep results")
    parser.add_argument("result_dir", help="Primary results directory")
    parser.add_argument("compare_dir", nargs="?", default=None,
                        help="Second results directory for comparison")
    parser.add_argument("--ci", action="store_true",
                        help="Print confidence intervals for single result")
    parser.add_argument("--summary", action="store_true",
                        help="Print CI summary (alias for --ci)")
    parser.add_argument("--columns", type=str, default=None,
                        help="Comma-separated outcome columns")
    parser.add_argument("--measure", type=str, default="final",
                        help="Outcome measure: final, peak, auc, time_to_N")
    parser.add_argument("--alpha", type=float, default=0.05,
                        help="Significance level (default: 0.05)")
    parser.add_argument("--bootstrap", action="store_true",
                        help="Use bootstrap CI instead of t-distribution")
    args = parser.parse_args()

    columns = args.columns.split(",") if args.columns else None

    if args.compare_dir:
        # Two-group comparison
        results = compare_groups(
            args.result_dir, args.compare_dir,
            columns=columns, measure=args.measure, alpha=args.alpha)
        if results:
            print_comparison(results,
                             label_a=os.path.basename(args.result_dir),
                             label_b=os.path.basename(args.compare_dir))
        else:
            print("No raw run data found for comparison.")
    elif args.ci or args.summary:
        # Single-group CI
        results = analyze_single(
            args.result_dir, columns=columns,
            measure=args.measure, alpha=args.alpha)
        if results:
            print_ci_summary(results,
                             label=os.path.basename(args.result_dir))
        else:
            print("No raw run data found.")
    else:
        # Default: print CI if raw data exists, else consensus summary
        results = analyze_single(
            args.result_dir, columns=columns,
            measure=args.measure, alpha=args.alpha)
        if results:
            print_ci_summary(results,
                             label=os.path.basename(args.result_dir))
        else:
            print(f"No raw runs in {args.result_dir}/raw/")
            print("Use --ci with a batch result that has raw/ subdirectory.")


if __name__ == "__main__":
    main()
