"""Cure criteria for chronic conditions.

Defines measurable thresholds for disease resolution (not just symptom
management). Each condition has criteria that must ALL be met for a
simulation to be classified as "cured".

Usage:
    from literature.cure_criteria import evaluate_cure

    result = evaluate_cure(mean_data, condition="diabetic")
    print(result["verdict"])  # "CURED", "PARTIAL", or "UNCURED"
    for c in result["criteria"]:
        print(f"  {c['name']}: {'PASS' if c['met'] else 'FAIL'} ({c['detail']})")
"""

import math


# ---------------------------------------------------------------------------
# Cure thresholds per condition
# ---------------------------------------------------------------------------

CRITERIA = {
    "diabetic": {
        "label": "Diabetic Wound Resolution",
        "checks": [
            {
                "name": "Complete closure",
                "metric": "wound_closure_pct",
                "op": "final_gte",
                "threshold": 95.0,
                "unit": "%",
                "rationale": "Full re-epithelialization, not just 90% functional closure",
            },
            {
                "name": "Inflammation resolved",
                "metric": "mean_infl_wound",
                "op": "tail_mean_lt",
                "threshold": 0.05,
                "tail_frac": 0.25,
                "unit": "normalized",
                "rationale": "Sustained low inflammation for final 25% of sim (no chronic plateau)",
            },
            {
                "name": "MMP normalized",
                "metric": "mean_mmp_wound",
                "op": "tail_mean_lt",
                "threshold": 0.1,
                "tail_frac": 0.25,
                "unit": "normalized",
                "rationale": "MMP activity returns to baseline (not persistent degradation)",
            },
            {
                "name": "Collagen restored",
                "metric": "mean_collagen_wound",
                "op": "final_gte",
                "threshold": 0.5,
                "unit": "normalized",
                "rationale": "Collagen deposition reaches functional dermis level",
            },
            {
                "name": "Neutrophils cleared",
                "metric": "n_neutrophils",
                "op": "tail_mean_lt",
                "threshold": 50.0,
                "tail_frac": 0.15,
                "unit": "cells",
                "rationale": "Neutrophil infiltrate resolved (no persistent recruitment)",
            },
        ],
    },
    "pressure": {
        "label": "Pressure Ulcer Resolution",
        "checks": [
            {
                "name": "Wound closure",
                "metric": "wound_closure_pct",
                "op": "final_gte",
                "threshold": 90.0,
                "unit": "%",
                "rationale": "Stage II/III pressure ulcers close fully with adequate offloading",
            },
            {
                "name": "Inflammation resolved",
                "metric": "mean_infl_wound",
                "op": "tail_mean_lt",
                "threshold": 0.08,
                "tail_frac": 0.2,
                "unit": "normalized",
                "rationale": "Chronic inflammation loop broken (no ischemia-reperfusion cycling)",
            },
            {
                "name": "Macrophage resolution",
                "metric": "n_macrophages",
                "op": "tail_trend_negative",
                "threshold": 0.0,
                "tail_frac": 0.3,
                "unit": "trend",
                "rationale": "Macrophage count declining (M1-to-M2 transition occurring)",
            },
            {
                "name": "Tissue viability",
                "metric": "mean_tissue_viability",
                "op": "final_gte",
                "threshold": 0.7,
                "unit": "normalized",
                "rationale": "Tissue viability restored above ischemic threshold",
            },
        ],
    },
    "burn": {
        "label": "Burn Wound Resolution",
        "checks": [
            {
                "name": "Complete closure",
                "metric": "wound_closure_pct",
                "op": "final_gte",
                "threshold": 95.0,
                "unit": "%",
                "rationale": "Full closure including deep partial-thickness areas",
            },
            {
                "name": "Inflammation resolved",
                "metric": "mean_infl_wound",
                "op": "tail_mean_lt",
                "threshold": 0.06,
                "tail_frac": 0.2,
                "unit": "normalized",
                "rationale": "No persistent thermal inflammation",
            },
            {
                "name": "Scar minimized",
                "metric": "scar_magnitude",
                "op": "final_lt",
                "threshold": 0.3,
                "unit": "normalized",
                "rationale": "Hypertrophic scarring prevented (major burn morbidity)",
            },
            {
                "name": "Collagen quality",
                "metric": "mean_collagen_wound",
                "op": "final_gte",
                "threshold": 0.4,
                "unit": "normalized",
                "rationale": "Functional collagen matrix deposited in burn bed",
            },
        ],
    },
    "rheumatoid": {
        "label": "Rheumatoid Arthritis Remission",
        "checks": [
            {
                "name": "TNF-alpha suppressed",
                "metric": "mean_tnf_alpha_wound",
                "op": "tail_mean_lt",
                "threshold": 0.05,
                "tail_frac": 0.25,
                "unit": "normalized",
                "rationale": "TNF-alpha at near-zero sustained level (not just acute suppression)",
            },
            {
                "name": "IL-6 suppressed",
                "metric": "mean_il6_wound",
                "op": "tail_mean_lt",
                "threshold": 0.05,
                "tail_frac": 0.25,
                "unit": "normalized",
                "rationale": "IL-6 cascade resolved",
            },
            {
                "name": "Cartilage preserved",
                "metric": "mean_cartilage_wound",
                "op": "final_gte",
                "threshold": 0.85,
                "unit": "integrity",
                "rationale": "Minimal irreversible cartilage erosion (<15% loss)",
            },
            {
                "name": "Bone preserved",
                "metric": "mean_bone_wound",
                "op": "final_gte",
                "threshold": 0.90,
                "unit": "integrity",
                "rationale": "Subchondral bone structure maintained",
            },
            {
                "name": "Pannus regressed",
                "metric": "mean_synovial_wound",
                "op": "tail_trend_negative",
                "threshold": 0.0,
                "tail_frac": 0.3,
                "unit": "trend",
                "rationale": "Synovial pannus actively regressing (not just stabilized)",
            },
        ],
    },
}


# ---------------------------------------------------------------------------
# Evaluation engine
# ---------------------------------------------------------------------------

def _check_criterion(data, check):
    """Evaluate one criterion against simulation data. Returns (met, detail)."""
    metric = check["metric"]
    if metric not in data or not data[metric]:
        return False, f"{metric} not available"

    values = data[metric]
    op = check["op"]
    threshold = check["threshold"]

    if op == "final_gte":
        val = values[-1]
        met = val >= threshold
        return met, f"{val:.2f} {'>=':s} {threshold} {check['unit']}"

    if op == "final_lt":
        val = values[-1]
        met = val < threshold
        return met, f"{val:.2f} {'<':s} {threshold} {check['unit']}"

    if op == "tail_mean_lt":
        tail_frac = check.get("tail_frac", 0.25)
        n_tail = max(1, int(len(values) * tail_frac))
        tail = values[-n_tail:]
        mean_val = sum(tail) / len(tail)
        met = mean_val < threshold
        return met, f"tail mean {mean_val:.4f} {'<':s} {threshold}"

    if op == "tail_trend_negative":
        tail_frac = check.get("tail_frac", 0.3)
        n_tail = max(2, int(len(values) * tail_frac))
        tail = values[-n_tail:]
        # Simple linear regression slope
        n = len(tail)
        x_mean = (n - 1) / 2.0
        y_mean = sum(tail) / n
        num = sum((i - x_mean) * (v - y_mean) for i, v in enumerate(tail))
        den = sum((i - x_mean) ** 2 for i in range(n))
        slope = num / den if den > 0 else 0
        met = slope < 0
        return met, f"slope {slope:.6f} ({'declining' if met else 'rising/flat'})"

    return False, f"unknown op: {op}"


def evaluate_cure(data, condition):
    """Evaluate cure criteria for a condition against simulation data.

    Args:
        data: dict of metric_name -> list of values (from load_csv or aggregate)
        condition: one of "diabetic", "pressure", "burn", "rheumatoid"

    Returns:
        dict with keys:
            label: human-readable condition label
            verdict: "CURED", "PARTIAL", or "UNCURED"
            criteria: list of dicts with name, met (bool), detail, rationale
            score: fraction of criteria met (0.0 to 1.0)
    """
    if condition not in CRITERIA:
        return {
            "label": condition,
            "verdict": "N/A",
            "criteria": [],
            "score": 0.0,
        }

    spec = CRITERIA[condition]
    results = []
    for check in spec["checks"]:
        met, detail = _check_criterion(data, check)
        results.append({
            "name": check["name"],
            "met": met,
            "detail": detail,
            "rationale": check["rationale"],
        })

    n_met = sum(1 for r in results if r["met"])
    n_total = len(results)
    score = n_met / n_total if n_total > 0 else 0.0

    if n_met == n_total:
        verdict = "CURED"
    elif n_met >= n_total * 0.6:
        verdict = "PARTIAL"
    else:
        verdict = "UNCURED"

    return {
        "label": spec["label"],
        "verdict": verdict,
        "criteria": results,
        "score": score,
    }


def print_cure_assessment(result):
    """Print a formatted cure assessment to stdout."""
    label = result["label"]
    verdict = result["verdict"]
    score = result["score"]

    verdict_marker = {
        "CURED": "[CURED]",
        "PARTIAL": "[PARTIAL]",
        "UNCURED": "[UNCURED]",
        "N/A": "[N/A]",
    }

    print(f"\n  {label}: {verdict_marker.get(verdict, verdict)} ({score:.0%})")
    for c in result["criteria"]:
        status = "PASS" if c["met"] else "FAIL"
        print(f"    [{status}] {c['name']}: {c['detail']}")
