"""Compute mean RMSE across replicates and compare WITH vs WITHOUT mechanisms.

Reads /tmp/replicate_runs/metrics_seed*.csv, runs validation on each,
prints mean and stddev per observable.
"""
import csv
import glob
import os
import statistics
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "literature"))

from lib import (load_csv, detect_condition, detect_modules,
                 validate_wound, validate_fibroblast,
                 validate_microenvironment, validate_ph)


def rmse_dict(path):
    sim = load_csv(path)
    days = [h / 24.0 for h in sim["time_h"]]
    has_wound, has_fibro, _, has_micro, has_ph, _ = detect_modules(sim)
    out = {}
    if has_wound:
        w = validate_wound(sim, days, "normal")
        out["closure"] = w["closure_rmse"]
        out["inflammation"] = w["inflammation_rmse"] * 100
        out["neutrophils"] = w["neut_rmse"] * 100
        out["macrophages"] = w["mac_rmse"] * 100
    if has_fibro:
        f = validate_fibroblast(sim, days)
        out["fibroblasts"] = f["fibro_rmse"] * 100
        out["myofibroblasts"] = f["myofib_rmse"] * 100
        out["collagen"] = f["collagen_rmse"] * 100
    if has_micro:
        m = validate_microenvironment(sim, days, "normal")
        out["tgfb"] = m["tgfb_rmse"] * 100
        out["vegf"] = m["vegf_rmse"] * 100
        out["fibronectin"] = m["fn_rmse"] * 100
        out["mmp"] = m["mmp_rmse"] * 100
    if has_ph:
        p = validate_ph(sim, days)
        out["ph"] = p["ph_rmse"] * 100
    return out


def main():
    paths = sorted(glob.glob("/tmp/replicate_runs/metrics_seed*.csv"))
    if not paths:
        print("No replicate CSVs found")
        return 1

    results = [rmse_dict(p) for p in paths]
    keys = sorted(results[0].keys())

    print(f"=== {len(results)} replicates ===")
    print(f"{'observable':<16} {'mean':>7} {'stddev':>7}  per-seed")
    for k in keys:
        vals = [r[k] for r in results]
        m = statistics.mean(vals)
        s = statistics.stdev(vals) if len(vals) > 1 else 0
        per = " ".join(f"{v:6.2f}" for v in vals)
        marker = " <-- > 15%" if m > 15 else ""
        print(f"{k:<16} {m:7.2f} {s:7.2f}  {per}{marker}")


if __name__ == "__main__":
    sys.exit(main() or 0)
