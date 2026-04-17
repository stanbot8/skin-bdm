"""Compare mean RMSE WITH mechanisms vs WITHOUT across replicate seeds."""
import glob
import os
import statistics
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "literature"))
from lib import (load_csv, detect_modules, validate_wound,
                 validate_fibroblast, validate_microenvironment, validate_ph)


def rmse(path):
    sim = load_csv(path)
    days = [h / 24.0 for h in sim["time_h"]]
    hw, hf, _, hm, hp, _ = detect_modules(sim)
    r = {}
    if hw:
        w = validate_wound(sim, days, "normal")
        r["closure"] = w["closure_rmse"]
        r["inflammation"] = w["inflammation_rmse"] * 100
        r["neutrophils"] = w["neut_rmse"] * 100
        r["macrophages"] = w["mac_rmse"] * 100
    if hf:
        f = validate_fibroblast(sim, days)
        r["fibroblasts"] = f["fibro_rmse"] * 100
        r["myofibroblasts"] = f["myofib_rmse"] * 100
        r["collagen"] = f["collagen_rmse"] * 100
    if hm:
        m = validate_microenvironment(sim, days, "normal")
        r["tgfb"] = m["tgfb_rmse"] * 100
        r["vegf"] = m["vegf_rmse"] * 100
        r["fibronectin"] = m["fn_rmse"] * 100
        r["mmp"] = m["mmp_rmse"] * 100
    if hp:
        p = validate_ph(sim, days)
        r["ph"] = p["ph_rmse"] * 100
    return r


def fleet_stats(dir_path):
    paths = sorted(glob.glob(os.path.join(dir_path, "metrics_seed*.csv")))
    if not paths:
        return None
    per_seed = [rmse(p) for p in paths]
    keys = sorted(per_seed[0].keys())
    return {k: statistics.mean(r[k] for r in per_seed) for k in keys}


def main():
    with_m = fleet_stats("/tmp/replicate_runs")
    no_m = fleet_stats("/tmp/baseline_runs")
    if with_m is None or no_m is None:
        print("missing replicate runs")
        return 1

    keys = sorted(with_m.keys())
    print(f"{'observable':<16} {'WITHOUT':>8} {'WITH':>8}  {'delta':>7}  verdict")
    for k in keys:
        w = with_m[k]
        b = no_m[k]
        d = w - b
        verdict = "better" if d < -0.3 else ("worse" if d > 0.3 else "~same")
        flag = " !" if w > 15 else ""
        print(f"{k:<16} {b:8.2f} {w:8.2f}  {d:+7.2f}  {verdict}{flag}")


if __name__ == "__main__":
    sys.exit(main() or 0)
