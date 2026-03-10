#!/usr/bin/env python3
"""SkiBiDy Dashboard: lightweight dark-theme GUI.

Pure Python stdlib server + vanilla JS + Chart.js (CDN).
No Streamlit, no heavy deps. Runs on any Python 3.6+.

Usage:  python3 scripts/dashboard.py [--port 8501]
"""

import csv
import glob
import http.server
import json
import os
import shutil
import subprocess
import sys
import threading
import urllib.parse
import webbrowser

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "literature"))
PORT = 8501

# User studies directory: SKIBIDY_USER_STUDIES env var, or ~/.skibidy/studies/
USER_STUDIES_DIR = os.environ.get(
    "SKIBIDY_USER_STUDIES",
    os.path.join(os.path.expanduser("~"), ".skibidy", "studies"))
ENGINE_STUDIES_DIR = os.path.join(ROOT, "studies")

# ---------------------------------------------------------------------------
# Data helpers
# ---------------------------------------------------------------------------

def load_csv(path):
    with open(path) as f:
        reader = csv.DictReader(row for row in f if not row.startswith("#"))
        rows = list(reader)
    if not rows:
        return {}
    data = {}
    for key in rows[0]:
        try:
            data[key] = [float(r[key]) for r in rows]
        except ValueError:
            data[key] = [r[key] for r in rows]
    return data


def _read_skibidy_meta(study_dir, name):
    """Read .skibidy project file and return metadata dict."""
    proj = os.path.join(study_dir, name + ".skibidy")
    if not os.path.isfile(proj):
        return {"name": name, "description": "", "category": ""}
    meta = {"name": name, "description": "", "category": ""}
    for line in open(proj):
        line = line.strip()
        if line.startswith("description"):
            meta["description"] = line.split("=", 1)[1].strip().strip('"')
        elif line.startswith("category"):
            meta["category"] = line.split("=", 1)[1].strip().strip('"')
    return meta


def _scan_study_dir(base_dir, scope):
    """Scan a directory for studies, returning list of study dicts."""
    if not os.path.isdir(base_dir):
        return []
    skip = {"scenarios", "__pycache__", "shared"}
    out = []
    for n in os.listdir(base_dir):
        p = os.path.join(base_dir, n)
        if not os.path.isdir(p) or n in skip or n.startswith("."):
            continue
        meta = _read_skibidy_meta(p, n)
        out.append({
            "name": n, "path": p, "scope": scope,
            "description": meta.get("description", ""),
            "category": meta.get("category", ""),
            "mtime": os.path.getmtime(p),
        })
    out.sort(key=lambda x: x["mtime"], reverse=True)
    return out


def find_studies():
    """Return dict {name: path} for all studies (engine + user)."""
    engine = _scan_study_dir(ENGINE_STUDIES_DIR, "engine")
    user = _scan_study_dir(USER_STUDIES_DIR, "user")
    combined = {}
    for s in engine + user:
        combined[s["name"]] = s["path"]
    return combined


def find_studies_detailed():
    """Return list of study dicts with scope, description, category."""
    engine = _scan_study_dir(ENGINE_STUDIES_DIR, "engine")
    user = _scan_study_dir(USER_STUDIES_DIR, "user")
    return {"engine": engine, "user": user}


def find_results(study_dir):
    results = []
    batch_dir = os.path.join(ROOT, "batch", "results")
    study_name = os.path.basename(study_dir)
    if os.path.isdir(batch_dir):
        for d in sorted(os.listdir(batch_dir), reverse=True):
            full = os.path.join(batch_dir, d)
            if not os.path.isdir(full) or study_name not in d:
                continue
            cp = os.path.join(full, "metrics.csv")
            if os.path.isfile(cp):
                results.append({"type": "batch", "label": d, "csv": cp, "dir": full})
    res_dir = os.path.join(study_dir, "results")
    if os.path.isdir(res_dir):
        cp = os.path.join(res_dir, "metrics.csv")
        if os.path.isfile(cp):
            results.append({"type": "study", "label": "latest", "csv": cp, "dir": res_dir})
        exp_dir = os.path.join(res_dir, "experiments")
        if os.path.isdir(exp_dir):
            for exp in sorted(os.listdir(exp_dir)):
                ep = os.path.join(exp_dir, exp)
                if not os.path.isdir(ep):
                    continue
                for cfg in sorted(os.listdir(ep)):
                    cfgp = os.path.join(ep, cfg)
                    for cn in ("consensus_metrics.csv", "metrics.csv"):
                        cp = os.path.join(cfgp, cn)
                        if os.path.isfile(cp):
                            results.append({"type": "experiment",
                                            "label": f"{exp}/{cfg}",
                                            "csv": cp, "dir": cfgp})
                            break
    return results


def find_experiments(study_dir):
    d = os.path.join(study_dir, "experiments")
    if not os.path.isdir(d):
        return []
    out = []
    for f in sorted(glob.glob(os.path.join(d, "*.toml"))):
        name = os.path.splitext(os.path.basename(f))[0]
        with open(f) as fh:
            out.append({"name": name, "content": fh.read()})
    return out


def find_treatments(study_dir):
    out = []
    for t in sorted(glob.glob(os.path.join(ROOT, "treatments", "*.toml"))):
        name = os.path.splitext(os.path.basename(t))[0]
        with open(t) as f:
            out.append({"name": name, "scope": "engine", "content": f.read()})
    for t in sorted(glob.glob(os.path.join(study_dir, "treatments", "*.toml"))):
        name = os.path.splitext(os.path.basename(t))[0]
        with open(t) as f:
            out.append({"name": name, "scope": "study", "content": f.read()})
    return out


def find_modules():
    out = []
    d = os.path.join(ROOT, "modules")
    if not os.path.isdir(d):
        return out
    for mod in sorted(os.listdir(d)):
        cfg = os.path.join(d, mod, "config.toml")
        if os.path.isfile(cfg):
            with open(cfg) as f:
                out.append({"name": mod, "content": f.read()})
    return out


def find_profiles():
    d = os.path.join(ROOT, "profiles")
    if not os.path.isdir(d):
        return []
    out = []
    for f in sorted(os.listdir(d)):
        if f.endswith(".toml") and f != "TEMPLATE.toml":
            out.append(os.path.splitext(f)[0])
    return out


def create_study(name, duration, profile, modules_list, description=""):
    """Create study directory scaffold in user studies directory."""
    safe = name.strip().lower().replace(" ", "-")
    safe = "".join(c for c in safe if c.isalnum() or c in "-_")
    if not safe:
        return {"ok": False, "error": "Invalid study name"}
    # Check both engine and user dirs for conflicts
    if os.path.exists(os.path.join(ENGINE_STUDIES_DIR, safe)):
        return {"ok": False, "error": f"Engine study '{safe}' already exists"}
    study_dir = os.path.join(USER_STUDIES_DIR, safe)
    if os.path.exists(study_dir):
        return {"ok": False, "error": f"Study '{safe}' already exists"}
    os.makedirs(os.path.join(study_dir, "experiments"), exist_ok=True)
    os.makedirs(os.path.join(study_dir, "treatments"), exist_ok=True)
    os.makedirs(os.path.join(study_dir, "results"), exist_ok=True)
    # Write .skibidy project file
    with open(os.path.join(study_dir, safe + ".skibidy"), "w") as f:
        f.write("[project]\n")
        f.write(f'name = "{safe}"\n')
        f.write(f'description = "{description}"\n')
        f.write('category = "custom"\n')
    # Write preset.toml
    lines = [f"# Study: {safe}\n"]
    lines.append(f"\n[skin]\n")
    lines.append(f"duration_days = {int(duration)}\n")
    for mod in modules_list:
        section = f"skin.{mod}"
        lines.append(f"\n[{section}]\n")
        lines.append("enabled = true\n")
    with open(os.path.join(study_dir, "preset.toml"), "w") as f:
        f.writelines(lines)
    return {"ok": True, "study": safe, "profile": profile, "scope": "user"}


def duplicate_study(source_name, new_name):
    """Duplicate an engine study into user studies directory."""
    safe = new_name.strip().lower().replace(" ", "-")
    safe = "".join(c for c in safe if c.isalnum() or c in "-_")
    if not safe:
        return {"ok": False, "error": "Invalid study name"}
    # Find the source study
    all_studies = find_studies()
    if source_name not in all_studies:
        return {"ok": False, "error": f"Source study '{source_name}' not found"}
    src_dir = all_studies[source_name]
    dest_dir = os.path.join(USER_STUDIES_DIR, safe)
    if os.path.exists(dest_dir):
        return {"ok": False, "error": f"Study '{safe}' already exists"}
    os.makedirs(USER_STUDIES_DIR, exist_ok=True)
    # Copy everything except results
    shutil.copytree(src_dir, dest_dir, ignore=shutil.ignore_patterns(
        "results", "__pycache__", "*.pyc"))
    os.makedirs(os.path.join(dest_dir, "results"), exist_ok=True)
    # Update .skibidy project file
    old_proj = os.path.join(dest_dir, source_name + ".skibidy")
    new_proj = os.path.join(dest_dir, safe + ".skibidy")
    if os.path.isfile(old_proj):
        with open(old_proj) as f:
            content = f.read()
        content = content.replace(f'name = "{source_name}"', f'name = "{safe}"')
        with open(new_proj, "w") as f:
            f.write(content)
        if old_proj != new_proj:
            os.remove(old_proj)
    else:
        with open(new_proj, "w") as f:
            f.write("[project]\n")
            f.write(f'name = "{safe}"\n')
            f.write(f'description = "Fork of {source_name}"\n')
            f.write('category = "custom"\n')
    return {"ok": True, "study": safe, "source": source_name}


def create_experiment(study_name, exp_name, description, profile, runs,
                      configs):
    """Create an experiment TOML in the study's experiments/ dir."""
    studies = find_studies()
    if study_name not in studies:
        return {"ok": False, "error": f"Study '{study_name}' not found"}
    safe = exp_name.strip().lower().replace(" ", "_")
    safe = "".join(c for c in safe if c.isalnum() or c == "_")
    if not safe:
        return {"ok": False, "error": "Invalid experiment name"}
    exp_dir = os.path.join(studies[study_name], "experiments")
    os.makedirs(exp_dir, exist_ok=True)
    path = os.path.join(exp_dir, safe + ".toml")
    if os.path.exists(path):
        return {"ok": False, "error": f"Experiment '{safe}' already exists"}
    lines = [f'[experiment]\n']
    lines.append(f'name = "{exp_name}"\n')
    if description:
        lines.append(f'description = "{description}"\n')
    if profile:
        lines.append(f'profile = "{profile}"\n')
    lines.append(f'study = "{study_name}"\n')
    lines.append(f'runs_per_config = {int(runs)}\n')
    for cfg in configs:
        lines.append(f'\n[[experiment.configs]]\n')
        lines.append(f'label = "{cfg.get("label", "config")}"\n')
        treats = cfg.get("treatments", [])
        if treats:
            lines.append(f'treatments = {json.dumps(treats)}\n')
        overrides = cfg.get("overrides", {})
        if overrides:
            lines.append(f'[experiment.configs.overrides]\n')
            for k, v in overrides.items():
                if isinstance(v, str):
                    lines.append(f'"{k}" = "{v}"\n')
                else:
                    lines.append(f'"{k}" = {v}\n')
    with open(path, "w") as f:
        f.writelines(lines)
    return {"ok": True, "name": safe}


def create_treatment(study_name, treat_name, content):
    """Create a treatment TOML in the study's treatments/ dir."""
    studies = find_studies()
    if study_name not in studies:
        return {"ok": False, "error": f"Study '{study_name}' not found"}
    safe = treat_name.strip().lower().replace(" ", "_")
    safe = "".join(c for c in safe if c.isalnum() or c == "_")
    if not safe:
        return {"ok": False, "error": "Invalid treatment name"}
    treat_dir = os.path.join(studies[study_name], "treatments")
    os.makedirs(treat_dir, exist_ok=True)
    path = os.path.join(treat_dir, safe + ".toml")
    if os.path.exists(path):
        return {"ok": False, "error": f"Treatment '{safe}' already exists"}
    with open(path, "w") as f:
        f.write(content)
    return {"ok": True, "name": safe}


def run_validation(csv_path):
    try:
        from lib import (detect_modules, validate_wound, validate_fibroblast,
                         validate_microenvironment, validate_ph, validate_ra)
        data = load_csv(csv_path)
        if not data:
            return {"error": "Empty CSV"}
        sim_days = [h / 24.0 for h in data["time_h"]]
        hw, hf, ht, hm, hp, hr = detect_modules(data)
        res = {}
        if hw:
            w = validate_wound(data, sim_days)
            res["Wound closure"] = w["closure_rmse"] / 100.0
            res["Inflammation"] = w["inflammation_rmse"]
            res["Neutrophils"] = w["neut_rmse"]
            res["Macrophages"] = w["mac_rmse"]
        if hf:
            f = validate_fibroblast(data, sim_days)
            res["Fibroblasts"] = f["fibro_rmse"]
            res["Myofibroblasts"] = f["myofib_rmse"]
            res["Collagen"] = f["collagen_rmse"]
        if hm:
            m = validate_microenvironment(data, sim_days)
            res[u"TGF-\u03b2"] = m["tgfb_rmse"]
            res["VEGF"] = m["vegf_rmse"]
            res["Fibronectin"] = m["fn_rmse"]
            res["MMP"] = m["mmp_rmse"]
        if hp:
            res["Wound pH"] = validate_ph(data, sim_days)["ph_rmse"]
        if hr:
            a = validate_ra(data, sim_days)
            res[u"TNF-\u03b1"] = a["tnf_rmse"]
            res["IL-6"] = a["il6_rmse"]
            res["Cartilage"] = a["cart_rmse"]
            if a.get("has_bone"):
                res["Bone erosion"] = a["bone_rmse"]
            if a.get("has_tcell"):
                res["T cells"] = a["tcell_rmse"]
            if a.get("has_syn"):
                res["Synovial pannus"] = a["syn_rmse"]
        return res
    except Exception as e:
        return {"error": str(e)}


def _find_paraview():
    pv = shutil.which("paraview")
    if pv:
        return pv
    pv_dir = os.environ.get("ParaView_DIR", "")
    if pv_dir:
        c = os.path.join(pv_dir, "bin", "paraview")
        if os.path.isfile(c):
            return c
    return None


def launch_paraview(output_dir):
    pv = _find_paraview()
    if not pv:
        return False
    # Look for .pvd files in the output directory tree
    pvd = os.path.join(output_dir, "output.pvd")
    if not os.path.isfile(pvd):
        # Search for any .pvd file
        pvds = glob.glob(os.path.join(output_dir, "**", "*.pvd"), recursive=True)
        if pvds:
            pvd = pvds[0]
        else:
            return False
    subprocess.Popen([pv, pvd], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return True


# ---------------------------------------------------------------------------
# API handler
# ---------------------------------------------------------------------------

_run_state = {"proc": None, "log": ""}


class DashHandler(http.server.BaseHTTPRequestHandler):
    def log_message(self, fmt, *args):
        pass  # silence request logs

    def _json(self, obj, code=200):
        data = json.dumps(obj).encode()
        self.send_response(code)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", len(data))
        self.end_headers()
        self.wfile.write(data)

    def _html(self, content):
        data = content.encode()
        self.send_response(200)
        self.send_header("Content-Type", "text/html; charset=utf-8")
        self.send_header("Content-Length", len(data))
        self.end_headers()
        self.wfile.write(data)

    def do_GET(self):
        parsed = urllib.parse.urlparse(self.path)
        path = parsed.path
        qs = dict(urllib.parse.parse_qsl(parsed.query))

        if path == "/":
            self._html(HTML_PAGE)
        elif path == "/api/studies":
            self._json(find_studies_detailed())
        elif path == "/api/user-studies-dir":
            self._json({"path": USER_STUDIES_DIR})
        elif path == "/api/results":
            study = qs.get("study", "")
            studies = find_studies()
            if study not in studies:
                self._json([])
                return
            self._json(find_results(studies[study]))
        elif path == "/api/csv":
            csv_path = qs.get("path", "")
            if not csv_path or not os.path.isfile(csv_path):
                self._json({"error": "not found"}, 404)
                return
            data = load_csv(csv_path)
            # Only send columns that have data, keep payload small
            skip = set()
            for k, v in data.items():
                if isinstance(v, list) and all(x == 0 for x in v):
                    skip.add(k)
            trimmed = {k: v for k, v in data.items() if k not in skip}
            self._json(trimmed)
        elif path == "/api/validate":
            csv_path = qs.get("path", "")
            if not csv_path or not os.path.isfile(csv_path):
                self._json({"error": "not found"}, 404)
                return
            self._json(run_validation(csv_path))
        elif path == "/api/experiments":
            study = qs.get("study", "")
            studies = find_studies()
            if study not in studies:
                self._json([])
                return
            self._json(find_experiments(studies[study]))
        elif path == "/api/treatments":
            study = qs.get("study", "")
            studies = find_studies()
            if study not in studies:
                self._json([])
                return
            self._json(find_treatments(studies[study]))
        elif path == "/api/modules":
            self._json(find_modules())
        elif path == "/api/profiles":
            self._json(find_profiles())
        elif path == "/api/paraview":
            d = qs.get("dir", "")
            ok = launch_paraview(d) if d else False
            self._json({"ok": ok})
        elif path == "/api/build":
            rc = subprocess.run(
                ["bash", "-c",
                 "source ~/biodynamo/build/bin/thisbdm.sh && biodynamo build 2>&1"],
                capture_output=True, text=True, cwd=ROOT, timeout=300)
            self._json({"ok": rc.returncode == 0,
                        "log": (rc.stdout + rc.stderr)[-3000:]})
        elif path == "/api/run":
            study = qs.get("study", "")
            n = qs.get("n", "1")
            skin = qs.get("skin", "")
            treatment = qs.get("treatment", "")
            if _run_state["proc"] and _run_state["proc"].poll() is None:
                self._json({"ok": False, "error": "A run is already in progress"})
                return
            cmd = ("source ~/biodynamo/build/bin/thisbdm.sh && "
                   f"python3 batch/batch.py -n {n}")
            if study:
                cmd += f" --study {study}"
            if skin:
                cmd += f" --skin {skin}"
            if treatment:
                cmd += f" --treatment {treatment}"
            _run_state["log"] = ""
            _run_state["proc"] = subprocess.Popen(
                ["bash", "-c", cmd],
                cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                text=True)
            # Read output in background thread
            def _reader(proc):
                for line in proc.stdout:
                    _run_state["log"] += line
                proc.wait()
            threading.Thread(target=_reader, args=(_run_state["proc"],),
                             daemon=True).start()
            self._json({"ok": True})
        elif path == "/api/run/status":
            proc = _run_state["proc"]
            running = proc is not None and proc.poll() is None
            rc = proc.returncode if proc and not running else None
            self._json({"running": running, "returncode": rc,
                        "log": _run_state["log"][-4000:]})
        else:
            self.send_error(404)

    def do_POST(self):
        parsed = urllib.parse.urlparse(self.path)
        path = parsed.path
        length = int(self.headers.get("Content-Length", 0))
        body = json.loads(self.rfile.read(length)) if length else {}

        if path == "/api/create-study":
            res = create_study(
                body.get("name", ""),
                body.get("duration", 30),
                body.get("profile", ""),
                body.get("modules", []),
                body.get("description", ""))
            self._json(res)
        elif path == "/api/duplicate-study":
            res = duplicate_study(
                body.get("source", ""),
                body.get("name", ""))
            self._json(res)
        elif path == "/api/create-experiment":
            res = create_experiment(
                body.get("study", ""),
                body.get("name", ""),
                body.get("description", ""),
                body.get("profile", ""),
                body.get("runs", 5),
                body.get("configs", []))
            self._json(res)
        elif path == "/api/create-treatment":
            res = create_treatment(
                body.get("study", ""),
                body.get("name", ""),
                body.get("content", ""))
            self._json(res)
        else:
            self.send_error(404)


# ---------------------------------------------------------------------------
# HTML (single page, dark theme, Chart.js)
# ---------------------------------------------------------------------------

HTML_PAGE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>SkiBiDy Dashboard</title>
<script src="https://cdn.jsdelivr.net/npm/chart.js@4/dist/chart.umd.min.js"></script>
<style>
*{margin:0;padding:0;box-sizing:border-box}
:root{
  --bg:#0d1117;--bg2:#161b22;--bg3:#21262d;
  --border:#30363d;--text:#c9d1d9;--text2:#8b949e;
  --blue:#58a6ff;--green:#3fb950;--yellow:#d29922;
  --red:#f85149;--purple:#bc8cff;--white:#f0f6fc;
}
body{background:var(--bg);color:var(--text);font:14px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",Helvetica,Arial,sans-serif}
a{color:var(--blue);text-decoration:none}

/* Layout */
.shell{display:flex;height:100vh;overflow:hidden}
.sidebar{width:220px;min-width:220px;background:var(--bg2);border-right:1px solid var(--border);padding:16px;overflow-y:auto;display:flex;flex-direction:column}
.main{flex:1;overflow-y:auto;padding:24px 28px}

/* Sidebar */
.logo{font-size:18px;font-weight:600;color:var(--white);margin-bottom:2px}
.logo-sub{font-size:11px;color:var(--text2);margin-bottom:20px}
.study-section{margin-bottom:14px}
.study-section-header{display:flex;justify-content:space-between;align-items:center;padding:4px 0;cursor:pointer;user-select:none}
.study-section-header span{font-size:10px;text-transform:uppercase;letter-spacing:.5px;color:var(--text2);font-weight:600}
.study-section-toggle{color:var(--text2);font-size:10px;transition:transform .15s;background:none;border:none;cursor:pointer}
.study-section.collapsed .study-list{display:none}
.study-section.collapsed .study-section-toggle{transform:rotate(-90deg)}
.study-list{list-style:none;margin:0;padding:0}
.study-list li{padding:5px 10px;border-radius:6px;cursor:pointer;font-size:13px;color:var(--text);transition:background .1s;display:flex;align-items:center;justify-content:space-between;gap:4px}
.study-list li:hover{background:var(--bg3)}
.study-list li.active{background:var(--blue);color:#fff}
.study-list li .s-name{flex:1;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}
.study-list li .s-desc{display:none}
.study-list li:hover .s-desc{display:block;font-size:10px;color:var(--text2);position:absolute;left:230px;background:var(--bg2);border:1px solid var(--border);padding:4px 8px;border-radius:4px;white-space:nowrap;z-index:50}
.btn-dup{background:none;border:none;color:var(--text2);cursor:pointer;font-size:11px;padding:2px 4px;border-radius:3px;opacity:0;transition:all .12s}
.study-list li:hover .btn-dup{opacity:1}
.btn-dup:hover{color:var(--blue);background:var(--bg)}
.user-dir-label{font-size:10px;color:var(--text2);padding:6px 0 2px;word-break:break-all}
.btn{padding:6px 12px;border-radius:6px;border:1px solid var(--border);background:var(--bg3);color:var(--text);cursor:pointer;font-size:12px;transition:all .15s;text-align:center;flex:1}
.btn:hover{background:var(--blue);color:#fff;border-color:var(--blue)}
.btn.sm{padding:4px 8px;font-size:11px}

/* Action panel (bottom right) */
.action-panel{position:fixed;bottom:20px;right:20px;background:var(--bg2);border:1px solid var(--border);border-radius:10px;width:260px;z-index:100;box-shadow:0 4px 20px rgba(0,0,0,.4);overflow:hidden;transition:all .2s}
.action-panel.collapsed{width:44px;height:44px;border-radius:50%;cursor:pointer}
.action-panel.collapsed .ap-body{display:none}
.action-panel.collapsed .ap-toggle{transform:rotate(180deg)}
.ap-header{display:flex;justify-content:space-between;align-items:center;padding:10px 14px;cursor:pointer;user-select:none}
.ap-header .ap-title{font-size:11px;color:var(--text2);text-transform:uppercase;letter-spacing:.5px;margin:0}
.ap-toggle{color:var(--text2);font-size:14px;transition:transform .2s;background:none;border:none;cursor:pointer;padding:0}
.action-panel.collapsed .ap-header{padding:10px 12px;justify-content:center}
.action-panel.collapsed .ap-title{display:none}
.ap-body{padding:0 14px 14px}
.ap-body .ap-row{display:flex;gap:6px;margin-bottom:6px}
.ap-body .ap-row:last-child{margin-bottom:0}
.action-panel select,.action-panel input[type=text]{font-size:12px;padding:5px 8px}
.action-panel .btn{font-size:12px}

/* Tabs */
.tabs{display:flex;gap:0;background:var(--bg2);border-radius:8px;padding:3px;margin-bottom:20px}
.tab{padding:7px 18px;border-radius:6px;cursor:pointer;font-size:13px;color:var(--text2);transition:all .12s;user-select:none}
.tab:hover{color:var(--text)}
.tab.active{background:var(--blue);color:#fff}
.tab-panel{display:none}
.tab-panel.active{display:block}

/* Cards */
.cards{display:grid;grid-template-columns:repeat(auto-fill,minmax(180px,1fr));gap:12px;margin-bottom:20px}
.card{background:var(--bg2);border:1px solid var(--border);border-radius:8px;padding:14px}
.card-label{font-size:12px;color:var(--blue);margin-bottom:2px}
.card-value{font-size:24px;font-weight:600;color:var(--white)}
.card-sub{font-size:11px;color:var(--text2)}

/* Charts */
.chart-grid{display:grid;grid-template-columns:repeat(3,1fr);gap:16px;margin-bottom:20px}
.chart-box{background:var(--bg2);border:1px solid var(--border);border-radius:8px;padding:12px}
.chart-title{font-size:12px;color:var(--text2);margin-bottom:8px;text-transform:capitalize}
canvas{width:100%!important;height:200px!important}

/* RMSE table */
.rmse-table{background:var(--bg2);border:1px solid var(--border);border-radius:8px;overflow:hidden;margin-bottom:20px}
.rmse-row{display:flex;justify-content:space-between;align-items:center;padding:8px 14px;border-bottom:1px solid var(--border);font-size:13px}
.rmse-row:last-child{border-bottom:none}
.rmse-row:hover{background:var(--bg3)}
.pill{display:inline-block;padding:1px 8px;border-radius:10px;font-size:11px;font-weight:500;margin-left:8px}
.pill-green{background:#238636;color:#fff}
.pill-yellow{background:#9e6a03;color:#fff}
.pill-red{background:#da3633;color:#fff}
.rmse-good{color:var(--green)}.rmse-ok{color:var(--yellow)}.rmse-bad{color:var(--red)}

/* Select / multiselect */
select,input[type=text]{background:var(--bg);border:1px solid var(--border);color:var(--text);border-radius:6px;padding:6px 10px;font-size:13px;width:100%}
select:focus,input:focus{outline:none;border-color:var(--blue)}
.field{margin-bottom:12px}
.field label{font-size:12px;color:var(--text2);display:block;margin-bottom:4px}
.chip-bar{display:flex;flex-wrap:wrap;gap:4px;margin-bottom:16px}
.chip{padding:3px 10px;border-radius:12px;font-size:11px;cursor:pointer;border:1px solid var(--border);background:var(--bg3);color:var(--text);transition:all .12s;user-select:none}
.chip.on{background:var(--blue);color:#fff;border-color:var(--blue)}

/* Toml code */
.toml-block{background:var(--bg);border:1px solid var(--border);border-radius:6px;padding:12px;font-family:"Fira Code",monospace;font-size:12px;white-space:pre-wrap;color:var(--text);max-height:300px;overflow-y:auto;margin-bottom:8px}
.exp-header{padding:10px 14px;background:var(--bg2);border:1px solid var(--border);border-radius:6px;cursor:pointer;margin-bottom:6px;font-size:13px;display:flex;justify-content:space-between;align-items:center}
.exp-header:hover{background:var(--bg3)}
.exp-body{display:none;margin-bottom:10px}
.exp-body.open{display:block}

/* Run selector */
.run-bar{display:flex;gap:8px;align-items:center;margin-bottom:16px}
.run-bar select{max-width:400px}

/* Info */
.info-box{background:var(--bg2);border:1px solid var(--border);border-radius:8px;padding:20px;text-align:center;color:var(--text2)}

/* Toast */
.toast{position:fixed;bottom:200px;right:20px;background:var(--green);color:#fff;padding:8px 16px;border-radius:6px;font-size:13px;opacity:0;transition:opacity .3s;pointer-events:none;z-index:999}
.toast.show{opacity:1}

/* Config columns */
.config-grid{display:grid;grid-template-columns:1fr 1fr;gap:20px}
.config-col h3{font-size:14px;color:var(--blue);margin-bottom:10px}

/* Scrollbar */
::-webkit-scrollbar{width:8px}
::-webkit-scrollbar-track{background:var(--bg)}
::-webkit-scrollbar-thumb{background:var(--bg3);border-radius:4px}

/* Loading */
.spinner{display:inline-block;width:14px;height:14px;border:2px solid var(--border);border-top-color:var(--blue);border-radius:50%;animation:spin .6s linear infinite}
@keyframes spin{to{transform:rotate(360deg)}}

/* Builder */
.builder-section{background:var(--bg2);border:1px solid var(--border);border-radius:8px;padding:18px;margin-bottom:16px}
.builder-section h3{font-size:14px;color:var(--blue);margin-bottom:14px}
.builder-grid{display:grid;grid-template-columns:1fr 1fr;gap:12px}
.builder-grid.three{grid-template-columns:1fr 1fr 1fr}
.builder-full{grid-column:1/-1}
.module-grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(140px,1fr));gap:6px}
.module-toggle{display:flex;align-items:center;gap:6px;padding:5px 8px;border:1px solid var(--border);border-radius:6px;cursor:pointer;font-size:12px;transition:all .12s;user-select:none}
.module-toggle:hover{border-color:var(--blue)}
.module-toggle.on{background:var(--blue);color:#fff;border-color:var(--blue)}
.module-toggle input{display:none}
.cfg-card{background:var(--bg);border:1px solid var(--border);border-radius:6px;padding:12px;margin-bottom:8px}
.cfg-card .cfg-header{display:flex;justify-content:space-between;align-items:center;margin-bottom:8px}
.cfg-card .cfg-header span{font-size:13px;font-weight:500;color:var(--white)}
.btn-remove{background:none;border:none;color:var(--red);cursor:pointer;font-size:16px;padding:0 4px}
.btn-add{display:inline-flex;align-items:center;gap:4px;padding:6px 14px;border-radius:6px;border:1px dashed var(--border);background:transparent;color:var(--blue);cursor:pointer;font-size:12px;transition:all .12s}
.btn-add:hover{border-color:var(--blue);background:var(--bg3)}
textarea{background:var(--bg);border:1px solid var(--border);color:var(--text);border-radius:6px;padding:8px 10px;font-family:"Fira Code",monospace;font-size:12px;width:100%;resize:vertical;min-height:100px}
textarea:focus{outline:none;border-color:var(--blue)}
.builder-tabs{display:flex;gap:0;margin-bottom:16px}
.builder-tab{padding:8px 20px;cursor:pointer;font-size:13px;color:var(--text2);border-bottom:2px solid transparent;transition:all .12s;user-select:none}
.builder-tab:hover{color:var(--text)}
.builder-tab.active{color:var(--blue);border-bottom-color:var(--blue)}
.builder-pane{display:none}
.builder-pane.active{display:block}

@media(max-width:900px){.chart-grid{grid-template-columns:1fr 1fr}.sidebar{width:180px;min-width:180px}}
@media(max-width:600px){.chart-grid{grid-template-columns:1fr}.shell{flex-direction:column}.sidebar{width:100%;min-width:0;flex-direction:row;height:auto;overflow-x:auto}.study-list{display:flex;gap:4px;margin:0}}
</style>
</head>
<body>
<div class="shell">
  <div class="sidebar">
    <div class="logo">&#129516; SkiBiDy</div>
    <div class="logo-sub">Skin BioDynaMo</div>
    <div class="study-section" id="engineSection">
      <div class="study-section-header" onclick="this.parentElement.classList.toggle('collapsed')">
        <span>Engine Studies</span>
        <button class="study-section-toggle">&#9660;</button>
      </div>
      <ul class="study-list" id="engineList"></ul>
    </div>
    <div class="study-section" id="userSection">
      <div class="study-section-header" onclick="this.parentElement.classList.toggle('collapsed')">
        <span>User Studies</span>
        <button class="study-section-toggle">&#9660;</button>
      </div>
      <ul class="study-list" id="userList"></ul>
      <div class="user-dir-label" id="userDirLabel"></div>
    </div>
  </div>
  <div class="main">
    <div class="tabs" id="tabBar">
      <div class="tab active" data-tab="overview">Overview</div>
      <div class="tab" data-tab="metrics">Metrics</div>
      <div class="tab" data-tab="compare">Compare</div>
      <div class="tab" data-tab="validation">Validation</div>
      <div class="tab" data-tab="experiments">Experiments</div>
      <div class="tab" data-tab="config">Config</div>
      <div class="tab" data-tab="builder">Study Builder</div>
    </div>

    <div class="tab-panel active" id="p-overview"></div>
    <div class="tab-panel" id="p-metrics"></div>
    <div class="tab-panel" id="p-compare"></div>
    <div class="tab-panel" id="p-validation"></div>
    <div class="tab-panel" id="p-experiments"></div>
    <div class="tab-panel" id="p-config"></div>
    <div class="tab-panel" id="p-builder"></div>
  </div>
</div>
<div class="action-panel" id="actionPanel">
  <div class="ap-header" onclick="toggleActionPanel()">
    <span class="ap-title">Actions</span>
    <button class="ap-toggle" title="Collapse">&#9660;</button>
  </div>
  <div class="ap-body">
    <div class="ap-row">
      <select id="runCount" style="width:56px">
        <option value="1">1</option><option value="3">3</option>
        <option value="5">5</option><option value="10" selected>10</option>
        <option value="20">20</option>
      </select>
      <input id="runSkin" type="text" placeholder="skin" style="flex:1" />
    </div>
    <div class="ap-row">
      <input id="runTreatment" type="text" placeholder="treatment (optional)" style="flex:1" />
    </div>
    <div class="ap-row">
      <button class="btn" onclick="doRun()" id="runBtn">&#9654; Run</button>
      <button class="btn" onclick="doBuild()">&#128296; Build</button>
      <button class="btn" onclick="doParaview()">&#9654; PV</button>
    </div>
    <div id="runStatus" style="font-size:11px;color:var(--text2);display:none;margin-top:6px"></div>
  </div>
</div>
<div class="toast" id="toast"></div>

<script>
// State
let S = {study:'', studies:[], engineStudies:[], userStudies:[], results:[], csvCache:{}, charts:[], activeTab:'overview'};
const PALETTE = ['#58a6ff','#f85149','#3fb950','#d29922','#bc8cff','#f778ba','#79c0ff','#ffa657'];
const KEY_COLS = ['wound_closure_pct','mean_infl_wound','n_neutrophils','n_macrophages',
  'mean_collagen_wound','mean_mmp_wound','mean_tnf_alpha_wound','mean_il6_wound',
  'mean_cartilage_wound','mean_bone_wound','mean_scab_wound','mean_o2_wound'];

// Helpers
const $ = s => document.querySelector(s);
const $$ = s => document.querySelectorAll(s);
const api = (path,qs={}) => {
  let u = '/api/'+path+'?'+new URLSearchParams(qs);
  return fetch(u).then(r=>r.json());
};
function toast(msg){let t=$('#toast');t.textContent=msg;t.classList.add('show');setTimeout(()=>t.classList.remove('show'),2500)}
function prettify(s){return s.replace(/^mean_/,'').replace(/_wound$/,'').replace(/_/g,' ')}

// Tabs
$('#tabBar').addEventListener('click', e=>{
  let tab = e.target.dataset.tab;
  if(!tab) return;
  $$('.tab').forEach(t=>t.classList.toggle('active',t.dataset.tab===tab));
  $$('.tab-panel').forEach(p=>p.classList.toggle('active',p.id==='p-'+tab));
  S.activeTab = tab;
  renderTab(tab);
});

// Study list
function makeStudyLi(s, showDup){
  let li = document.createElement('li');
  li.dataset.study = s.name;
  let span = document.createElement('span');
  span.className = 's-name';
  span.textContent = s.name.replace(/-/g,' ');
  if(s.description){span.title = s.description}
  li.appendChild(span);
  if(showDup){
    let btn = document.createElement('button');
    btn.className = 'btn-dup';
    btn.title = 'Duplicate to user studies';
    btn.textContent = '\u2398';
    btn.onclick = (e) => {e.stopPropagation(); doDuplicate(s.name)};
    li.appendChild(btn);
  }
  li.onclick = () => selectStudy(s.name);
  return li;
}

async function refreshStudyList(){
  let data = await api('studies');
  S.engineStudies = data.engine || [];
  S.userStudies = data.user || [];
  S.studies = [...S.engineStudies, ...S.userStudies].map(s=>s.name);
  let eul = $('#engineList');
  let uul = $('#userList');
  eul.innerHTML = '';
  uul.innerHTML = '';
  S.engineStudies.forEach(s => eul.appendChild(makeStudyLi(s, true)));
  S.userStudies.forEach(s => uul.appendChild(makeStudyLi(s, false)));
  // Show/hide user section
  if(!S.userStudies.length){
    let lbl = $('#userDirLabel');
    if(lbl) lbl.textContent = 'No user studies yet';
  }
  // Highlight active
  $$('.study-list li').forEach(li=>li.classList.toggle('active',li.dataset.study===S.study));
  // Show user dir path
  api('user-studies-dir').then(r=>{
    let lbl=$('#userDirLabel');
    if(lbl && r.path) lbl.textContent = r.path.replace(os_home(),'~');
  });
}

function os_home(){
  // Best effort home dir detection for display
  let parts = (USER_STUDIES_DIR||'').split('/');
  return parts.length>2 ? '/'+parts[1]+'/'+parts[2] : '';
}
let USER_STUDIES_DIR = '';

async function init(){
  let data = await api('studies');
  S.engineStudies = data.engine || [];
  S.userStudies = data.user || [];
  S.studies = [...S.engineStudies, ...S.userStudies].map(s=>s.name);
  let eul = $('#engineList');
  let uul = $('#userList');
  S.engineStudies.forEach(s => eul.appendChild(makeStudyLi(s, true)));
  S.userStudies.forEach(s => uul.appendChild(makeStudyLi(s, false)));
  // Get user studies dir for display
  api('user-studies-dir').then(r=>{
    USER_STUDIES_DIR = r.path||'';
    let lbl=$('#userDirLabel');
    if(lbl) lbl.textContent = r.path||'';
  });
  // Auto-select study from URL query param (?study=wound)
  let params = new URLSearchParams(window.location.search);
  let target = params.get('study');
  if(target && S.studies.includes(target)){
    selectStudy(target);
  } else if(S.studies.length){
    selectStudy(S.studies[0]);
  }
}

async function selectStudy(name){
  S.study = name;
  $$('.study-list li').forEach(li=>li.classList.toggle('active', li.dataset.study===name));
  S.results = await api('results',{study:name});
  S.csvCache = {};
  renderTab(S.activeTab);
}

async function doDuplicate(sourceName){
  let newName = prompt('New study name (fork of ' + sourceName + '):', sourceName + '-custom');
  if(!newName) return;
  let r = await postApi('duplicate-study', {source: sourceName, name: newName});
  if(r.ok){
    toast('Duplicated "' + sourceName + '" as "' + r.study + '"');
    await refreshStudyList();
    selectStudy(r.study);
  } else {
    toast(r.error || 'Failed to duplicate');
  }
}

async function getCSV(path){
  if(S.csvCache[path]) return S.csvCache[path];
  let d = await api('csv',{path});
  S.csvCache[path] = d;
  return d;
}

// Destroy old charts
function clearCharts(){S.charts.forEach(c=>c.destroy());S.charts=[]}

function makeChart(canvas, label, xData, yData, color){
  let ctx = canvas.getContext('2d');
  let c = new Chart(ctx, {
    type:'line',
    data:{labels:xData, datasets:[{label, data:yData, borderColor:color||PALETTE[0],
      backgroundColor:'transparent', borderWidth:1.5, pointRadius:0, tension:0.3}]},
    options:{responsive:true,maintainAspectRatio:false,animation:false,
      plugins:{legend:{display:false}},
      scales:{x:{ticks:{color:'#8b949e',maxTicksToLimit:6,font:{size:10}},grid:{color:'#21262d'}},
              y:{ticks:{color:'#8b949e',maxTicksToLimit:5,font:{size:10}},grid:{color:'#21262d'}}}}
  });
  S.charts.push(c);
  return c;
}

function makeMultiChart(canvas, xData, datasets){
  let ctx = canvas.getContext('2d');
  let c = new Chart(ctx, {
    type:'line',
    data:{labels:xData, datasets:datasets.map((ds,i)=>({
      label:ds.label, data:ds.data, borderColor:PALETTE[i%PALETTE.length],
      backgroundColor:'transparent', borderWidth:1.5, pointRadius:0, tension:0.3
    }))},
    options:{responsive:true,maintainAspectRatio:false,animation:false,
      plugins:{legend:{labels:{color:'#c9d1d9',font:{size:11}}}},
      scales:{x:{ticks:{color:'#8b949e',maxTicksToLimit:6,font:{size:10}},grid:{color:'#21262d'}},
              y:{ticks:{color:'#8b949e',maxTicksToLimit:5,font:{size:10}},grid:{color:'#21262d'}}}}
  });
  S.charts.push(c);
}

// ---------------------------------------------------------------------------
// Tab renderers
// ---------------------------------------------------------------------------

async function renderTab(tab){
  clearCharts();
  if(tab==='overview') renderOverview();
  else if(tab==='metrics') renderMetrics();
  else if(tab==='compare') renderCompare();
  else if(tab==='validation') renderValidation();
  else if(tab==='experiments') renderExperiments();
  else if(tab==='config') renderConfig();
  else if(tab==='builder') renderBuilder();
}

async function renderOverview(){
  let el = $('#p-overview');
  if(!S.results.length){el.innerHTML='<div class="info-box">No results yet. Run a simulation first.</div>';return}
  let r = S.results[0];
  let data = await getCSV(r.csv);
  if(data.error){el.innerHTML='<div class="info-box">Error loading CSV</div>';return}

  let days = data.time_days||[];
  let last = days.length ? days[days.length-1] : 0;
  let agents = data.n_agents ? Math.round(data.n_agents[data.n_agents.length-1]) : 0;
  let closure = data.wound_closure_pct ? data.wound_closure_pct[data.wound_closure_pct.length-1].toFixed(1)+'%' : null;
  let cart = data.mean_cartilage_wound ? data.mean_cartilage_wound[data.mean_cartilage_wound.length-1].toFixed(3) : null;

  let html = '<div class="cards">';
  html += `<div class="card"><div class="card-label">Study</div><div class="card-value">${S.study.replace(/_/g,' ')}</div><div class="card-sub">${S.results.length} result(s)</div></div>`;
  html += `<div class="card"><div class="card-label">Duration</div><div class="card-value">${last.toFixed(1)}d</div><div class="card-sub">${days.length} steps</div></div>`;
  if(closure) html += `<div class="card"><div class="card-label">Wound Closure</div><div class="card-value">${closure}</div></div>`;
  else if(cart) html += `<div class="card"><div class="card-label">Cartilage</div><div class="card-value">${cart}</div></div>`;
  html += `<div class="card"><div class="card-label">Agents</div><div class="card-value">${agents.toLocaleString()}</div></div>`;
  html += '</div>';

  // Charts
  let cols = KEY_COLS.filter(c=>data[c]&&data[c].some(v=>v!==0));
  cols = cols.slice(0,6);
  if(cols.length){
    html += '<div class="chart-grid">';
    cols.forEach((c,i) => {
      html += `<div class="chart-box"><div class="chart-title">${prettify(c)}</div><canvas id="ov-${i}"></canvas></div>`;
    });
    html += '</div>';
  }
  el.innerHTML = html;
  // Draw charts after DOM update
  requestAnimationFrame(()=>{
    cols.forEach((c,i)=>{
      let cv = document.getElementById('ov-'+i);
      if(cv) makeChart(cv, prettify(c), days, data[c]);
    });
  });
}

async function renderMetrics(){
  let el = $('#p-metrics');
  if(!S.results.length){el.innerHTML='<div class="info-box">No data.</div>';return}

  let html = '<div class="run-bar"><select id="m-run">';
  S.results.forEach((r,i)=>{html+=`<option value="${i}">[${r.type}] ${r.label}</option>`});
  html += '</select></div><div id="m-chips" class="chip-bar"></div><div id="m-charts" class="chart-grid"></div>';
  el.innerHTML = html;

  let selEl = $('#m-run');
  let loadRun = async () => {
    clearCharts();
    let r = S.results[selEl.value];
    let data = await getCSV(r.csv);
    let days = data.time_days||[];
    let skip = new Set(['step','time_h','time_days']);
    let cols = Object.keys(data).filter(k=>!skip.has(k)&&typeof data[k][0]==='number');

    // Chips
    let chipEl = $('#m-chips');
    let active = new Set(cols.filter(c=>KEY_COLS.includes(c)).slice(0,6));
    if(!active.size) cols.slice(0,6).forEach(c=>active.add(c));

    function renderChips(){
      chipEl.innerHTML = '';
      cols.forEach(c=>{
        let chip = document.createElement('span');
        chip.className = 'chip'+(active.has(c)?' on':'');
        chip.textContent = prettify(c);
        chip.onclick = ()=>{active.has(c)?active.delete(c):active.add(c);renderChips();drawCharts()};
        chipEl.appendChild(chip);
      });
    }
    function drawCharts(){
      clearCharts();
      let cg = $('#m-charts');
      cg.innerHTML='';
      [...active].forEach((c,i)=>{
        let box = document.createElement('div');box.className='chart-box';
        box.innerHTML=`<div class="chart-title">${prettify(c)}</div><canvas id="mc-${i}"></canvas>`;
        cg.appendChild(box);
      });
      requestAnimationFrame(()=>{
        [...active].forEach((c,i)=>{
          let cv=document.getElementById('mc-'+i);
          if(cv) makeChart(cv, prettify(c), days, data[c]);
        });
      });
    }
    renderChips();drawCharts();
  };
  selEl.onchange = loadRun;
  loadRun();
}

async function renderCompare(){
  let el = $('#p-compare');
  if(S.results.length<2){el.innerHTML='<div class="info-box">Need 2+ results to compare.</div>';return}

  let html = '<div class="field"><label>Select runs (hold Ctrl/Cmd)</label><select id="cmp-runs" multiple size="4">';
  S.results.forEach((r,i)=>{let sel=i<2?' selected':'';html+=`<option value="${i}"${sel}>[${r.type}] ${r.label}</option>`});
  html += '</select></div>';
  html += '<div class="field"><label>Observable</label><select id="cmp-col"></select></div>';
  html += '<div id="cmp-chart-area"></div>';
  el.innerHTML = html;

  let runsEl = $('#cmp-runs'), colEl = $('#cmp-col'), chartArea = $('#cmp-chart-area');

  // Populate column list from first result
  let d0 = await getCSV(S.results[0].csv);
  let skip = new Set(['step','time_h','time_days']);
  let cols = Object.keys(d0).filter(k=>!skip.has(k)&&typeof d0[k][0]==='number');
  cols.forEach(c=>{let o=document.createElement('option');o.value=c;o.textContent=prettify(c);colEl.appendChild(o)});

  async function draw(){
    clearCharts();
    chartArea.innerHTML='';
    let idxs = [...runsEl.selectedOptions].map(o=>+o.value);
    let col = colEl.value;
    if(idxs.length<2||!col) return;

    chartArea.innerHTML='<div class="chart-box" style="grid-column:1/-1"><canvas id="cmp-cv"></canvas></div>';
    let datasets=[], xData=null;
    for(let idx of idxs){
      let d = await getCSV(S.results[idx].csv);
      if(!d[col]) continue;
      if(!xData) xData = d.time_days||[];
      datasets.push({label:S.results[idx].label.slice(0,30), data:d[col]});
    }
    requestAnimationFrame(()=>{
      let cv=document.getElementById('cmp-cv');
      if(cv&&datasets.length) makeMultiChart(cv, xData, datasets);
    });
  }
  runsEl.onchange=draw; colEl.onchange=draw;
  draw();
}

async function renderValidation(){
  let el = $('#p-validation');
  if(!S.results.length){el.innerHTML='<div class="info-box">No results.</div>';return}

  let html = '<div class="run-bar"><select id="v-run">';
  S.results.forEach((r,i)=>{html+=`<option value="${i}">[${r.type}] ${r.label}</option>`});
  html += '</select></div><div id="v-body"><div class="spinner"></div> Validating...</div>';
  el.innerHTML = html;

  let selEl = $('#v-run');
  async function validate(){
    let r = S.results[selEl.value];
    let body = $('#v-body');
    body.innerHTML = '<div class="spinner"></div> Validating...';
    let val = await api('validate',{path:r.csv});
    if(val.error){body.innerHTML=`<div class="info-box" style="color:var(--red)">${val.error}</div>`;return}
    let entries = Object.entries(val);
    if(!entries.length){body.innerHTML='<div class="info-box">No validatable modules detected.</div>';return}

    let avg = entries.reduce((s,[_,v])=>s+v,0)/entries.length;
    let passed = entries.filter(([_,v])=>v<0.10).length;
    let worst = entries.reduce((a,b)=>b[1]>a[1]?b:a);

    let h = '<div class="cards">';
    h += `<div class="card"><div class="card-label">Mean RMSE</div><div class="card-value">${(avg*100).toFixed(1)}%</div><div class="card-sub">across ${entries.length} observables</div></div>`;
    h += `<div class="card"><div class="card-label">Passing</div><div class="card-value">${passed}/${entries.length}</div><div class="card-sub">&lt; 10% RMSE</div></div>`;
    h += `<div class="card"><div class="card-label">Worst</div><div class="card-value">${worst[0]}</div><div class="card-sub">${(worst[1]*100).toFixed(1)}%</div></div>`;
    h += '</div>';

    h += '<div class="rmse-table">';
    entries.forEach(([name,val])=>{
      let cls = val<0.10?'rmse-good':val<0.20?'rmse-ok':'rmse-bad';
      let pcls = val<0.10?'pill-green':val<0.20?'pill-yellow':'pill-red';
      let plbl = val<0.10?'PASS':val<0.20?'OK':'HIGH';
      h += `<div class="rmse-row"><span>${name}</span><span><span class="${cls}">${(val*100).toFixed(1)}%</span> <span class="pill ${pcls}">${plbl}</span></span></div>`;
    });
    h += '</div>';
    body.innerHTML = h;
  }
  selEl.onchange = validate;
  validate();
}

async function renderExperiments(){
  let el = $('#p-experiments');
  let exps = await api('experiments',{study:S.study});
  if(!exps.length){el.innerHTML=`<div class="info-box">No experiments for ${S.study}.</div>`;return}
  let html = '';
  exps.forEach((e,i)=>{
    html += `<div class="exp-header" onclick="this.nextElementSibling.classList.toggle('open')">${e.name} <span style="color:var(--text2)">&rsaquo;</span></div>`;
    html += `<div class="exp-body"><div class="toml-block">${esc(e.content)}</div></div>`;
  });
  el.innerHTML = html;
}

async function renderConfig(){
  let el = $('#p-config');
  let [treats, mods] = await Promise.all([api('treatments',{study:S.study}), api('modules')]);
  let html = '<div class="config-grid"><div class="config-col"><h3>Treatments</h3>';
  treats.forEach(t=>{
    html += `<div class="exp-header" onclick="this.nextElementSibling.classList.toggle('open')">${t.name} <span style="font-size:11px;color:var(--text2)">${t.scope}</span></div>`;
    html += `<div class="exp-body"><div class="toml-block">${esc(t.content)}</div></div>`;
  });
  html += '</div><div class="config-col"><h3>Modules</h3>';
  mods.forEach(m=>{
    html += `<div class="exp-header" onclick="this.nextElementSibling.classList.toggle('open')">${m.name}</div>`;
    html += `<div class="exp-body"><div class="toml-block">${esc(m.content)}</div></div>`;
  });
  html += '</div></div>';
  el.innerHTML = html;
}

function esc(s){return s.replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;')}

function postApi(path, body){
  return fetch('/api/'+path, {method:'POST', headers:{'Content-Type':'application/json'}, body:JSON.stringify(body)}).then(r=>r.json());
}

// Builder state
let B = {profiles:[], modules:[], selectedMods:new Set(), expConfigs:[], builderSubTab:'study'};

async function renderBuilder(){
  let el = $('#p-builder');
  // Load profiles and modules if not cached
  if(!B.profiles.length) B.profiles = await api('profiles');
  if(!B.modules.length){
    let mods = await api('modules');
    B.modules = mods.map(m=>m.name);
  }

  let html = '';
  // Sub-tabs: New Study, New Experiment, New Treatment
  html += '<div class="builder-tabs">';
  html += `<div class="builder-tab ${B.builderSubTab==='study'?'active':''}" onclick="switchBuilderTab('study')">New Study</div>`;
  html += `<div class="builder-tab ${B.builderSubTab==='experiment'?'active':''}" onclick="switchBuilderTab('experiment')">New Experiment</div>`;
  html += `<div class="builder-tab ${B.builderSubTab==='treatment'?'active':''}" onclick="switchBuilderTab('treatment')">New Treatment</div>`;
  html += '</div>';

  // -- New Study pane --
  html += `<div class="builder-pane ${B.builderSubTab==='study'?'active':''}" id="bp-study">`;
  html += '<div class="builder-section"><h3>Create Study</h3>';
  html += '<div class="builder-grid">';
  html += '<div class="field"><label>Study name</label><input type="text" id="b-study-name" placeholder="e.g. venous-ulcer" /></div>';
  html += '<div class="field"><label>Description</label><input type="text" id="b-study-desc" placeholder="Short description" /></div>';
  html += '<div class="field"><label>Duration (days)</label><input type="text" id="b-study-dur" value="30" /></div>';
  html += `<div class="field"><label>Skin profile</label><select id="b-study-profile"><option value="">(none)</option>${B.profiles.map(p=>'<option value="'+p+'">'+p+'</option>').join('')}</select></div>`;
  html += '<div class="field"><label>Steps (auto)</label><input type="text" id="b-study-steps" disabled /></div>';
  html += '</div>';
  html += '<div class="field" style="margin-top:12px"><label>Enable modules</label>';
  html += '<div class="module-grid" id="b-mod-grid">';
  B.modules.forEach(m=>{
    let on = B.selectedMods.has(m);
    html += `<div class="module-toggle ${on?'on':''}" onclick="toggleMod(this,'${m}')">${m}</div>`;
  });
  html += '</div></div>';
  html += '<div style="margin-top:16px;display:flex;gap:8px">';
  html += '<button class="btn" onclick="doCreateStudy()" style="flex:none;padding:8px 24px">Create Study</button>';
  html += '</div>';
  html += '</div></div>';

  // -- New Experiment pane --
  html += `<div class="builder-pane ${B.builderSubTab==='experiment'?'active':''}" id="bp-experiment">`;
  html += '<div class="builder-section"><h3>Create Experiment</h3>';
  html += '<p style="font-size:12px;color:var(--text2);margin-bottom:12px">Creates an experiment TOML in the currently selected study: <strong style="color:var(--blue)">' + S.study + '</strong></p>';
  html += '<div class="builder-grid">';
  html += '<div class="field"><label>Experiment name</label><input type="text" id="b-exp-name" placeholder="e.g. dose_response" /></div>';
  html += `<div class="field"><label>Profile override</label><select id="b-exp-profile"><option value="">(study default)</option>${B.profiles.map(p=>'<option value="'+p+'">'+p+'</option>').join('')}</select></div>`;
  html += '<div class="field"><label>Description</label><input type="text" id="b-exp-desc" placeholder="Brief description" /></div>';
  html += '<div class="field"><label>Runs per config</label><input type="text" id="b-exp-runs" value="5" /></div>';
  html += '</div>';

  // Configs
  html += '<div style="margin-top:14px"><label style="font-size:12px;color:var(--text2);display:block;margin-bottom:8px">Configurations</label>';
  html += '<div id="b-exp-configs">';
  B.expConfigs.forEach((cfg,i) => {
    html += buildCfgCard(cfg, i);
  });
  html += '</div>';
  html += '<button class="btn-add" onclick="addExpConfig()">+ Add configuration</button>';
  html += '</div>';

  html += '<div style="margin-top:16px"><button class="btn" onclick="doCreateExperiment()" style="padding:8px 24px">Create Experiment</button></div>';
  html += '</div></div>';

  // -- New Treatment pane --
  html += `<div class="builder-pane ${B.builderSubTab==='treatment'?'active':''}" id="bp-treatment">`;
  html += '<div class="builder-section"><h3>Create Treatment</h3>';
  html += '<p style="font-size:12px;color:var(--text2);margin-bottom:12px">Creates a treatment TOML in the currently selected study: <strong style="color:var(--blue)">' + S.study + '</strong></p>';
  html += '<div class="builder-grid">';
  html += '<div class="field"><label>Treatment name</label><input type="text" id="b-treat-name" placeholder="e.g. collagenase" /></div>';
  html += '<div class="field"><label>Based on</label><select id="b-treat-template"><option value="">(blank)</option></select></div>';
  html += '</div>';
  html += '<div class="field" style="margin-top:8px"><label>TOML content</label>';
  html += '<textarea id="b-treat-content" rows="12" placeholder="# Treatment parameters\n[skin.diabetic]\nprolif_factor = 0.8"></textarea>';
  html += '</div>';
  html += '<div style="margin-top:12px"><button class="btn" onclick="doCreateTreatment()" style="padding:8px 24px">Create Treatment</button></div>';
  html += '</div></div>';

  el.innerHTML = html;

  // Auto-calc steps from duration
  let durEl = document.getElementById('b-study-dur');
  let stepsEl = document.getElementById('b-study-steps');
  if(durEl && stepsEl){
    let calc = ()=>{let d=parseFloat(durEl.value)||30; stepsEl.value=Math.round(d*24/0.1)+' steps'};
    durEl.oninput = calc;
    calc();
  }

  // Load treatment templates for "based on" dropdown
  if(S.study){
    let treats = await api('treatments',{study:S.study});
    let sel = document.getElementById('b-treat-template');
    if(sel){
      treats.forEach(t=>{
        let o = document.createElement('option');
        o.value = t.name;
        o.textContent = t.name + ' (' + t.scope + ')';
        o.dataset.content = t.content;
        sel.appendChild(o);
      });
      sel.onchange = ()=>{
        let opt = sel.selectedOptions[0];
        let ta = document.getElementById('b-treat-content');
        if(opt && opt.dataset.content) ta.value = opt.dataset.content;
      };
    }
  }
}

function switchBuilderTab(tab){
  B.builderSubTab = tab;
  document.querySelectorAll('.builder-tab').forEach(t=>t.classList.toggle('active',t.textContent.toLowerCase().includes(tab)));
  document.querySelectorAll('.builder-pane').forEach(p=>p.classList.toggle('active',p.id==='bp-'+tab));
}

function toggleMod(el, name){
  if(B.selectedMods.has(name)){B.selectedMods.delete(name);el.classList.remove('on')}
  else{B.selectedMods.add(name);el.classList.add('on')}
}

function buildCfgCard(cfg, idx){
  let h = `<div class="cfg-card" data-idx="${idx}">`;
  h += '<div class="cfg-header"><span>Config ' + (idx+1) + '</span><button class="btn-remove" onclick="removeExpConfig('+idx+')">&times;</button></div>';
  h += '<div class="builder-grid">';
  h += `<div class="field"><label>Label</label><input type="text" class="cfg-label" value="${esc(cfg.label||'')}" placeholder="e.g. low dose" onchange="B.expConfigs[${idx}].label=this.value" /></div>`;
  h += `<div class="field"><label>Treatments (comma-separated)</label><input type="text" class="cfg-treats" value="${(cfg.treatments||[]).join(', ')}" placeholder="e.g. anti_tnf, methotrexate" onchange="B.expConfigs[${idx}].treatments=this.value.split(',').map(s=>s.trim()).filter(Boolean)" /></div>`;
  h += `<div class="field builder-full"><label>Overrides (key = value, one per line)</label><textarea class="cfg-overrides" rows="3" placeholder='"skin.diabetic.prolif_factor" = 0.5\n"skin.immune.macrophage_m1_duration_h" = 72' onchange="parseOverrides(${idx},this.value)">${formatOverrides(cfg.overrides||{})}</textarea></div>`;
  h += '</div></div>';
  return h;
}

function formatOverrides(obj){
  return Object.entries(obj).map(([k,v])=>'"'+k+'" = '+JSON.stringify(v)).join('\n');
}

function parseOverrides(idx, text){
  let obj = {};
  text.split('\n').forEach(line=>{
    line = line.trim();
    if(!line) return;
    let m = line.match(/^"([^"]+)"\s*=\s*(.+)$/);
    if(m){
      let val = m[2].trim();
      try{val = JSON.parse(val)}catch(e){}
      obj[m[1]] = val;
    }
  });
  B.expConfigs[idx].overrides = obj;
}

function addExpConfig(){
  B.expConfigs.push({label:'', treatments:[], overrides:{}});
  let container = document.getElementById('b-exp-configs');
  container.innerHTML = '';
  B.expConfigs.forEach((cfg,i)=>{container.innerHTML += buildCfgCard(cfg,i)});
}

function removeExpConfig(idx){
  B.expConfigs.splice(idx, 1);
  let container = document.getElementById('b-exp-configs');
  container.innerHTML = '';
  B.expConfigs.forEach((cfg,i)=>{container.innerHTML += buildCfgCard(cfg,i)});
}

async function doCreateStudy(){
  let name = document.getElementById('b-study-name').value.trim();
  let dur = parseInt(document.getElementById('b-study-dur').value) || 30;
  let profile = document.getElementById('b-study-profile').value;
  let mods = [...B.selectedMods];
  let desc = document.getElementById('b-study-desc')?.value.trim()||'';
  if(!name){toast('Enter a study name');return}
  let r = await postApi('create-study', {name, duration:dur, profile, modules:mods, description:desc});
  if(r.ok){
    toast('Study "'+r.study+'" created in user studies');
    await refreshStudyList();
    selectStudy(r.study);
  } else {
    toast(r.error||'Failed');
  }
}

async function doCreateExperiment(){
  let name = document.getElementById('b-exp-name').value.trim();
  let desc = document.getElementById('b-exp-desc').value.trim();
  let profile = document.getElementById('b-exp-profile').value;
  let runs = parseInt(document.getElementById('b-exp-runs').value) || 5;
  if(!name){toast('Enter an experiment name');return}
  if(!B.expConfigs.length){toast('Add at least one configuration');return}
  let r = await postApi('create-experiment', {study:S.study, name, description:desc, profile, runs, configs:B.expConfigs});
  if(r.ok){
    toast('Experiment "'+r.name+'" created');
    B.expConfigs = [];
  } else {
    toast(r.error||'Failed');
  }
}

async function doCreateTreatment(){
  let name = document.getElementById('b-treat-name').value.trim();
  let content = document.getElementById('b-treat-content').value;
  if(!name){toast('Enter a treatment name');return}
  if(!content.trim()){toast('Enter treatment TOML content');return}
  let r = await postApi('create-treatment', {study:S.study, name, content});
  if(r.ok){
    toast('Treatment "'+r.name+'" created');
  } else {
    toast(r.error||'Failed');
  }
}

// Action panel toggle
function toggleActionPanel(){
  let p = document.getElementById('actionPanel');
  p.classList.toggle('collapsed');
}

// Actions
async function doBuild(){
  toast('Building...');
  let r = await api('build');
  toast(r.ok?'Build succeeded':'Build failed');
}

async function doParaview(){
  if(!S.results.length){toast('No results');return}
  let r = await api('paraview',{dir:S.results[0].dir});
  toast(r.ok?'ParaView launched':'ParaView not found');
}

async function doRun(){
  let n = $('#runCount').value;
  let skin = $('#runSkin').value.trim();
  let treatment = $('#runTreatment').value.trim();
  let params = {study:S.study, n};
  if(skin) params.skin = skin;
  if(treatment) params.treatment = treatment;
  let r = await api('run', params);
  if(!r.ok){toast(r.error||'Failed to start');return}
  toast(`Running ${n} sim(s)...`);
  $('#runBtn').disabled = true;
  $('#runBtn').textContent = '\u23f3 Running...';
  let status = $('#runStatus');
  status.style.display = 'block';
  // Poll status
  let poll = setInterval(async()=>{
    let s = await api('run/status');
    let lines = s.log.trim().split('\n');
    let last = lines.filter(l=>l.trim()).slice(-3).join('\n');
    status.textContent = last;
    if(!s.running){
      clearInterval(poll);
      $('#runBtn').disabled = false;
      $('#runBtn').textContent = '\u25b6 Run Batch';
      toast(s.returncode===0?'Batch complete':'Batch finished with errors');
      // Refresh results
      S.results = await api('results',{study:S.study});
      S.csvCache = {};
      renderTab(S.activeTab);
      setTimeout(()=>{status.style.display='none'},5000);
    }
  }, 2000);
}

init();
</script>
</body>
</html>
"""

# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------

def main():
    port = PORT
    study = ""
    args = sys.argv[1:]
    i = 0
    while i < len(args):
        if args[i] == "--port" and i + 1 < len(args):
            port = int(args[i + 1])
            i += 2
        elif args[i] == "--study" and i + 1 < len(args):
            study = args[i + 1]
            i += 2
        elif args[i].endswith(".skibidy"):
            # Double-clicked .skibidy file: extract study name from path
            proj_path = os.path.abspath(args[i])
            study = os.path.basename(os.path.dirname(proj_path))
            i += 1
        else:
            i += 1

    server = http.server.HTTPServer(("0.0.0.0", port), DashHandler)
    url = f"http://localhost:{port}"
    if study:
        url += f"?study={study}"
    print(f"SkiBiDy Dashboard running at {url}")

    # Open browser
    threading.Timer(0.5, lambda: webbrowser.open(url)).start()

    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nShutting down.")
        server.server_close()


if __name__ == "__main__":
    main()
