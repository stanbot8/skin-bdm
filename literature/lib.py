"""Shared utilities for validation scripts.

Provides CSV loading, interpolation, normalization, error metrics,
module-level validate/plot/print functions used by validate_all.py
and the individual compare_*.py wrappers.
"""

import csv
import math
import os


# ---------------------------------------------------------------------------
# Low-level helpers
# ---------------------------------------------------------------------------

def load_csv(path):
    """Load a CSV file into a dict of lists, skipping comment lines."""
    with open(path) as f:
        reader = csv.DictReader(
            (row for row in f if not row.startswith("#"))
        )
        rows = list(reader)
    data = {}
    for key in rows[0]:
        try:
            data[key] = [float(r[key]) for r in rows]
        except ValueError:
            data[key] = [r[key] for r in rows]
    return data


def interpolate(x_ref, y_ref, x_query):
    """Linear interpolation of (x_ref, y_ref) at x_query points."""
    result = []
    for xq in x_query:
        if xq <= x_ref[0]:
            result.append(y_ref[0])
        elif xq >= x_ref[-1]:
            result.append(y_ref[-1])
        else:
            for i in range(len(x_ref) - 1):
                if x_ref[i] <= xq <= x_ref[i + 1]:
                    t = (xq - x_ref[i]) / (x_ref[i + 1] - x_ref[i])
                    result.append(y_ref[i] + t * (y_ref[i + 1] - y_ref[i]))
                    break
    return result


def peak_normalize(values):
    """Normalize a list by its maximum value. Returns (normalized, peak)."""
    peak = max(values)
    if peak > 0:
        return [v / peak for v in values], peak
    return values, 0


def end_normalize(values):
    """Normalize a list by its final value. Returns (normalized, final)."""
    final = values[-1] if values else 0
    if final > 0:
        return [v / final for v in values], final
    return values, 0


def compute_rmse(sim, ref):
    """Compute root-mean-square error between two lists."""
    if not sim:
        return 0.0
    errors = [s - r for s, r in zip(sim, ref)]
    return math.sqrt(sum(e ** 2 for e in errors) / len(errors))


def phase_rmse(sim_days, sim_vals, ref_at_sim, day_start, day_end):
    """Compute RMSE over a day range."""
    indices = [i for i, d in enumerate(sim_days)
               if day_start <= d <= day_end]
    if not indices:
        return 0.0
    errors = [sim_vals[i] - ref_at_sim[i] for i in indices]
    return math.sqrt(sum(e ** 2 for e in errors) / len(errors))


def surface_fraction(n_cells, cell_diameter=4.0, packing=0.5):
    """Estimate surface-shell fraction for a sphere of n_cells."""
    if n_cells <= 1:
        return 1.0
    R = (cell_diameter / 2.0) * (n_cells / packing) ** (1.0 / 3.0)
    ratio = max(0, 1 - cell_diameter / R)
    return 1.0 - ratio ** 3


# ---------------------------------------------------------------------------
# Reference data resolution (consensus CSVs live inside each module)
# ---------------------------------------------------------------------------

_MODULES_DIR = os.path.join(os.path.dirname(__file__), os.pardir, "modules")

_FILE_MODULE_MAP = {
    "closure_kinetics_punch_biopsy.csv": "wound",
    "inflammation_timecourse.csv": "inflammation",
    "immune_cell_kinetics.csv": "immune",
    "myofibroblast_kinetics.csv": "fibroblast",
    "collagen_deposition.csv": "fibroblast",
    "tgfb_kinetics.csv": "fibroblast",
    "fibroblast_kinetics.csv": "fibroblast",
    "fibronectin_kinetics.csv": "fibronectin",
    "mmp_kinetics.csv": "mmp",
    "vegf_kinetics.csv": "angiogenesis",
    "diabetic_closure_kinetics.csv": "diabetic",
    "diabetic_inflammation_timecourse.csv": "diabetic",
    "diabetic_immune_cell_kinetics.csv": "diabetic",
    "tumor_doubling_time.csv": "tumor",
    "tumor_growth_rate.csv": "tumor",
    "tumor_proliferation_index.csv": "tumor",
}


def _ref_path(filename):
    """Resolve a consensus CSV filename to its module data path."""
    module = _FILE_MODULE_MAP[filename]
    return os.path.join(_MODULES_DIR, module, "data", filename)


# ---------------------------------------------------------------------------
# Condition detection (normal vs diabetic)
# ---------------------------------------------------------------------------

def detect_condition(config_path=None):
    """Detect simulation condition from bdm.toml.

    Looks for [skin.diabetic] mode = true in the config file.
    Returns "diabetic" or "normal".
    """
    if config_path is None:
        config_path = os.path.join(os.path.dirname(__file__), os.pardir,
                                   "bdm.toml")
    try:
        with open(config_path) as f:
            in_diabetic_section = False
            for line in f:
                stripped = line.strip()
                if stripped.startswith("["):
                    in_diabetic_section = (stripped == "[skin.diabetic]")
                elif in_diabetic_section and "mode" in stripped:
                    if "=" in stripped:
                        val = stripped.split("=", 1)[1].strip().lower()
                        if val.startswith("true"):
                            return "diabetic"
    except FileNotFoundError:
        pass
    return "normal"


# ---------------------------------------------------------------------------
# Plot style constants
# ---------------------------------------------------------------------------

SIM_COLOR = "#D46664"
REF_COLOR = "#4A90D9"
REF_KW = dict(color=REF_COLOR, linewidth=1.5, linestyle="--",
              marker="o", markersize=4, label="Literature")
SIM_KW = dict(color=SIM_COLOR, linewidth=2, label="Simulation")


# ---------------------------------------------------------------------------
# Module detection
# ---------------------------------------------------------------------------

def detect_modules(sim):
    """Return (has_wound, has_fibroblast, has_tumor, has_microenv) booleans."""
    has_wound = ("wound_closure_pct" in sim
                 and max(sim["wound_closure_pct"]) > 0)
    has_fibroblast = (has_wound and "n_myofibroblasts" in sim
                      and max(sim["n_myofibroblasts"]) > 0)
    has_tumor = ("n_tumor_cells" in sim
                 and (max(sim["n_tumor_cells"]) > 0
                      or ("tumor_field_cells" in sim
                          and max(sim["tumor_field_cells"]) > 0)))
    has_microenv = (has_wound
                    and "mean_tgfb_wound" in sim
                    and max(sim["mean_tgfb_wound"]) > 0)
    return has_wound, has_fibroblast, has_tumor, has_microenv


# ---------------------------------------------------------------------------
# Validate functions (pure computation, no I/O)
# ---------------------------------------------------------------------------

def validate_wound(sim, sim_days, condition="normal"):
    """Compute wound validation metrics. Returns result dict.

    When condition="diabetic", loads diabetic-specific reference curves
    derived from literature scaling factors (Mirza & Koh 2011, Galiano 2004,
    Michaels 2007, Wetzler 2000, Khanna 2010).
    """
    if condition == "diabetic":
        ref_closure = load_csv(_ref_path("diabetic_closure_kinetics.csv"))
        ref_infl = load_csv(_ref_path("diabetic_inflammation_timecourse.csv"))
        ref_immune = load_csv(_ref_path("diabetic_immune_cell_kinetics.csv"))
    else:
        ref_closure = load_csv(_ref_path("closure_kinetics_punch_biopsy.csv"))
        ref_infl = load_csv(_ref_path("inflammation_timecourse.csv"))
        ref_immune = load_csv(_ref_path("immune_cell_kinetics.csv"))

    sim_closure = sim["wound_closure_pct"]
    ref_closure_at_sim = interpolate(
        ref_closure["day"], ref_closure["closure_pct"], sim_days)
    closure_rmse = compute_rmse(sim_closure, ref_closure_at_sim)
    closure_max = max(abs(s - r) for s, r in zip(sim_closure, ref_closure_at_sim))

    infl_rmse = phase_rmse(sim_days, sim_closure, ref_closure_at_sim, 0, 3)
    prolif_rmse = phase_rmse(sim_days, sim_closure, ref_closure_at_sim, 3, 14)
    remod_rmse = phase_rmse(sim_days, sim_closure, ref_closure_at_sim, 14, 28)

    sim_infl, infl_peak = peak_normalize(sim["mean_infl_wound"])
    ref_infl_at_sim = interpolate(
        ref_infl["day"], ref_infl["inflammation_normalized"], sim_days)
    inflammation_rmse = compute_rmse(sim_infl, ref_infl_at_sim)

    sim_neut, neut_peak = peak_normalize(sim["n_neutrophils"])
    ref_neut_at_sim = interpolate(
        ref_immune["day"], ref_immune["neutrophils_normalized"], sim_days)
    neut_rmse = compute_rmse(sim_neut, ref_neut_at_sim)

    sim_mac, mac_peak = peak_normalize(sim["n_macrophages"])
    ref_mac_at_sim = interpolate(
        ref_immune["day"], ref_immune["macrophages_normalized"], sim_days)
    mac_rmse = compute_rmse(sim_mac, ref_mac_at_sim)

    return dict(
        # Series for plotting
        condition=condition,
        sim_closure=sim_closure, ref_closure=ref_closure,
        ref_closure_at_sim=ref_closure_at_sim,
        sim_infl=sim_infl, ref_infl=ref_infl,
        ref_infl_at_sim=ref_infl_at_sim, infl_peak=infl_peak,
        sim_neut=sim_neut, ref_immune=ref_immune,
        ref_neut_at_sim=ref_neut_at_sim, neut_peak=neut_peak,
        sim_mac=sim_mac, ref_mac_at_sim=ref_mac_at_sim, mac_peak=mac_peak,
        # Metrics
        closure_rmse=closure_rmse, closure_max=closure_max,
        infl_rmse=infl_rmse, prolif_rmse=prolif_rmse, remod_rmse=remod_rmse,
        inflammation_rmse=inflammation_rmse,
        neut_rmse=neut_rmse, mac_rmse=mac_rmse,
    )


def validate_fibroblast(sim, sim_days):
    """Compute fibroblast/collagen validation metrics. Returns result dict."""
    ref_myofib = load_csv(_ref_path("myofibroblast_kinetics.csv"))
    ref_collagen = load_csv(_ref_path("collagen_deposition.csv"))
    ref_fibro = load_csv(_ref_path("fibroblast_kinetics.csv"))

    sim_myofib, myofib_peak = peak_normalize(sim["n_myofibroblasts"])
    ref_myofib_at_sim = interpolate(
        ref_myofib["day"], ref_myofib["myofibroblasts_normalized"], sim_days)
    myofib_rmse = compute_rmse(sim_myofib, ref_myofib_at_sim)

    sim_collagen, collagen_final = end_normalize(sim["mean_collagen_wound"])
    ref_collagen_at_sim = interpolate(
        ref_collagen["day"], ref_collagen["collagen_normalized"], sim_days)
    collagen_rmse = compute_rmse(sim_collagen, ref_collagen_at_sim)

    sim_fibro, fibro_peak = peak_normalize(sim["n_fibroblasts"])
    ref_fibro_at_sim = interpolate(
        ref_fibro["day"], ref_fibro["fibroblasts_normalized"], sim_days)
    fibro_rmse = compute_rmse(sim_fibro, ref_fibro_at_sim)

    return dict(
        sim_myofib=sim_myofib, ref_myofib=ref_myofib,
        ref_myofib_at_sim=ref_myofib_at_sim, myofib_peak=myofib_peak,
        sim_collagen=sim_collagen, ref_collagen=ref_collagen,
        ref_collagen_at_sim=ref_collagen_at_sim, collagen_final=collagen_final,
        sim_fibro=sim_fibro, ref_fibro=ref_fibro,
        ref_fibro_at_sim=ref_fibro_at_sim, fibro_peak=fibro_peak,
        myofib_rmse=myofib_rmse, collagen_rmse=collagen_rmse,
        fibro_rmse=fibro_rmse,
    )


def validate_tumor(sim, sim_days):
    """Compute tumor validation metrics. Returns result dict."""
    ref_td = load_csv(_ref_path("tumor_doubling_time.csv"))
    ref_ki67 = load_csv(_ref_path("tumor_proliferation_index.csv"))

    bcc_doubling_days = 148.0
    for i, t in enumerate(ref_td["tumor_type"]):
        if t == "bcc_mean":
            bcc_doubling_days = ref_td["doubling_time_days"][i]
            break

    bcc_ki67_pct = 27.4
    for i, t in enumerate(ref_ki67["tumor_type"]):
        if t == "bcc_all":
            bcc_ki67_pct = ref_ki67["ki67_pct"][i]
            break

    sim_agents = sim["n_tumor_cells"]
    if "tumor_field_cells" in sim:
        sim_tumor = [a + f for a, f in zip(sim_agents, sim["tumor_field_cells"])]
    else:
        sim_tumor = sim_agents

    nonzero = [(d, n) for d, n in zip(sim_days, sim_tumor) if n > 0]
    observed_doubling = float("inf")
    if len(nonzero) >= 2:
        d0, n0 = nonzero[0]
        d1, n1 = nonzero[-1]
        if n1 > n0 and n0 > 0:
            observed_doubling = (d1 - d0) * math.log(2) / math.log(n1 / n0)

    ref_tumor_exp = []
    if nonzero:
        n_init = nonzero[0][1]
        d_start = nonzero[0][0]
        for d in sim_days:
            if d >= d_start and n_init > 0:
                ref_tumor_exp.append(
                    n_init * 2 ** ((d - d_start) / bcc_doubling_days))
            else:
                ref_tumor_exp.append(0)

    has_cycling = "n_tumor_cycling" in sim
    sim_ki67 = []
    mean_ki67 = float("nan")
    if has_cycling:
        sim_ki67 = [100.0 * c / n if n > 0 else 0
                    for c, n in zip(sim["n_tumor_cycling"], sim_tumor)]
        half = len(sim_ki67) // 2
        tail = sim_ki67[half:]
        mean_ki67 = sum(tail) / max(1, len(tail))

    sim_span = sim_days[-1] - (nonzero[0][0] if nonzero else 0)
    n_init = nonzero[0][1] if nonzero else 0
    ref_final = n_init * 2 ** (sim_span / bcc_doubling_days) if n_init > 0 else 0
    obs_final = sim_tumor[-1]
    obs_x = obs_final / n_init if n_init > 0 else 0
    ref_x = ref_final / n_init if n_init > 0 else 0
    sf = surface_fraction(obs_final)

    return dict(
        sim_tumor=sim_tumor, ref_tumor_exp=ref_tumor_exp,
        nonzero=nonzero, observed_doubling=observed_doubling,
        bcc_doubling_days=bcc_doubling_days, bcc_ki67_pct=bcc_ki67_pct,
        has_cycling=has_cycling, sim_ki67=sim_ki67, mean_ki67=mean_ki67,
        sim_span=sim_span, n_init=n_init,
        ref_final=ref_final, obs_final=obs_final,
        obs_x=obs_x, ref_x=ref_x, sf=sf,
    )


def validate_microenvironment(sim, sim_days):
    """Compute microenvironment validation metrics (TGF-b, VEGF, fibronectin, MMP)."""
    ref_tgfb = load_csv(_ref_path("tgfb_kinetics.csv"))
    ref_vegf = load_csv(_ref_path("vegf_kinetics.csv"))
    ref_fn = load_csv(_ref_path("fibronectin_kinetics.csv"))
    ref_mmp = load_csv(_ref_path("mmp_kinetics.csv"))

    sim_tgfb, tgfb_peak = peak_normalize(sim["mean_tgfb_wound"])
    ref_tgfb_at_sim = interpolate(
        ref_tgfb["day"], ref_tgfb["tgfb_normalized"], sim_days)
    tgfb_rmse = compute_rmse(sim_tgfb, ref_tgfb_at_sim)

    sim_vegf, vegf_peak = peak_normalize(sim["mean_vegf_wound"])
    ref_vegf_at_sim = interpolate(
        ref_vegf["day"], ref_vegf["vegf_normalized"], sim_days)
    vegf_rmse = compute_rmse(sim_vegf, ref_vegf_at_sim)

    sim_fn, fn_peak = peak_normalize(sim["mean_fibronectin_wound"])
    ref_fn_at_sim = interpolate(
        ref_fn["day"], ref_fn["fibronectin_normalized"], sim_days)
    fn_rmse = compute_rmse(sim_fn, ref_fn_at_sim)

    sim_mmp, mmp_peak = peak_normalize(sim["mean_mmp_wound"])
    ref_mmp_at_sim = interpolate(
        ref_mmp["day"], ref_mmp["mmp_normalized"], sim_days)
    mmp_rmse = compute_rmse(sim_mmp, ref_mmp_at_sim)

    return dict(
        sim_tgfb=sim_tgfb, ref_tgfb=ref_tgfb,
        ref_tgfb_at_sim=ref_tgfb_at_sim, tgfb_peak=tgfb_peak,
        sim_vegf=sim_vegf, ref_vegf=ref_vegf,
        ref_vegf_at_sim=ref_vegf_at_sim, vegf_peak=vegf_peak,
        sim_fn=sim_fn, ref_fn=ref_fn,
        ref_fn_at_sim=ref_fn_at_sim, fn_peak=fn_peak,
        sim_mmp=sim_mmp, ref_mmp=ref_mmp,
        ref_mmp_at_sim=ref_mmp_at_sim, mmp_peak=mmp_peak,
        tgfb_rmse=tgfb_rmse, vegf_rmse=vegf_rmse,
        fn_rmse=fn_rmse, mmp_rmse=mmp_rmse,
    )


# ---------------------------------------------------------------------------
# Plot functions (draw into provided axes)
# ---------------------------------------------------------------------------

def plot_wound_panels(r, sim_days, axes):
    """Draw 4 wound panels into a 2x2 axes array."""
    cond = r.get("condition", "normal")
    ref_label = "Lit. (diabetic)" if cond == "diabetic" else "Literature"
    ref_kw = dict(REF_KW, label=ref_label)

    ax = axes[0, 0]
    ax.plot(sim_days, r["sim_closure"], **SIM_KW)
    ax.plot(r["ref_closure"]["day"], r["ref_closure"]["closure_pct"], **ref_kw)
    ax.axhline(90, color="gray", linestyle=":", linewidth=0.8, alpha=0.6)
    ax.set_ylabel("Wound closure (%)")
    ax.set_title("Wound Closure" + (" (diabetic)" if cond == "diabetic" else ""))
    ax.set_ylim(-2, 105)
    ax.set_xlim(0, 30)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.text(0.98, 0.05, f"RMSE = {r['closure_rmse']:.1f}%",
            transform=ax.transAxes, ha="right", va="bottom",
            fontsize=8, color="gray")

    ax = axes[0, 1]
    ax.plot(sim_days, r["sim_infl"], **SIM_KW)
    ax.plot(r["ref_infl"]["day"], r["ref_infl"]["inflammation_normalized"], **ref_kw)
    ax.set_ylabel("Inflammation (normalized)")
    ax.set_title("Inflammation Timecourse")
    ax.set_ylim(-0.05, 1.15)
    ax.set_xlim(0, 30)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.text(0.98, 0.95, f"RMSE = {r['inflammation_rmse'] * 100:.1f}%",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=8, color="gray")

    ax = axes[1, 0]
    ax.plot(sim_days, r["sim_neut"], **SIM_KW)
    ax.plot(r["ref_immune"]["day"], r["ref_immune"]["neutrophils_normalized"], **ref_kw)
    ax.set_ylabel("Neutrophils (normalized)")
    ax.set_title("Neutrophil Kinetics")
    ax.set_ylim(-0.05, 1.15)
    ax.set_xlim(0, 30)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.text(0.98, 0.85, f"RMSE = {r['neut_rmse'] * 100:.1f}%",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=8, color="gray")

    ax = axes[1, 1]
    ax.plot(sim_days, r["sim_mac"], **SIM_KW)
    ax.plot(r["ref_immune"]["day"], r["ref_immune"]["macrophages_normalized"], **ref_kw)
    ax.set_ylabel("Macrophages (normalized)")
    ax.set_title("Macrophage Kinetics")
    ax.set_ylim(-0.05, 1.15)
    ax.set_xlim(0, 30)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.text(0.98, 0.85, f"RMSE = {r['mac_rmse'] * 100:.1f}%",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=8, color="gray")

    for a in axes[-1]:
        a.set_xlabel("Time (days)")


def plot_fibroblast_panels(r, sim_days, axes):
    """Draw 3 fibroblast panels into axes array."""
    ax = axes[0]
    ax.plot(sim_days, r["sim_fibro"], **SIM_KW)
    ax.plot(r["ref_fibro"]["day"], r["ref_fibro"]["fibroblasts_normalized"], **REF_KW)
    ax.set_ylabel("Fibroblasts (normalized)")
    ax.set_title("Fibroblast Kinetics")
    ax.set_ylim(-0.05, 1.15)
    ax.set_xlim(0, 30)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.text(0.98, 0.85, f"RMSE = {r['fibro_rmse'] * 100:.1f}%",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=8, color="gray")

    ax = axes[1]
    ax.plot(sim_days, r["sim_myofib"], **SIM_KW)
    ax.plot(r["ref_myofib"]["day"], r["ref_myofib"]["myofibroblasts_normalized"], **REF_KW)
    ax.set_ylabel("Myofibroblasts (normalized)")
    ax.set_title("Myofibroblast Kinetics")
    ax.set_ylim(-0.05, 1.15)
    ax.set_xlim(0, 30)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.text(0.98, 0.85, f"RMSE = {r['myofib_rmse'] * 100:.1f}%",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=8, color="gray")

    ax = axes[2]
    ax.plot(sim_days, r["sim_collagen"], **SIM_KW)
    ax.plot(r["ref_collagen"]["day"], r["ref_collagen"]["collagen_normalized"], **REF_KW)
    ax.set_ylabel("Collagen (normalized)")
    ax.set_title("Collagen Deposition")
    ax.set_ylim(-0.05, 1.15)
    ax.set_xlim(0, 30)
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(True, alpha=0.3)
    ax.text(0.98, 0.85, f"RMSE = {r['collagen_rmse'] * 100:.1f}%",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=8, color="gray")

    axes[-1].set_xlabel("Time (days)")


def plot_tumor_panels(r, sim_days, axes):
    """Draw tumor panels into axes (2 or 3 element array)."""
    td = r["bcc_doubling_days"]

    # Left: log scale growth
    ax = axes[0]
    ax.plot(sim_days, r["sim_tumor"], **SIM_KW)
    if r["ref_tumor_exp"]:
        ax.plot(sim_days, r["ref_tumor_exp"], color=REF_COLOR, linewidth=1.5,
                linestyle="--", label=f"Ref (Td={td:.0f}d)")
    ax.set_ylabel("Tumor cells")
    ax.set_title("Tumor Growth")
    ax.set_yscale("log")
    ax.set_ylim(bottom=1)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    if r["observed_doubling"] < float("inf"):
        ax.text(0.98, 0.05, f"Obs Td = {r['observed_doubling']:.0f}d",
                transform=ax.transAxes, ha="right", va="bottom",
                fontsize=8, color="gray")

    # Right: Ki-67 or linear growth
    ax = axes[1]
    if r["has_cycling"]:
        ax.plot(sim_days, r["sim_ki67"], **SIM_KW)
        ax.axhline(r["bcc_ki67_pct"], color=REF_COLOR, linewidth=1.5,
                   linestyle="--",
                   label=f"BCC Ki-67 = {r['bcc_ki67_pct']:.1f}% (established)")
        ki67_scale = [surface_fraction(n) * 100 for n in r["sim_tumor"]]
        ax.plot(sim_days, ki67_scale, color="#7B9F35", linewidth=1.2,
                linestyle=":", label="Scale-adjusted est.")
        ax.set_ylabel("Cycling fraction (%)")
        ax.set_title("Ki-67 Proxy")
        ax.set_ylim(0, 105)
        ax.text(0.98, 0.85, f"Mean = {r['mean_ki67']:.1f}%",
                transform=ax.transAxes, ha="right", va="top",
                fontsize=8, color="gray")
    else:
        ax.plot(sim_days, r["sim_tumor"], **SIM_KW)
        if r["ref_tumor_exp"]:
            ax.plot(sim_days, r["ref_tumor_exp"], color=REF_COLOR, linewidth=1.5,
                    linestyle="--", label=f"Ref (Td={td:.0f}d)")
        ax.set_ylabel("Tumor cells")
        ax.set_title("Tumor Growth (linear)")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    axes[-1].set_xlabel("Time (days)")


def plot_microenvironment_panels(r, sim_days, axes):
    """Draw 4 microenvironment panels into a 2x2 axes array."""
    ax = axes[0, 0]
    ax.plot(sim_days, r["sim_tgfb"], **SIM_KW)
    ax.plot(r["ref_tgfb"]["day"], r["ref_tgfb"]["tgfb_normalized"], **REF_KW)
    ax.set_ylabel("TGF-b (normalized)")
    ax.set_title("TGF-b1 Kinetics")
    ax.set_ylim(-0.05, 1.15)
    ax.set_xlim(0, 30)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.text(0.98, 0.85, f"RMSE = {r['tgfb_rmse'] * 100:.1f}%",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=8, color="gray")

    ax = axes[0, 1]
    ax.plot(sim_days, r["sim_vegf"], **SIM_KW)
    ax.plot(r["ref_vegf"]["day"], r["ref_vegf"]["vegf_normalized"], **REF_KW)
    ax.set_ylabel("VEGF (normalized)")
    ax.set_title("VEGF Kinetics")
    ax.set_ylim(-0.05, 1.15)
    ax.set_xlim(0, 30)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.text(0.98, 0.85, f"RMSE = {r['vegf_rmse'] * 100:.1f}%",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=8, color="gray")

    ax = axes[1, 0]
    ax.plot(sim_days, r["sim_fn"], **SIM_KW)
    ax.plot(r["ref_fn"]["day"], r["ref_fn"]["fibronectin_normalized"], **REF_KW)
    ax.set_ylabel("Fibronectin (normalized)")
    ax.set_title("Fibronectin Kinetics")
    ax.set_ylim(-0.05, 1.15)
    ax.set_xlim(0, 30)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.text(0.98, 0.85, f"RMSE = {r['fn_rmse'] * 100:.1f}%",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=8, color="gray")

    ax = axes[1, 1]
    ax.plot(sim_days, r["sim_mmp"], **SIM_KW)
    ax.plot(r["ref_mmp"]["day"], r["ref_mmp"]["mmp_normalized"], **REF_KW)
    ax.set_ylabel("MMP (normalized)")
    ax.set_title("MMP Activity")
    ax.set_ylim(-0.05, 1.15)
    ax.set_xlim(0, 30)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.text(0.98, 0.85, f"RMSE = {r['mmp_rmse'] * 100:.1f}%",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=8, color="gray")

    for a in axes[-1]:
        a.set_xlabel("Time (days)")


# ---------------------------------------------------------------------------
# Print summary
# ---------------------------------------------------------------------------

def print_summary(wound=None, fibroblast=None, tumor=None, microenv=None):
    """Print one unified validation summary block.

    Only shows modules that produced results. Disabled modules are
    silently omitted.
    """
    cond = wound.get("condition", "normal") if wound else "normal"
    cond_label = f" [{cond}]" if cond != "normal" else ""
    print("=" * 60)
    print(f"  skibidy Validation Summary{cond_label}")
    print("=" * 60)
    if wound:
        w = wound
        print(f"  Wound closure       RMSE = {w['closure_rmse']:.2f} %   max = {w['closure_max']:.2f} %")
        print(f"    Inflammatory 0-3d RMSE = {w['infl_rmse']:.2f} %")
        print(f"    Proliferative 3-14d     = {w['prolif_rmse']:.2f} %")
        print(f"    Remodeling 14-28d       = {w['remod_rmse']:.2f} %")
        print(f"  Inflammation        RMSE = {w['inflammation_rmse'] * 100:.2f} %")
        print(f"  Neutrophils         RMSE = {w['neut_rmse'] * 100:.2f} %")
        print(f"  Macrophages         RMSE = {w['mac_rmse'] * 100:.2f} %")
    if fibroblast:
        f = fibroblast
        print(f"  Fibroblasts         RMSE = {f['fibro_rmse'] * 100:.2f} %")
        print(f"  Myofibroblasts      RMSE = {f['myofib_rmse'] * 100:.2f} %")
        print(f"  Collagen            RMSE = {f['collagen_rmse'] * 100:.2f} %")
    if microenv:
        m = microenv
        print(f"  TGF-b               RMSE = {m['tgfb_rmse'] * 100:.2f} %")
        print(f"  VEGF                RMSE = {m['vegf_rmse'] * 100:.2f} %")
        print(f"  Fibronectin         RMSE = {m['fn_rmse'] * 100:.2f} %")
        print(f"  MMP                 RMSE = {m['mmp_rmse'] * 100:.2f} %")
    if tumor:
        t = tumor
        print(f"  Tumor span          {t['sim_span']:.0f}d  Td: Obs={t['observed_doubling']:.0f}d Ref={t['bcc_doubling_days']:.0f}d")
        print(f"  Tumor growth        Obs={t['obs_final']:.0f} ({t['obs_x']:.1f}x)  Ref={t['ref_final']:.0f} ({t['ref_x']:.1f}x)")
        if t["has_cycling"]:
            print(f"  Tumor Ki-67 proxy   Sim={t['mean_ki67']:.1f}%  Ref={t['bcc_ki67_pct']:.1f}%  Scale-adj ~{t['sf']*100:.0f}%")
    print("=" * 60)
