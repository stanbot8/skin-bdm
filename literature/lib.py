"""Shared utilities for validation scripts.

Provides CSV loading, interpolation, normalization, error metrics,
module-level validate/plot/print functions used by validate_all.py
and the individual compare_*.py wrappers.
"""

import csv
import math
import os
from batch.lib import load_csv


# ---------------------------------------------------------------------------
# Low-level helpers
# ---------------------------------------------------------------------------

def plots_dir(csv_path):
    """Derive plots output directory from a metrics CSV path.

    Returns <csv_parent>/plots/, e.g. studies/wound/results/plots/.
    """
    return os.path.join(os.path.dirname(os.path.abspath(csv_path)), "plots")


# load_csv imported from batch.lib (canonical version)
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


def compute_rmse_ci(sim, ref):
    """Compute RMSE with 95% confidence interval via delta method.

    The CI quantifies uncertainty in the RMSE estimate arising from
    finite timepoint sampling. Useful for distinguishing "excellent fit"
    from "mediocre but within noise".

    Returns dict with rmse, se, ci_lo, ci_hi.
    """
    if not sim:
        return dict(rmse=0.0, se=0.0, ci_lo=0.0, ci_hi=0.0)
    sq_errors = [(s - r) ** 2 for s, r in zip(sim, ref)]
    n = len(sq_errors)
    mse = sum(sq_errors) / n
    rmse = math.sqrt(mse)
    # Delta method: SE(RMSE) = std(squared_errors) / (2 * RMSE * sqrt(n))
    var_sq = sum((e - mse) ** 2 for e in sq_errors) / max(n - 1, 1)
    se = math.sqrt(var_sq) / (2 * max(rmse, 1e-10) * math.sqrt(n))
    ci_lo = max(0, rmse - 1.96 * se)
    ci_hi = rmse + 1.96 * se
    return dict(rmse=rmse, se=se, ci_lo=ci_lo, ci_hi=ci_hi)


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
_STUDIES_DIR = os.path.join(os.path.dirname(__file__), os.pardir, "studies")

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
    "diabetic_mmp_kinetics.csv": "diabetic",
    "diabetic_tgfb_kinetics.csv": "diabetic",
    "tumor_doubling_time.csv": "tumor",
    "tumor_growth_rate.csv": "tumor",
    "tumor_proliferation_index.csv": "tumor",
    "ph_kinetics.csv": "ph",
    "burn_closure_kinetics.csv": "burn",
    "burn_inflammation_timecourse.csv": "burn",
    "burn_immune_cell_kinetics.csv": "burn",
    "pressure_closure_kinetics.csv": "pressure",
    "pressure_inflammation_timecourse.csv": "pressure",
    "pressure_immune_cell_kinetics.csv": "pressure",
    "surgical_closure_kinetics.csv": "wound",
    "surgical_inflammation_timecourse.csv": "wound",
    "surgical_immune_cell_kinetics.csv": "wound",
}


def _ref_path(filename):
    """Resolve a consensus CSV filename to its module data path."""
    module = _FILE_MODULE_MAP[filename]
    return os.path.join(_MODULES_DIR, module, "data", filename)


def _study_ref_path(study, module, filename):
    """Resolve a consensus CSV filename for a study-scoped module."""
    return os.path.join(_STUDIES_DIR, study, "modules", module, "data", filename)


# ---------------------------------------------------------------------------
# Condition detection (normal vs diabetic)
# ---------------------------------------------------------------------------

def detect_condition(config_path=None):
    """Detect simulation condition from bdm.toml.

    Checks for condition-specific sections ([skin.diabetic], [skin.burn],
    [skin.pressure]) with mode = true, or study = "surgical".
    Returns "diabetic", "burn", "pressure", "surgical", or "normal".
    """
    if config_path is None:
        config_path = os.path.join(os.path.dirname(__file__), os.pardir,
                                   "bdm.toml")
    try:
        with open(config_path) as f:
            current_section = ""
            for line in f:
                stripped = line.strip()
                if stripped.startswith("["):
                    current_section = stripped
                elif ("mode" in stripped or "enabled" in stripped) and "=" in stripped:
                    key = stripped.split("=", 1)[0].strip()
                    if key not in ("mode", "enabled"):
                        continue
                    val = stripped.split("=", 1)[1].strip().lower()
                    if val.startswith("true"):
                        if current_section == "[skin.diabetic]":
                            return "diabetic"
                        if current_section == "[skin.burn]":
                            return "burn"
                        if current_section == "[skin.pressure]":
                            return "pressure"
                        if current_section == "[skin.ra]":
                            return "rheumatoid"
                elif "output_dir" in stripped and "=" in stripped:
                    val = stripped.split("=", 1)[1].strip().strip('"').strip("'").lower()
                    if "surgical" in val:
                        return "surgical"
                elif ("study" in stripped and "=" in stripped
                      and current_section in ("[skin]", "[simulation]")):
                    val = stripped.split("=", 1)[1].strip().strip('"').strip("'").lower()
                    if val == "surgical":
                        return "surgical"
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
    """Return (has_wound, has_fibroblast, has_tumor, has_microenv, has_ph, has_ra) booleans."""
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
    has_ph = (has_wound
              and "mean_ph_wound" in sim
              and max(sim["mean_ph_wound"]) > 0)
    has_ra = ("mean_tnf_alpha_wound" in sim
              and max(sim["mean_tnf_alpha_wound"]) > 0)
    return has_wound, has_fibroblast, has_tumor, has_microenv, has_ph, has_ra


# ---------------------------------------------------------------------------
# Validate functions (pure computation, no I/O)
# ---------------------------------------------------------------------------

def validate_wound(sim, sim_days, condition="normal"):
    """Compute wound validation metrics. Returns result dict.

    Supports conditions: normal, diabetic, burn, pressure, surgical.
    Each loads condition-specific reference curves from the corresponding
    module data directory.
    """
    if condition == "diabetic":
        ref_closure = load_csv(_ref_path("diabetic_closure_kinetics.csv"))
        ref_infl = load_csv(_ref_path("diabetic_inflammation_timecourse.csv"))
        ref_immune = load_csv(_ref_path("diabetic_immune_cell_kinetics.csv"))
    elif condition == "burn":
        ref_closure = load_csv(_ref_path("burn_closure_kinetics.csv"))
        ref_infl = load_csv(_ref_path("burn_inflammation_timecourse.csv"))
        ref_immune = load_csv(_ref_path("burn_immune_cell_kinetics.csv"))
    elif condition == "pressure":
        ref_closure = load_csv(_ref_path("pressure_closure_kinetics.csv"))
        ref_infl = load_csv(_ref_path("pressure_inflammation_timecourse.csv"))
        ref_immune = load_csv(_ref_path("pressure_immune_cell_kinetics.csv"))
    elif condition == "surgical":
        ref_closure = load_csv(_ref_path("surgical_closure_kinetics.csv"))
        ref_infl = load_csv(_ref_path("surgical_inflammation_timecourse.csv"))
        ref_immune = load_csv(_ref_path("surgical_immune_cell_kinetics.csv"))
    else:
        ref_closure = load_csv(_ref_path("closure_kinetics_punch_biopsy.csv"))
        ref_infl = load_csv(_ref_path("inflammation_timecourse.csv"))
        ref_immune = load_csv(_ref_path("immune_cell_kinetics.csv"))

    sim_closure = sim["wound_closure_pct"]
    ref_closure_at_sim = interpolate(
        ref_closure["day"], ref_closure["closure_pct"], sim_days)
    closure_rmse = compute_rmse(sim_closure, ref_closure_at_sim)
    closure_ci = compute_rmse_ci(sim_closure, ref_closure_at_sim)
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
        closure_ci=closure_ci,
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


def validate_microenvironment(sim, sim_days, condition="normal"):
    """Compute microenvironment validation metrics (TGF-b, VEGF, fibronectin, MMP).

    When condition="diabetic", uses diabetic-specific reference curves for
    MMP (sustained elevation; Lobmann et al. 2002) and TGF-b (delayed peak;
    Mirza & Koh 2011).
    """
    if condition == "diabetic":
        ref_tgfb = load_csv(_ref_path("diabetic_tgfb_kinetics.csv"))
        ref_mmp = load_csv(_ref_path("diabetic_mmp_kinetics.csv"))
    else:
        ref_tgfb = load_csv(_ref_path("tgfb_kinetics.csv"))
        ref_mmp = load_csv(_ref_path("mmp_kinetics.csv"))
    ref_vegf = load_csv(_ref_path("vegf_kinetics.csv"))
    ref_fn = load_csv(_ref_path("fibronectin_kinetics.csv"))

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


def validate_ph(sim, sim_days):
    """Compute wound pH (alkalinity) validation metrics. Returns result dict.

    pH is modeled as normalized alkalinity: 1.0 = fresh wound (pH 7.4),
    0.0 = healed acidic skin (pH 5.5). Reference: Schneider et al. 2007.
    """
    ref_ph = load_csv(_ref_path("ph_kinetics.csv"))

    sim_ph = sim["mean_ph_wound"]
    ref_ph_at_sim = interpolate(
        ref_ph["day"], ref_ph["ph_alkalinity_normalized"], sim_days)
    ph_rmse = compute_rmse(sim_ph, ref_ph_at_sim)

    return dict(
        sim_ph=sim_ph, ref_ph=ref_ph,
        ref_ph_at_sim=ref_ph_at_sim,
        ph_rmse=ph_rmse,
    )


def validate_ra(sim, sim_days):
    """Compute RA validation metrics. Returns result dict.

    Six observables: TNF-alpha, IL-6 (peak-normalized), cartilage and bone
    (absolute integrity), T cell density and synovial pannus (peak-normalized).
    Reference kinetics derived from Feldmann & Maini 2003, Kishimoto 2005,
    McInnes & Schett 2011, Schett & Gravallese 2012, Firestein 2003.
    """
    ref_tnf = load_csv(_study_ref_path("rheumatoid", "rheumatoid",
                                        "tnf_alpha_kinetics.csv"))
    ref_il6 = load_csv(_study_ref_path("rheumatoid", "rheumatoid",
                                        "il6_kinetics.csv"))
    ref_cart = load_csv(_study_ref_path("rheumatoid", "rheumatoid",
                                         "cartilage_erosion.csv"))
    ref_bone = load_csv(_study_ref_path("rheumatoid", "rheumatoid",
                                         "bone_erosion.csv"))
    ref_tcell = load_csv(_study_ref_path("rheumatoid", "rheumatoid",
                                          "tcell_density.csv"))
    ref_syn = load_csv(_study_ref_path("rheumatoid", "rheumatoid",
                                        "synovial_pannus.csv"))

    # TNF-alpha: peak-normalize simulation, compare to reference
    sim_tnf, tnf_peak = peak_normalize(sim["mean_tnf_alpha_wound"])
    ref_tnf_at_sim = interpolate(
        ref_tnf["day"], ref_tnf["tnf_alpha_normalized"], sim_days)
    tnf_rmse = compute_rmse(sim_tnf, ref_tnf_at_sim)

    # IL-6: peak-normalize simulation
    sim_il6, il6_peak = peak_normalize(sim["mean_il6_wound"])
    ref_il6_at_sim = interpolate(
        ref_il6["day"], ref_il6["il6_normalized"], sim_days)
    il6_rmse = compute_rmse(sim_il6, ref_il6_at_sim)

    # Cartilage: normalize to initial value (wound cylinder includes non-cartilage voxels)
    sim_cart_raw = sim["mean_cartilage_wound"]
    cart_init = sim_cart_raw[0] if sim_cart_raw[0] > 1e-6 else 1.0
    sim_cart = [v / cart_init for v in sim_cart_raw]
    ref_cart_at_sim = interpolate(
        ref_cart["day"], ref_cart["cartilage_integrity"], sim_days)
    cart_rmse = compute_rmse(sim_cart, ref_cart_at_sim)

    # Phase-specific RMSE for TNF
    tnf_flare_rmse = phase_rmse(sim_days, sim_tnf, ref_tnf_at_sim, 0, 7)
    tnf_chronic_rmse = phase_rmse(sim_days, sim_tnf, ref_tnf_at_sim, 7, 30)

    # Bone: absolute integrity comparison (slower erosion than cartilage)
    has_bone = ("mean_bone_wound" in sim
                and max(sim["mean_bone_wound"]) > 0)
    bone_rmse = 0.0
    sim_bone = []
    ref_bone_at_sim = []
    if has_bone:
        sim_bone_raw = sim["mean_bone_wound"]
        bone_init = sim_bone_raw[0] if sim_bone_raw[0] > 1e-6 else 1.0
        sim_bone = [v / bone_init for v in sim_bone_raw]
        ref_bone_at_sim = interpolate(
            ref_bone["day"], ref_bone["bone_integrity"], sim_days)
        bone_rmse = compute_rmse(sim_bone, ref_bone_at_sim)

    # T cell density: peak-normalize
    has_tcell = ("mean_tcell_wound" in sim
                 and max(sim["mean_tcell_wound"]) > 0)
    tcell_rmse = 0.0
    sim_tcell = []
    tcell_peak = 0.0
    ref_tcell_at_sim = []
    if has_tcell:
        sim_tcell, tcell_peak = peak_normalize(sim["mean_tcell_wound"])
        ref_tcell_at_sim = interpolate(
            ref_tcell["day"], ref_tcell["tcell_normalized"], sim_days)
        tcell_rmse = compute_rmse(sim_tcell, ref_tcell_at_sim)

    # Synovial pannus: peak-normalize
    has_syn = ("mean_synovial_wound" in sim
               and max(sim["mean_synovial_wound"]) > 0)
    syn_rmse = 0.0
    sim_syn = []
    syn_peak = 0.0
    ref_syn_at_sim = []
    if has_syn:
        sim_syn, syn_peak = peak_normalize(sim["mean_synovial_wound"])
        ref_syn_at_sim = interpolate(
            ref_syn["day"], ref_syn["synovial_normalized"], sim_days)
        syn_rmse = compute_rmse(sim_syn, ref_syn_at_sim)

    return dict(
        sim_tnf=sim_tnf, ref_tnf=ref_tnf,
        ref_tnf_at_sim=ref_tnf_at_sim, tnf_peak=tnf_peak,
        sim_il6=sim_il6, ref_il6=ref_il6,
        ref_il6_at_sim=ref_il6_at_sim, il6_peak=il6_peak,
        sim_cart=sim_cart, ref_cart=ref_cart,
        ref_cart_at_sim=ref_cart_at_sim,
        tnf_rmse=tnf_rmse, il6_rmse=il6_rmse, cart_rmse=cart_rmse,
        tnf_flare_rmse=tnf_flare_rmse, tnf_chronic_rmse=tnf_chronic_rmse,
        has_bone=has_bone, sim_bone=sim_bone, ref_bone=ref_bone,
        ref_bone_at_sim=ref_bone_at_sim, bone_rmse=bone_rmse,
        has_tcell=has_tcell, sim_tcell=sim_tcell, ref_tcell=ref_tcell,
        ref_tcell_at_sim=ref_tcell_at_sim, tcell_peak=tcell_peak,
        tcell_rmse=tcell_rmse,
        has_syn=has_syn, sim_syn=sim_syn, ref_syn=ref_syn,
        ref_syn_at_sim=ref_syn_at_sim, syn_peak=syn_peak,
        syn_rmse=syn_rmse,
    )


# ---------------------------------------------------------------------------
# Plot functions (draw into provided axes)
# ---------------------------------------------------------------------------

def plot_wound_panels(r, sim_days, axes):
    """Draw 4 wound panels into a 2x2 axes array."""
    cond = r.get("condition", "normal")
    ref_label = f"Lit. ({cond})" if cond != "normal" else "Literature"
    ref_kw = dict(REF_KW, label=ref_label)

    ax = axes[0, 0]
    ax.plot(sim_days, r["sim_closure"], **SIM_KW)
    ax.plot(r["ref_closure"]["day"], r["ref_closure"]["closure_pct"], **ref_kw)
    ax.axhline(90, color="gray", linestyle=":", linewidth=0.8, alpha=0.6)
    ax.set_ylabel("Wound closure (%)")
    ax.set_title("Wound Closure" + (f" ({cond})" if cond != "normal" else ""))
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


def plot_ph_panel(r, sim_days, ax):
    """Draw wound pH (alkalinity) panel into a single axes."""
    ax.plot(sim_days, r["sim_ph"], **SIM_KW)
    ax.plot(r["ref_ph"]["day"], r["ref_ph"]["ph_alkalinity_normalized"], **REF_KW)
    ax.set_ylabel("Wound alkalinity (normalized)")
    ax.set_title("Wound pH Recovery")
    ax.set_ylim(-0.05, 1.15)
    ax.set_xlim(0, 30)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.text(0.98, 0.85, f"RMSE = {r['ph_rmse'] * 100:.1f}%",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=8, color="gray")
    ax.set_xlabel("Time (days)")


def plot_ra_panels(r, sim_days, axes):
    """Draw RA panels into axes array.

    Expects 3 axes for core panels (TNF, IL-6, cartilage).
    If extended data is available (bone, T cell, synovial), expects 6 axes.
    """
    ax = axes[0]
    ax.plot(sim_days, r["sim_tnf"], **SIM_KW)
    ax.plot(r["ref_tnf"]["day"], r["ref_tnf"]["tnf_alpha_normalized"], **REF_KW)
    ax.set_ylabel("TNF-alpha (normalized)")
    ax.set_title("TNF-alpha Kinetics")
    ax.set_ylim(-0.05, 1.15)
    ax.set_xlim(0, 30)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.text(0.98, 0.85, f"RMSE = {r['tnf_rmse'] * 100:.1f}%",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=8, color="gray")

    ax = axes[1]
    ax.plot(sim_days, r["sim_il6"], **SIM_KW)
    ax.plot(r["ref_il6"]["day"], r["ref_il6"]["il6_normalized"], **REF_KW)
    ax.set_ylabel("IL-6 (normalized)")
    ax.set_title("IL-6 Kinetics")
    ax.set_ylim(-0.05, 1.15)
    ax.set_xlim(0, 30)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.text(0.98, 0.85, f"RMSE = {r['il6_rmse'] * 100:.1f}%",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=8, color="gray")

    ax = axes[2]
    ax.plot(sim_days, r["sim_cart"], **SIM_KW)
    ax.plot(r["ref_cart"]["day"], r["ref_cart"]["cartilage_integrity"], **REF_KW)
    ax.set_ylabel("Cartilage integrity")
    ax.set_title("Cartilage Erosion")
    ax.set_ylim(-0.05, 1.15)
    ax.set_xlim(0, 30)
    ax.legend(fontsize=8, loc="lower left")
    ax.grid(True, alpha=0.3)
    ax.text(0.98, 0.15, f"RMSE = {r['cart_rmse'] * 100:.1f}%",
            transform=ax.transAxes, ha="right", va="bottom",
            fontsize=8, color="gray")

    # Extended panels (bone, T cell, synovial) if axes available
    if len(axes) < 6:
        axes[-1].set_xlabel("Time (days)")
        return

    if r.get("has_bone") and r["sim_bone"]:
        ax = axes[3]
        ax.plot(sim_days, r["sim_bone"], **SIM_KW)
        ax.plot(r["ref_bone"]["day"], r["ref_bone"]["bone_integrity"], **REF_KW)
        ax.set_ylabel("Bone integrity")
        ax.set_title("Subchondral Bone Erosion")
        ax.set_ylim(-0.05, 1.15)
        ax.set_xlim(0, 30)
        ax.legend(fontsize=8, loc="lower left")
        ax.grid(True, alpha=0.3)
        ax.text(0.98, 0.15, f"RMSE = {r['bone_rmse'] * 100:.1f}%",
                transform=ax.transAxes, ha="right", va="bottom",
                fontsize=8, color="gray")
    else:
        axes[3].set_visible(False)

    if r.get("has_tcell") and r["sim_tcell"]:
        ax = axes[4]
        ax.plot(sim_days, r["sim_tcell"], **SIM_KW)
        ax.plot(r["ref_tcell"]["day"], r["ref_tcell"]["tcell_normalized"], **REF_KW)
        ax.set_ylabel("T cell density (normalized)")
        ax.set_title("T Cell Infiltration")
        ax.set_ylim(-0.05, 1.15)
        ax.set_xlim(0, 30)
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
        ax.text(0.98, 0.85, f"RMSE = {r['tcell_rmse'] * 100:.1f}%",
                transform=ax.transAxes, ha="right", va="top",
                fontsize=8, color="gray")
    else:
        axes[4].set_visible(False)

    if r.get("has_syn") and r["sim_syn"]:
        ax = axes[5]
        ax.plot(sim_days, r["sim_syn"], **SIM_KW)
        ax.plot(r["ref_syn"]["day"], r["ref_syn"]["synovial_normalized"], **REF_KW)
        ax.set_ylabel("Pannus density (normalized)")
        ax.set_title("Synovial Pannus Growth")
        ax.set_ylim(-0.05, 1.15)
        ax.set_xlim(0, 30)
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
        ax.text(0.98, 0.85, f"RMSE = {r['syn_rmse'] * 100:.1f}%",
                transform=ax.transAxes, ha="right", va="top",
                fontsize=8, color="gray")
    else:
        axes[5].set_visible(False)

    axes[-1].set_xlabel("Time (days)")


# ---------------------------------------------------------------------------
# Print summary
# ---------------------------------------------------------------------------

def print_summary(wound=None, fibroblast=None, tumor=None, microenv=None,
                  ph=None, ra=None):
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
        ci = w.get("closure_ci")
        if ci:
            print(f"  Wound closure       RMSE = {w['closure_rmse']:.2f} %   "
                  f"95% CI [{ci['ci_lo']:.2f}, {ci['ci_hi']:.2f}]   max = {w['closure_max']:.2f} %")
        else:
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
    if ph:
        print(f"  Wound pH            RMSE = {ph['ph_rmse'] * 100:.2f} %")
    if tumor:
        t = tumor
        print(f"  Tumor span          {t['sim_span']:.0f}d  Td: Obs={t['observed_doubling']:.0f}d Ref={t['bcc_doubling_days']:.0f}d")
        print(f"  Tumor growth        Obs={t['obs_final']:.0f} ({t['obs_x']:.1f}x)  Ref={t['ref_final']:.0f} ({t['ref_x']:.1f}x)")
        if t["has_cycling"]:
            print(f"  Tumor Ki-67 proxy   Sim={t['mean_ki67']:.1f}%  Ref={t['bcc_ki67_pct']:.1f}%  Scale-adj ~{t['sf']*100:.0f}%")
    if ra:
        a = ra
        print(f"  TNF-alpha           RMSE = {a['tnf_rmse'] * 100:.2f} %")
        print(f"    Flare 0-7d        RMSE = {a['tnf_flare_rmse'] * 100:.2f} %")
        print(f"    Chronic 7-30d     RMSE = {a['tnf_chronic_rmse'] * 100:.2f} %")
        print(f"  IL-6                RMSE = {a['il6_rmse'] * 100:.2f} %")
        print(f"  Cartilage integrity RMSE = {a['cart_rmse'] * 100:.2f} %")
        if a.get("has_bone"):
            print(f"  Bone integrity      RMSE = {a['bone_rmse'] * 100:.2f} %")
        if a.get("has_tcell"):
            print(f"  T cell density      RMSE = {a['tcell_rmse'] * 100:.2f} %")
        if a.get("has_syn"):
            print(f"  Synovial pannus     RMSE = {a['syn_rmse'] * 100:.2f} %")
    print("=" * 60)
