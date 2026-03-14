#!/usr/bin/env python3
"""Clinical Calibration Tool for Skin BioDynaMo (SkiBiDy).

Maps real patient clinical measurements to simulation parameters, generating
a calibrated profile TOML file. Each mapping is grounded in evidence-based
relationships between clinical observables and wound healing biology.

Usage:
    python3 scripts/clinical_calibrate.py --age 72 --bmi 32 --hba1c 7.2 \\
        --wound-type pressure --wound-size 30 --location sacrum --recommend

Author: stanbot8
"""

import argparse
import sys
from datetime import datetime
from pathlib import Path

# ---------------------------------------------------------------------------
# Default simulation baseline values (healthy adult)
# ---------------------------------------------------------------------------
DEFAULTS = {
    "prolif_factor": 1.0,
    "migration_factor": 1.0,
    "collagen_factor": 1.0,
    "perfusion_basal": 1.0,
    "angio_rate": 0.005,
    "stem_fraction": 0.12,
    "aging_factor": 1.0,
    "m1_duration_factor": 1.0,
    "resolution_factor": 1.0,
    "neutrophil_factor": 1.0,
    "efferocytosis_factor": 1.0,
    "inflammation_baseline": 0.0,
    "inflammation_sensitivity": 1.0,
    "ros_baseline": 0.0,
    "vegf_factor": 1.0,
    "mmp_factor": 1.0,
    "fibroblast_activation_factor": 1.0,
    "fibroblast_lifespan_factor": 1.0,
    "biofilm_seed_delay_h": None,       # None means no biofilm
    "biofilm_growth_rate": None,
    "g1_duration": 7.0,
    "s_duration": 6.0,
    "water_recovery_rate": 0.05,
    # Blood module
    "hemoglobin": 1.0,
    "anticoagulant_level": 0.0,
    "platelet_count": 1.0,
    "bleed_rate": 0.02,
}


def parse_args():
    """Parse command line arguments for patient clinical data."""
    p = argparse.ArgumentParser(
        description="Generate a calibrated simulation profile from patient data."
    )
    p.add_argument("--age", type=float, default=50,
                   help="Patient age in years (default: 50)")
    p.add_argument("--bmi", type=float, default=25,
                   help="Body mass index in kg/m2 (default: 25)")
    p.add_argument("--hba1c", type=float, default=5.5,
                   help="HbA1c percentage (default: 5.5, diabetic >6.5)")
    p.add_argument("--wound-type", choices=["punch", "incision", "pressure", "burn"],
                   default="punch", help="Wound type (default: punch)")
    p.add_argument("--wound-size", type=float, default=12,
                   help="Wound diameter in mm (default: 12)")
    p.add_argument("--wound-depth", choices=["superficial", "partial", "full"],
                   default="partial", help="Wound depth (default: partial)")
    p.add_argument("--smoking", action="store_true",
                   help="Patient is a current smoker")
    p.add_argument("--immunosuppressed", action="store_true",
                   help="Patient is on immunosuppressive therapy")
    p.add_argument("--infection", choices=["none", "low", "moderate", "severe"],
                   default="none", help="Wound infection level (default: none)")
    p.add_argument("--location", choices=["sacrum", "heel", "forearm", "abdomen", "face"],
                   default="forearm", help="Wound anatomical location (default: forearm)")
    p.add_argument("--albumin", type=float, default=4.0,
                   help="Serum albumin in g/dL (default: 4.0, malnutrition <3.5)")
    p.add_argument("--wbc", type=float, default=7.0,
                   help="White blood cell count in 10^3/uL (default: 7.0)")
    p.add_argument("--hemoglobin", type=float, default=14.0,
                   help="Hemoglobin in g/dL (default: 14.0, anemia <10)")
    p.add_argument("--platelets", type=float, default=250,
                   help="Platelet count in 10^3/uL (default: 250, low <100)")
    p.add_argument("--anticoagulant", choices=["none", "prophylactic", "therapeutic"],
                   default="none", help="Anticoagulant therapy level")
    p.add_argument("--output", type=str, default="profiles/patient.toml",
                   help="Output TOML path (default: profiles/patient.toml)")
    p.add_argument("--recommend", action="store_true",
                   help="Print treatment recommendations after generating profile")
    return p.parse_args()


# ---------------------------------------------------------------------------
# Calibration mappings: each function modifies params dict in place and
# returns a list of annotation strings documenting what changed.
# ---------------------------------------------------------------------------

def apply_age(params, age):
    """Age-related decline in healing capacity (Sgonc & Gruber 2013)."""
    notes = []
    if age > 30:
        params["aging_factor"] = max(0.3, 1.0 - (age - 30) * 0.01)
        params["stem_fraction"] = max(0.04, 0.12 - (age - 30) * 0.001)
        notes.append(f"Age {age}: aging_factor={params['aging_factor']:.2f}, "
                     f"stem_fraction={params['stem_fraction']:.3f}")
    if age > 60:
        params["g1_duration"] = 9.0
        params["s_duration"] = 8.0
        params["water_recovery_rate"] = 0.03
        notes.append("Age >60: slower cell cycle, reduced hydration")
    return notes


def apply_bmi(params, bmi):
    """Obesity impairs proliferation and perfusion (Pierpont et al. 2014)."""
    notes = []
    if bmi > 30:
        params["prolif_factor"] *= 0.8
        params["perfusion_basal"] *= 0.85
        notes.append(f"BMI {bmi}: prolif x0.8, perfusion x0.85 (obese)")
    if bmi > 35:
        params["prolif_factor"] *= 0.875    # cumulative with above: 0.7
        params["perfusion_basal"] *= 0.824  # cumulative: 0.7
        params["inflammation_baseline"] += 0.0003
        notes.append(f"BMI {bmi}: severe obesity adjustments, chronic inflammation +0.0003")
    return notes


def apply_hba1c(params, hba1c):
    """Diabetic cascade from glycemic impairment (Brownlee 2005, Mirza & Koh 2011)."""
    notes = []
    if hba1c <= 5.5:
        return notes

    excess = hba1c - 5.5
    severity = excess / 4.0  # normalized 0 to ~1

    params["m1_duration_factor"] = 1.0 + excess * 0.5
    params["prolif_factor"] *= max(0.3, 1.0 - excess * 0.15)
    params["collagen_factor"] *= max(0.3, 1.0 - excess * 0.12)
    params["inflammation_sensitivity"] = 1.0 + severity * 1.0
    notes.append(f"HbA1c {hba1c}%: M1 duration x{params['m1_duration_factor']:.1f}, "
                 f"immune dysfunction severity {severity:.2f}")

    if hba1c > 6.5:
        # Full diabetic profile activation
        params["diabetic_mode"] = True
        params["resolution_factor"] = max(0.3, 1.0 - severity * 0.7)
        params["efferocytosis_factor"] = max(0.3, 1.0 - severity * 0.5)
        params["neutrophil_factor"] = 1.0 + severity * 0.8
        params["vegf_factor"] = max(0.3, 1.0 - severity * 0.6)
        params["mmp_factor"] = 1.0 + severity * 1.5
        params["fibroblast_activation_factor"] = 1.0 + severity * 1.0
        params["fibroblast_lifespan_factor"] = max(0.5, 1.0 - severity * 0.4)
        params["inflammation_baseline"] += severity * 0.001
        notes.append(f"Diabetic mode enabled: resolution={params['resolution_factor']:.2f}, "
                     f"VEGF={params['vegf_factor']:.2f}")
    return notes


def apply_smoking(params):
    """Smoking impairs perfusion and angiogenesis (Silverstein 1992, Sorensen 2012)."""
    params["perfusion_basal"] *= 0.7
    params["angio_rate"] *= 0.5
    params["ros_baseline"] += 0.001
    return ["Smoking: perfusion x0.7, angiogenesis x0.5, ROS +0.001"]


def apply_immunosuppression(params):
    """Immunosuppressive therapy reduces all immune functions (Guo & DiPietro 2010)."""
    for key in ["neutrophil_factor", "efferocytosis_factor", "resolution_factor"]:
        params[key] *= 0.5
    params["m1_duration_factor"] *= 0.5
    return ["Immunosuppressed: all immune factors x0.5"]


def apply_infection(params, level):
    """Biofilm seeding based on clinical infection severity (Bjarnsholt 2013)."""
    if level == "none":
        return []

    delays = {"low": 72, "moderate": 24, "severe": 0}
    rates = {"low": 0.001, "moderate": 0.005, "severe": 0.01}

    params["biofilm_seed_delay_h"] = delays[level]
    params["biofilm_growth_rate"] = rates[level]
    return [f"Infection ({level}): biofilm delay={delays[level]}h, "
            f"growth={rates[level]}"]


def apply_location(params, location):
    """Anatomical location affects perfusion (Guo & DiPietro 2010)."""
    multipliers = {
        "sacrum": 0.8,
        "heel": 0.7,
        "forearm": 1.0,
        "abdomen": 0.9,
        "face": 1.3,
    }
    params["perfusion_basal"] *= multipliers[location]
    extra = ""
    if location == "sacrum":
        extra = ", pressure/shear risk"
    elif location == "heel":
        extra = ", weight-bearing stress"
    elif location == "face":
        extra = ", excellent vascular supply"
    return [f"Location ({location}): perfusion x{multipliers[location]}{extra}"]


def apply_albumin(params, albumin):
    """Nutritional status affects collagen and proliferation (Stechmiller 2010)."""
    notes = []
    if albumin < 3.5:
        params["collagen_factor"] *= albumin / 3.5
        params["prolif_factor"] *= albumin / 4.0
        notes.append(f"Albumin {albumin} g/dL: collagen x{albumin/3.5:.2f}, "
                     f"prolif x{albumin/4.0:.2f} (malnourished)")
    if albumin < 3.0:
        # Severe malnutrition blanket penalty
        for key in ["prolif_factor", "collagen_factor", "migration_factor",
                     "angio_rate", "resolution_factor"]:
            params[key] *= 0.7
        notes.append(f"Albumin {albumin} g/dL: severe malnutrition, all factors x0.7")
    return notes


def apply_blood(params, hemoglobin, platelets, anticoagulant):
    """Blood parameters affect coagulation, O2 delivery, and bleeding risk."""
    notes = []

    # Hemoglobin: normalize to 14 g/dL baseline
    hb_norm = hemoglobin / 14.0
    params["hemoglobin"] = hb_norm
    if hemoglobin < 10.0:
        # Anemia impairs O2 delivery and healing
        params["perfusion_basal"] *= max(0.5, hb_norm)
        notes.append(f"Hemoglobin {hemoglobin} g/dL: anemia, perfusion x{max(0.5, hb_norm):.2f}")
    if hemoglobin < 7.0:
        # Severe anemia: critical O2 deficit
        params["prolif_factor"] *= 0.6
        params["collagen_factor"] *= 0.6
        notes.append(f"Hemoglobin {hemoglobin} g/dL: severe anemia, prolif and collagen x0.6")

    # Platelets: normalize to 250k/uL baseline
    plt_norm = platelets / 250.0
    params["platelet_count"] = plt_norm
    if platelets < 100:
        # Thrombocytopenia: impaired clotting
        params["bleed_rate"] *= 1.0 + (1.0 - plt_norm)
        notes.append(f"Platelets {platelets}k: thrombocytopenia, bleed rate x{1.0 + (1.0 - plt_norm):.2f}")
    if platelets < 50:
        # Severe: spontaneous bleeding risk
        params["bleed_rate"] *= 1.5
        notes.append(f"Platelets {platelets}k: severe thrombocytopenia, high bleed risk")

    # Anticoagulant therapy
    if anticoagulant == "prophylactic":
        params["anticoagulant_level"] = 0.5
        params["bleed_rate"] *= 1.3
        notes.append("Prophylactic anticoagulation: bleed rate x1.3")
    elif anticoagulant == "therapeutic":
        params["anticoagulant_level"] = 1.0
        params["bleed_rate"] *= 1.8
        notes.append("Therapeutic anticoagulation: bleed rate x1.8, monitor for hematoma")

    return notes


def apply_wbc(params, wbc):
    """WBC count affects immune capacity (Guo & DiPietro 2010)."""
    notes = []
    if wbc < 4.0:
        ratio = wbc / 7.0
        params["neutrophil_factor"] *= ratio
        notes.append(f"WBC {wbc}: leukopenia, neutrophil factor x{ratio:.2f}")
    if wbc > 11.0:
        params["inflammation_baseline"] += 0.0005
        notes.append(f"WBC {wbc}: leukocytosis, inflammation baseline +0.0005")
    return notes


# ---------------------------------------------------------------------------
# TOML generation
# ---------------------------------------------------------------------------

def generate_toml(args, params, notes):
    """Build profile TOML string with documentation header."""
    lines = []

    # Header comment block summarizing patient and clinical factors
    lines.append(f"# Profile: patient (clinically calibrated)")
    lines.append(f"# Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    lines.append(f"# Tool: scripts/clinical_calibrate.py")
    lines.append(f"#")
    lines.append(f"# Patient parameters:")
    lines.append(f"#   Age: {args.age} years | BMI: {args.bmi} kg/m2 | HbA1c: {args.hba1c}%")
    lines.append(f"#   Wound: {args.wound_type}, {args.wound_size}mm, {args.wound_depth} thickness")
    lines.append(f"#   Location: {args.location} | Albumin: {args.albumin} g/dL | WBC: {args.wbc}")

    flags = []
    if args.smoking:
        flags.append("smoker")
    if args.immunosuppressed:
        flags.append("immunosuppressed")
    if args.infection != "none":
        flags.append(f"infection:{args.infection}")
    if flags:
        lines.append(f"#   Flags: {', '.join(flags)}")

    lines.append(f"#")
    lines.append(f"# Calibration factors applied:")
    for note in notes:
        lines.append(f"#   {note}")
    lines.append("")

    # [skin] section
    lines.append("[skin]")
    lines.append(f"stem_fraction = {params['stem_fraction']:.3f}")
    if params["g1_duration"] != DEFAULTS["g1_duration"]:
        lines.append(f"g1_duration = {params['g1_duration']:.1f}")
    if params["s_duration"] != DEFAULTS["s_duration"]:
        lines.append(f"s_duration = {params['s_duration']:.1f}")
    if params["water_recovery_rate"] != DEFAULTS["water_recovery_rate"]:
        lines.append(f"water_recovery_rate = {params['water_recovery_rate']:.3f}")
    lines.append("")

    # [skin.perfusion] section
    lines.append("[skin.perfusion]")
    lines.append(f"basal = {params['perfusion_basal']:.3f}")
    lines.append(f"angio_rate = {params['angio_rate']:.4f}")
    lines.append("")

    # [skin.diabetic] section (only if diabetic mode triggered)
    if params.get("diabetic_mode"):
        lines.append("[skin.diabetic]")
        lines.append("mode = true")
        lines.append(f"m1_duration_factor = {params['m1_duration_factor']:.1f}")
        lines.append(f"resolution_factor = {params['resolution_factor']:.2f}")
        lines.append(f"efferocytosis_factor = {params['efferocytosis_factor']:.2f}")
        lines.append(f"neutrophil_factor = {params['neutrophil_factor']:.1f}")
        lines.append(f"prolif_factor = {params['prolif_factor']:.2f}")
        lines.append(f"migration_factor = {params['migration_factor']:.2f}")
        lines.append(f"collagen_factor = {params['collagen_factor']:.2f}")
        lines.append(f"fibroblast_activation_factor = {params['fibroblast_activation_factor']:.1f}")
        lines.append(f"fibroblast_lifespan_factor = {params['fibroblast_lifespan_factor']:.2f}")
        lines.append(f"baseline_inflammation = {params['inflammation_baseline']:.4f}")
        lines.append(f"inflammation_sensitivity = {params['inflammation_sensitivity']:.1f}")
        lines.append(f"vegf_factor = {params['vegf_factor']:.2f}")
        lines.append(f"mmp_factor = {params['mmp_factor']:.1f}")
        lines.append("")
    elif params["m1_duration_factor"] != 1.0:
        # Pre-diabetic range: some immune impairment without full diabetic mode
        lines.append("[skin.diabetic]")
        lines.append("mode = false")
        lines.append(f"m1_duration_factor = {params['m1_duration_factor']:.1f}")
        lines.append(f"prolif_factor = {params['prolif_factor']:.2f}")
        lines.append(f"collagen_factor = {params['collagen_factor']:.2f}")
        lines.append(f"inflammation_sensitivity = {params['inflammation_sensitivity']:.1f}")
        lines.append("")

    # [skin.blood] section
    blood_changed = (params["hemoglobin"] != 1.0 or params["anticoagulant_level"] != 0.0
                     or params["platelet_count"] != 1.0 or params["bleed_rate"] != 0.02)
    if blood_changed:
        lines.append("[skin.blood]")
        lines.append(f"hemoglobin = {params['hemoglobin']:.2f}")
        lines.append(f"platelet_count = {params['platelet_count']:.2f}")
        lines.append(f"bleed_rate = {params['bleed_rate']:.4f}")
        if params["anticoagulant_level"] > 0:
            lines.append(f"anticoagulant_level = {params['anticoagulant_level']:.1f}")
        lines.append("")

    # Biofilm section if infection present
    if params["biofilm_seed_delay_h"] is not None:
        lines.append("[skin.biofilm]")
        lines.append(f"seed_delay_h = {params['biofilm_seed_delay_h']}")
        lines.append(f"growth_rate = {params['biofilm_growth_rate']:.4f}")
        lines.append("")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Treatment recommendations
# ---------------------------------------------------------------------------

def get_recommendations(args):
    """Return evidence-based clinical recommendations for the patient profile."""
    recs = []

    if args.hba1c > 7.0:
        recs.append("Consider glycemic control before elective surgery. "
                     "HbA1c >7% is associated with 2 to 3x wound complication risk "
                     "(Endara et al. 2013).")
    elif args.hba1c > 6.5:
        recs.append("Monitor blood glucose closely. Pre-diabetic HbA1c levels "
                     "may impair wound healing.")

    if args.albumin < 3.5:
        recs.append("Nutritional supplementation recommended (protein, zinc, "
                     "vitamin C). Hypoalbuminemia is a strong predictor of "
                     "delayed healing (Stechmiller 2010).")
    if args.albumin < 3.0:
        recs.append("Severe malnutrition detected. Consider enteral or parenteral "
                     "nutrition support before surgical intervention.")

    if args.smoking:
        recs.append("Smoking cessation advised. Smokers have 2 to 4x higher wound "
                     "complication rates (Sorensen 2012). Minimum 4 weeks abstinence "
                     "before elective procedures.")

    if args.bmi > 35:
        recs.append("Bariatric optimization may improve wound outcomes. Morbid "
                     "obesity impairs perfusion and increases surgical site infection "
                     "risk (Pierpont et al. 2014).")
    elif args.bmi > 30:
        recs.append("Obesity noted. Monitor for signs of impaired healing and "
                     "consider weight management.")

    if args.infection == "severe":
        recs.append("Systemic antibiotics indicated. Consider wound cultures and "
                     "sensitivity testing before empiric therapy. Surgical "
                     "debridement may be required (Bjarnsholt 2013).")
    elif args.infection == "moderate":
        recs.append("Topical antimicrobial therapy recommended. Monitor for "
                     "progression to systemic infection.")

    if args.wbc < 4.0:
        recs.append("Leukopenia detected. Evaluate for underlying cause. "
                     "Impaired neutrophil recruitment may delay inflammatory phase.")
    if args.wbc > 11.0:
        recs.append("Leukocytosis detected. Evaluate for systemic infection or "
                     "chronic inflammatory condition.")

    if args.hemoglobin < 10.0:
        recs.append("Anemia detected. Consider iron supplementation or transfusion "
                     "if hemoglobin <7 g/dL. Anemia impairs oxygen delivery to the "
                     "wound bed (Hopf et al. 1997).")
    if args.hemoglobin < 7.0:
        recs.append("Severe anemia: transfusion likely indicated. Tissue hypoxia "
                     "severely compromises all phases of wound healing.")

    if args.platelets < 100:
        recs.append(f"Thrombocytopenia ({args.platelets}k/uL). Evaluate underlying "
                     "cause. Impaired hemostasis may lead to wound hematoma and "
                     "delayed healing.")
    if args.platelets < 50:
        recs.append("Severe thrombocytopenia: high spontaneous bleeding risk. "
                     "Consider platelet transfusion before surgical intervention.")

    if args.anticoagulant == "therapeutic":
        recs.append("Therapeutic anticoagulation increases bleeding and hematoma risk. "
                     "Coordinate with prescribing physician about bridging strategy "
                     "for elective procedures (Douketis et al. 2015).")

    if args.age > 75:
        recs.append("Advanced age increases healing time significantly. Consider "
                     "negative pressure wound therapy for complex wounds.")

    if args.immunosuppressed:
        recs.append("Immunosuppression may delay all phases of wound healing. "
                     "Coordinate with prescribing physician regarding dose "
                     "adjustment during healing window.")

    if args.location in ("sacrum", "heel"):
        recs.append(f"Wound at {args.location} is at risk for pressure injury. "
                     "Implement offloading protocol and regular repositioning.")

    return recs


# ---------------------------------------------------------------------------
# Summary printing
# ---------------------------------------------------------------------------

def print_summary(args, params, notes):
    """Print a concise summary of the generated profile to stdout."""
    print("\n=== Clinical Calibration Summary ===\n")
    print(f"Patient: {args.age}y, BMI {args.bmi}, HbA1c {args.hba1c}%")
    print(f"Wound:   {args.wound_type} at {args.location}, "
          f"{args.wound_size}mm, {args.wound_depth}")

    diabetic_status = "diabetic" if params.get("diabetic_mode") else (
        "pre-diabetic" if args.hba1c > 5.7 else "normal")
    print(f"Status:  glycemic={diabetic_status}, albumin={args.albumin} g/dL, "
          f"WBC={args.wbc}")

    print("\nKey simulation parameters:")
    print(f"  Perfusion basal:      {params['perfusion_basal']:.3f}")
    print(f"  Angiogenesis rate:    {params['angio_rate']:.4f}")
    print(f"  Stem fraction:        {params['stem_fraction']:.3f}")
    print(f"  Proliferation factor: {params['prolif_factor']:.3f}")
    print(f"  Collagen factor:      {params['collagen_factor']:.3f}")
    print(f"  M1 duration factor:   {params['m1_duration_factor']:.2f}")
    if params["hemoglobin"] != 1.0:
        print(f"  Hemoglobin:           {params['hemoglobin']:.2f} (norm)")
    if params["platelet_count"] != 1.0:
        print(f"  Platelet count:       {params['platelet_count']:.2f} (norm)")
    if params["anticoagulant_level"] > 0:
        print(f"  Anticoagulant:        {params['anticoagulant_level']:.1f}")
    if params["bleed_rate"] != 0.02:
        print(f"  Bleed rate:           {params['bleed_rate']:.4f}")
    if params["biofilm_seed_delay_h"] is not None:
        print(f"  Biofilm delay:        {params['biofilm_seed_delay_h']}h")
        print(f"  Biofilm growth:       {params['biofilm_growth_rate']:.4f}")

    print(f"\n{len(notes)} calibration adjustments applied.")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    # Start from healthy baseline
    params = dict(DEFAULTS)
    all_notes = []

    # Apply each clinical factor in order. Each function returns
    # annotation strings and modifies params in place.
    all_notes.extend(apply_age(params, args.age))
    all_notes.extend(apply_bmi(params, args.bmi))
    all_notes.extend(apply_hba1c(params, args.hba1c))
    if args.smoking:
        all_notes.extend(apply_smoking(params))
    if args.immunosuppressed:
        all_notes.extend(apply_immunosuppression(params))
    all_notes.extend(apply_infection(params, args.infection))
    all_notes.extend(apply_location(params, args.location))
    all_notes.extend(apply_blood(params, args.hemoglobin, args.platelets, args.anticoagulant))
    all_notes.extend(apply_albumin(params, args.albumin))
    all_notes.extend(apply_wbc(params, args.wbc))

    # Generate TOML content
    toml_content = generate_toml(args, params, all_notes)

    # Write output file
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(toml_content)
    print(f"Profile written to: {output_path}")

    # Print summary
    print_summary(args, params, all_notes)

    # Treatment recommendations
    if args.recommend:
        recs = get_recommendations(args)
        if recs:
            print("\n=== Treatment Recommendations ===\n")
            for i, rec in enumerate(recs, 1):
                print(f"  {i}. {rec}")
        else:
            print("\nNo specific treatment recommendations for this patient profile.")

    print()


if __name__ == "__main__":
    main()
