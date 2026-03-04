# Pressure Ulcer Study

Pressure ulcer (bedsore) development, staging, prevention, and treatment overlaid on wound healing. Models ischemia-reperfusion injury from sustained mechanical load exceeding capillary closing pressure, with stage I through IV progression, moisture and shear co-factors, and repositioning protocols.

Pressure ulcers are the leading preventable hospital-acquired injury, affecting up to 2.5 million patients annually in the US. The simulation captures the core pathophysiology: tissue compression causes local ischemia, and subsequent reperfusion generates reactive oxygen species that amplify cell damage and inflammation.

## Running

```bash
# Baseline pressure ulcer (untreated, 10-run consensus)
python3 batch/batch.py -n 10 --skin pressure --study pressure-ulcer

# With pressure redistribution treatment
python3 batch/batch.py -n 10 --skin pressure --study pressure-ulcer --treatment pressure_redistribution

# Run experiments
python3 studies/run_experiments.py studies/pressure-ulcer/experiments/stage_progression.toml
python3 studies/run_experiments.py studies/pressure-ulcer/experiments/prevention_vs_treatment.toml
python3 studies/run_experiments.py studies/pressure-ulcer/experiments/treatment_comparison.toml
python3 studies/run_experiments.py studies/pressure-ulcer/experiments/moisture_shear.toml
```

## Configuration

- Duration: 42 days (6 weeks, typical healing window for Stage II/III)
- Wound: enabled (t=0)
- Pressure module: enabled (ischemia-reperfusion injury)
- Profile: `profiles/pressure.toml` (immobilized elderly patient, immunosenescence, reduced perfusion)

## Treatments

| Treatment | File | Mechanism |
|---|---|---|
| Pressure redistribution | [pressure_redistribution.toml](treatments/pressure_redistribution.toml) | Low-air-loss mattress with 2h repositioning protocol |
| Wound VAC | [wound_vac.toml](treatments/wound_vac.toml) | Negative pressure wound therapy for Stage III/IV |
| Nutrition | [nutrition.toml](treatments/nutrition.toml) | Protein, zinc, vitamin C supplementation |
| Silver dressing | [silver_dressing.toml](treatments/silver_dressing.toml) | Antimicrobial silver with moisture balance |
| Offloading | [offloading_protocol.toml](treatments/offloading_protocol.toml) | Complete pressure elimination (gold standard) |

## Experiments

| Experiment | File | Configs | Runs |
|---|---|---|---|
| Stage progression | [stage_progression.toml](experiments/stage_progression.toml) | 4 (Stage I through IV) | 5/config |
| Prevention vs treatment | [prevention_vs_treatment.toml](experiments/prevention_vs_treatment.toml) | 3 (none, early, late) | 5/config |
| Treatment comparison | [treatment_comparison.toml](experiments/treatment_comparison.toml) | 8 (untreated + 5 mono + 2 combos + all) | 5/config |
| Moisture and shear | [moisture_shear.toml](experiments/moisture_shear.toml) | 4 (2x2 factorial) | 5/config |
