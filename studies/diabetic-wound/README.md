# Diabetic Wound Study

Chronic diabetic ulcer simulation with therapeutic intervention screening. Tests 9 treatments individually and in combination, with timing sensitivity and severity spectrum analyses.

## Configuration

- Duration: 35 days
- Wound: enabled (t=0)
- pH, fibroblasts, hemostasis: enabled
- Diabetic mode: active

## Treatments

npwt, hbo, growth_factor, doxycycline, anti_inflammatory, msc, moisture, combination, senolytic

## Experiments (15 defined)

| Experiment | Description |
|---|---|
| aged_diabetic | Comorbidity across 4 phenotypes |
| biofilm_infection | Early/late infection with doxycycline |
| chronic_wound | 90-day chronic wounds with biofilm rescue |
| immune_tuning | Sensitivity analysis on immune axes |
| severity_spectrum | Diabetic severity sweep (0.5x to 1.5x) |
| wound_size | Radius scaling (1.5mm to 8mm) |
| treatment_timing | Delayed onset for top 3 treatments |
| schedule_* | Timing combinations for NPWT, HBO, MSC |

## Usage

```bash
# Single diabetic wound
python3 batch/batch.py -n 10 --study diabetic-wound --skin diabetic

# Treatment comparison (all 9 treatments)
python3 studies/diabetic-wound/treatment.py

# Adaptive combination search
python3 studies/diabetic-wound/adaptive.py --quick

# Run a specific experiment
python3 scripts/study/experiment_runner.py studies/diabetic-wound/experiments/severity_spectrum.toml
```
