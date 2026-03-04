# Surgical Site Infection (SSI) Study

Surgical wound healing and infection modeling. SSIs complicate 2 to 5 percent of all surgeries and are a leading cause of hospital readmission, morbidity, and death. This study models the CDC wound classification system (clean, clean-contaminated, contaminated, dirty) and tests evidence-based SSI prevention bundles.

## Configuration

- Duration: 30 days
- Wound: enabled (t=0)
- pH, fibroblasts, hemostasis, biofilm: enabled
- Profile: `profiles/surgical.toml` (clean incision, approximated edges)

## Treatments

prophylactic_antibiotics, chlorhexidine, negative_pressure, enhanced_recovery, antimicrobial_suture

## Experiments (4 defined)

| Experiment | Description |
|---|---|
| infection_risk | CDC wound classes (clean through dirty) |
| prevention_bundle | Additive SSI prevention measures |
| comorbidity_risk | Patient risk factors (obesity, diabetes, immunosuppression, age) |
| timing_optimization | Antibiotic prophylaxis timing window |

## Usage

```bash
# Single surgical wound (no infection)
python3 batch/batch.py -n 10 --study surgical --skin surgical

# With prophylactic antibiotics
python3 batch/batch.py -n 10 --study surgical --skin surgical --treatment prophylactic_antibiotics

# Run a specific experiment
python3 scripts/study/experiment_runner.py studies/surgical/experiments/infection_risk.toml
```
