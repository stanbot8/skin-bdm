# Rheumatoid Arthritis Study

RA-driven synovial inflammation, cartilage erosion, and bone destruction overlaid on wound healing. Uses a study-scoped module with TNF-alpha/IL-6 dual cytokine axis, T cell adaptive immunity, sigmoid flare dynamics, pannus hyperplasia, subchondral bone erosion, drug pharmacokinetics, and five treatment strategies.

## Running

```bash
# Baseline RA (untreated, 10-run consensus)
python3 batch/batch.py -n 10 --skin rheumatoid --study rheumatoid

# With anti-TNF treatment
python3 batch/batch.py -n 10 --skin rheumatoid --study rheumatoid --treatment anti_tnf

# Full treatment comparison study (baseline + 5 treatments)
python3 studies/rheumatoid/study.py

# Quick 5-run consensus
python3 studies/rheumatoid/study.py --quick

# Run experiments
python3 studies/run_experiments.py studies/rheumatoid/experiments/biologic_comparison.toml
python3 studies/run_experiments.py studies/rheumatoid/experiments/disease_severity.toml
python3 studies/run_experiments.py studies/rheumatoid/experiments/treatment_timing.toml
python3 studies/run_experiments.py studies/rheumatoid/experiments/cartilage_protection.toml
```

## Configuration

- Duration: 30 days
- Wound: enabled (t=0)
- RA module: enabled (study-scoped)
- Profile: `profiles/rheumatoid.toml` (M1 prolongation, chronic recruitment, FLS hyperplasia)

## Study Module

The RA module lives in `modules/rheumatoid/` (study-scoped, not engine-level):

| File | Description |
|---|---|
| tnf_alpha_pde.h | TNF-alpha cytokine field (17 kDa homotrimer) |
| il6_pde.h | IL-6 cytokine field (JAK/STAT3 axis) |
| cartilage_pde.h | Cartilage integrity (structural target) |
| synovial_fluid_pde.h | Synovial fluid / pannus tissue density |
| tcell_density_pde.h | T cell density (adaptive immunity proxy, Th1/Th17 infiltration) |
| bone_pde.h | Subchondral bone integrity (RANKL/osteoclast erosion target) |
| source_hook.h | Autoimmune + T cell + immune cell cytokine production, pannus growth, drug PK, treatment clearance |
| post_hook.h | TNF/IL-6 downstream effects, pannus-amplified cartilage and bone erosion |
| config.toml | All RA parameters with literature DOIs |

## Treatments

| Treatment | File | Mechanism |
|---|---|---|
| Anti-TNF | [anti_tnf.toml](treatments/anti_tnf.toml) | TNF-alpha neutralization (infliximab/adalimumab/etanercept) |
| Tocilizumab | [tocilizumab.toml](treatments/tocilizumab.toml) | IL-6R blockade (JAK/STAT3 suppression) |
| Methotrexate | [methotrexate.toml](treatments/methotrexate.toml) | Folate antagonist DMARD (broad immunosuppression, T cell proliferation block) |
| JAK inhibitor | [jak_inhibitor.toml](treatments/jak_inhibitor.toml) | JAK1/JAK3 kinase inhibition (tofacitinib/baricitinib, strongest T cell suppression) |
| Triple combination | [ra_combination.toml](treatments/ra_combination.toml) | Anti-TNF + tocilizumab + methotrexate |

## Experiments

| Experiment | File | Configs | Runs |
|---|---|---|---|
| Biologic comparison | [biologic_comparison.toml](experiments/biologic_comparison.toml) | 6 (untreated + 5 treatments) | 5/config |
| Disease severity | [disease_severity.toml](experiments/disease_severity.toml) | 5 (0.5x to 2.0x flare intensity) | 5/config |
| Treatment timing | [treatment_timing.toml](experiments/treatment_timing.toml) | 12 (anti-TNF + tocilizumab + JAK inhibitor x 4 windows) | 5/config |
| Cartilage and bone protection | [cartilage_protection.toml](experiments/cartilage_protection.toml) | 7 (untreated + 6 strategies) | 5/config |

## Validation

Three RA observables validated against literature reference curves:

| Observable | Normalization | Source |
|---|---|---|
| TNF-alpha | Peak-normalized | Feldmann & Maini 2003, McInnes & Schett 2011 |
| IL-6 | Peak-normalized | Kishimoto 2005, McInnes & Schett 2011 |
| Cartilage integrity | Absolute (1.0 = healthy) | McInnes & Schett 2011, Firestein 2003, Smolen 2018 |

Reference data in `modules/rheumatoid/data/`. Run standalone:

```bash
python3 literature/validators/compare.py ra studies/rheumatoid/results/metrics.csv
```

## Tests

Study-scoped tests in `tests/ra_test.cc`: field registration, cartilage initialization, sigmoid flare dynamics, pannus fibroblast boost, M1 prolongation, anti-TNF clearance.
