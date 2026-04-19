# Skin BioDynaMo (SkiBiDy)

**SkiBiDy** is a hybrid agent-continuum simulation of skin tissue biology built on [BioDynaMo](https://biodynamo.org). It models wound healing, immune response, fibroblast and collagen dynamics, scar formation, vascular perfusion, hemostasis, cellular senescence, cutaneous innervation, tumor growth, diabetic impairment, and therapeutic interventions, all from mechanistic first principles backed by DOI-linked source papers across 34 modules.

![Wound healing simulation](docs/skibidy.gif?v=2)

## How it works

Healthy skin runs as a composite field coupling 47 diffusion grids with no agents. When an event occurs (wound, infection, tumor), cells spawn from local field state, interact with the fields, and dissolve back once stable. This is the **UWYN** (Use What You Need) paradigm: the simulation only creates agents where the biology demands cellular resolution.

```
Corneum    [continuum]  barrier, desquamation
Granulosum [continuum]  keratohyalin, tight junctions
Spinosum   [continuum]  desmosomes; agents on event
Basale     [continuum]  stem/TA cells on event
---------- basement membrane ----------
Dermis     [continuum]  vasculature, O2/KGF, collagen, nerves; agents on event
```

## Quick start

```bash
source <path_to_biodynamo>/bin/thisbdm.sh
./run.sh                                     # interactive menu
./run.sh --study=wound                       # punch biopsy wound healing
./run.sh --study=diabetic-wound              # chronic diabetic ulcer
./run.sh --study=tumor                       # basal cell carcinoma growth
./run.sh --compare                           # normal vs diabetic side-by-side
./tests/test.sh                              # unit test suite (540 tests)
```

## Modules

34 modules, each self-contained in `modules/` with its own config, source, data, and README. See [docs/](docs/README.md) for the full module index with biology, parameters, coupling, and validation.

| Module | Default | Module | Default |
|--------|---------|--------|---------|
| [tissue](modules/tissue/) | on | [dermis](modules/dermis/) | on |
| [wound](modules/wound/) | on | [elastin](modules/elastin/) | off |
| [immune](modules/immune/) | on | [hyaluronan](modules/hyaluronan/) | off |
| [inflammation](modules/inflammation/) | on | [glucose](modules/glucose/) | on |
| [fibroblast](modules/fibroblast/) | on | [temperature](modules/temperature/) | on |
| [scar](modules/scar/) | on | [lactate](modules/lactate/) | on |
| [mmp](modules/mmp/) | on | [nitric_oxide](modules/nitric_oxide/) | on |
| [fibronectin](modules/fibronectin/) | on | [ph](modules/ph/) | off |
| [angiogenesis](modules/angiogenesis/) | on | [hemostasis](modules/hemostasis/) | off |
| [perfusion](modules/perfusion/) | on | [biofilm](modules/biofilm/) | off |
| [diabetic](modules/diabetic/) | off | [tumor](modules/tumor/) | off |
| [senescence](modules/senescence/) | on | [neuropathy](modules/neuropathy/) | on |
| [ros](modules/ros/) | on | [bioelectric](modules/bioelectric/) | on |
| [lymphatic](modules/lymphatic/) | on | [mechanotransduction](modules/mechanotransduction/) | on |
| [blood](modules/blood/) | off | [burn](modules/burn/) | off |
| [pressure](modules/pressure/) | off | [scab](modules/scab/) | on |
| [photon](modules/photon/) | off | [body_site](modules/body_site/) | off |

## Configuration

Config is layered TOML, merged at runtime:

```
bdm.core.toml              core tissue params
  + modules/*/config.toml   34 module configs (auto-merged)
  + profiles/*.toml          skin phenotype overlay
  + studies/*/preset.toml    study experiment overlay
  + studies/*/treatments/*.toml  therapeutic intervention overlay
  = bdm.toml                 runtime config (gitignored)
```

Skin profiles: `normal` (default), `aged`, `diabetic`, `aged_diabetic`, `keloid`, `hypertrophic`, `venous`, `scleroderma`, `psoriasis`, `rheumatoid`, `burn`, `pressure`, `surgical`, `tumor`, `tumor_wound`. Studies: `wound`, `diabetic-wound`, `tumor`, `tumor-wound`, `burn`, `pressure-ulcer`, `surgical`, `rheumatoid`, `baseline`, `full-model`.

## Treatment study

The diabetic treatment study compares 9 therapeutic interventions against an untreated diabetic baseline:

```bash
# Single treatment schedules (NPWT/HBO/MSC at different start days)
python3 scripts/study/gen_treatment_schedule.py --treatments npwt,hbo,msc --days 0,7,14,21

# Pairwise combinatorial scheduling
python3 scripts/study/gen_treatment_schedule.py --treatments npwt,hbo,msc --days 0,7,14,21 --combo

# Run all generated schedules
python3 studies/run_experiments.py studies/diabetic-wound/experiments/schedule_*.toml

# Full comparison study (all 9 treatments + baseline)
python3 studies/diabetic-wound/treatment.py --combos all
```

Available treatments: anti-inflammatory, HBO, NPWT, doxycycline, growth factor (PDGF-BB), MSC therapy, moisture dressings, senolytic (dasatinib+quercetin), and a rational combination targeting all four dysfunction axes. See [docs/treatments.md](docs/treatments.md) for mechanisms and references.

## Rheumatoid arthritis study

The RA study models synovial inflammation and cartilage erosion with a study-scoped TNF-alpha/IL-6 dual cytokine module:

```bash
# RA baseline (untreated, 10-run consensus)
python3 batch/batch.py -n 10 --skin rheumatoid --study rheumatoid

# Treatment comparison (baseline + 5 treatments)
python3 studies/rheumatoid/study.py

# Biologic comparison experiment
python3 studies/run_experiments.py studies/rheumatoid/experiments/biologic_comparison.toml
```

Available RA treatments: anti-TNF (infliximab/adalimumab), tocilizumab (anti-IL-6R), methotrexate (DMARD), JAK inhibitor, and triple combination. See [studies/rheumatoid/](studies/rheumatoid/) for details.

## Additional studies

| Study | Profile | Treatments | Experiments | Description |
|-------|---------|------------|-------------|-------------|
| [Burn](studies/burn/) | `burn` | 5 (cooling, debridement, silver sulfadiazine, skin substitute, pressure garment) | 4 | Jackson's three-zone thermal injury, stasis rescue, scar prevention |
| [Pressure ulcer](studies/pressure-ulcer/) | `pressure` | 5 (redistribution, wound VAC, nutrition, silver dressing, offloading) | 4 | Ischemia-reperfusion injury, stage I through IV, moisture/shear co-factors |
| [Surgical](studies/surgical/) | `surgical` | 5 (prophylactic antibiotics, chlorhexidine, negative pressure, enhanced recovery, antimicrobial suture) | 4 | SSI prevention bundles, CDC wound classes, comorbidity risk |

Shared treatments (usable across all studies) live in `studies/shared/treatments/`: triple antibiotic ointment.

## Validation

Validated against published literature across 12 observables with per-module `SOURCES.yaml` files providing DOI-linked citations. Four optional mechanistic toggles (gradient-driven recruitment, efferocytosis M1 to M2, constitutive + TGF-beta collagen, HIF-1alpha VEGF) can replace parametric models for biophysical fidelity testing.

```bash
python3 batch/batch.py -n 10 --study wound --validate          # parametric (default)
python3 batch/batch.py -n 10 --study full-model --validate      # mechanistic toggles
```

## Documentation

| Document | Description |
|----------|-------------|
| [Docs index](docs/README.md) | Module index, metrics columns, validation overview |
| [Guide](docs/guide.md) | Configuration, skin profiles, studies, running, visualization |
| [Architecture](docs/architecture.md) | UWYN hybrid agent-continuum design, CompositeField, field coupling |
| [Parameters](docs/parameters.md) | Parameter index with module links and config layering |
| [Treatments](docs/treatments.md) | Therapeutic interventions across 6 studies |
| [Studies](studies/README.md) | Packaged studies, experiments, example outputs |
| [Batch](batch/README.md) | Multi-run consensus, parameter sweeps |
| [Literature](literature/README.md) | Validation framework, RMSE dashboard, reference data |

## Third-party

| Library | Version | License |
|---------|---------|---------|
| [toml++](https://github.com/marzer/tomlplusplus) | 3.4.0 | [MIT](third_party/tomlplusplus/LICENSE) |
| [BioDynaMo](https://biodynamo.org) | 1.04+ (master) | Apache 2.0 (linked, not bundled) |

## License

[Apache 2.0](LICENSE)
