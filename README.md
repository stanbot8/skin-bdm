# Skin BioDynaMo (SkiBiDy)

**SkiBiDy** is a hybrid agent-continuum simulation of skin tissue biology built on [BioDynaMo](https://biodynamo.org). It models wound healing, immune response, fibroblast and collagen dynamics, scar formation, vascular perfusion, hemostasis, cellular senescence, cutaneous innervation, tumor growth, diabetic impairment, and therapeutic interventions, all from mechanistic first principles backed by 186 DOI-linked source papers.

![Wound healing simulation](docs/skibidy.gif?v=2)

## How it works

Healthy skin runs as a composite field coupling 33 PDEs with no agents. When an event occurs (wound, infection, tumor), cells spawn from local field state, interact with the fields, and dissolve back once stable. This is the **UWYN** (Use What You Need) paradigm: the simulation only creates agents where the biology demands cellular resolution.

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
./tests/test.sh                              # 130 unit tests
```

## Modules

24 modules, each self-contained in `modules/` with its own config, source, data, and README. See [docs/](docs/README.md) for the full module index with biology, parameters, coupling, and validation.

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

## Configuration

Config is layered TOML, merged at runtime:

```
bdm.core.toml              core tissue params
  + modules/*/config.toml   24 module configs (auto-merged)
  + profiles/*.toml          skin phenotype overlay
  + studies/*/preset.toml    study scenario overlay
  + treatments/*.toml        therapeutic intervention overlay
  = bdm.toml                 runtime config (gitignored)
```

Skin profiles: `normal` (default), `aged`, `diabetic`, `aged_diabetic`. Studies: `wound`, `diabetic-wound`, `tumor`, `tumor-wound`, `baseline`.

## Treatment study

The diabetic treatment study compares 8 therapeutic interventions against an untreated diabetic baseline:

```bash
# Single treatment schedules (NPWT/HBO/MSC at different start days)
python3 scripts/study/gen_treatment_schedule.py --treatments npwt,hbo,msc --days 0,7,14,21

# Pairwise combinatorial scheduling
python3 scripts/study/gen_treatment_schedule.py --treatments npwt,hbo,msc --days 0,7,14,21 --combo

# Run all generated schedules
./studies/run-scenarios.sh studies/diabetic-wound/scenarios/schedule_*.toml

# Full comparison study (all 8 treatments + baseline)
./studies/diabetic-wound/treatment.sh --combos=all
```

Available treatments: anti-inflammatory, HBO, NPWT, doxycycline, growth factor (PDGF-BB), MSC therapy, moisture dressings, and a rational combination targeting all four dysfunction axes. See [docs/treatments.md](docs/treatments.md) for mechanisms and references.


## Validation

Validated against published literature across 11 observables, backed by 191 DOI-linked source papers in per-module `SOURCES.yaml` files. Run `python3 literature/validate_all.py` for the RMSE dashboard or `python3 batch/batch.py -n 10 --study wound --validate` for a 10-run consensus with literature comparison.

## Documentation

| Document | Description |
|----------|-------------|
| [Docs index](docs/README.md) | Module index, metrics columns, validation overview |
| [Guide](docs/guide.md) | Configuration, skin profiles, studies, running, visualization |
| [Architecture](docs/architecture.md) | UWYN hybrid agent-continuum design, CompositeField, field coupling |
| [Parameters](docs/parameters.md) | Parameter index with module links and config layering |
| [Treatments](docs/treatments.md) | 8 therapeutic interventions for diabetic wounds |
| [Studies](studies/README.md) | Packaged experiments, scenarios, example outputs |
| [Batch](batch/README.md) | Multi-run consensus, parameter sweeps |
| [Literature](literature/README.md) | Validation framework, RMSE dashboard, reference data |

## Third-party

| Library | Version | License |
|---------|---------|---------|
| [toml++](https://github.com/marzer/tomlplusplus) | 3.4.0 | [MIT](third_party/tomlplusplus/LICENSE) |
| [BioDynaMo](https://biodynamo.org) | 1.04+ (master) | Apache 2.0 (linked, not bundled) |

## License

[Apache 2.0](LICENSE)
