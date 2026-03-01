# Skibidy

**Skin BioDynaMo (SkiBiDy)** is a hybrid agent-continuum simulation of skin tissue biology built on [BioDynaMo](https://biodynamo.org). Models wound healing, immune response, fibroblast/collagen dynamics, scar formation, vascular perfusion, hemostasis, tumor growth, and diabetic impairment.

![Wound healing simulation](docs/skibidy.gif?v=2)

## Overview

Healthy skin runs as a composite field coupling 23 PDEs with no agents. When an event occurs (wound, infection, tumor), cells spawn from local field state, interact with the fields, and dissolve back once stable. This is called **UWYN** (Use What You Need), meaning the simulation only creates agents where the biology demands resolution.

```
Corneum    [continuum]  barrier, desquamation
Granulosum [continuum]  keratohyalin, tight junctions
Spinosum   [continuum]  desmosomes; agents on event
Basale     [continuum]  stem/TA cells on event
---------- basement membrane ----------
Dermis     [continuum]  vasculature, O2/KGF, collagen; agents on event
```

## Quick start

```bash
source <path_to_biodynamo>/bin/thisbdm.sh
./run.sh                                     # interactive menu
./run.sh --preset=wound                      # punch biopsy wound healing
./run.sh --preset=diabetic_wound             # chronic diabetic ulcer
./run.sh --preset=tumor                      # basal cell carcinoma growth
./run.sh --compare                           # normal vs diabetic side-by-side
./tests/test.sh                              # 130 unit tests
```

## Modules

18 modules, each self-contained in `modules/` with its own config, source, data, and README. See [docs/](docs/README.md) for the module index with biology, parameters, coupling, and validation.

| Module | Default | Module | Default |
|--------|---------|--------|---------|
| [tissue](modules/tissue/) | on | [perfusion](modules/perfusion/) | on |
| [wound](modules/wound/) | on | [dermis](modules/dermis/) | on |
| [immune](modules/immune/) | on | [elastin](modules/elastin/) | off |
| [inflammation](modules/inflammation/) | on | [hyaluronan](modules/hyaluronan/) | off |
| [fibroblast](modules/fibroblast/) | on | [ph](modules/ph/) | off |
| [scar](modules/scar/) | on | [hemostasis](modules/hemostasis/) | off |
| [mmp](modules/mmp/) | on | [biofilm](modules/biofilm/) | off |
| [fibronectin](modules/fibronectin/) | on | [diabetic](modules/diabetic/) | off |
| [angiogenesis](modules/angiogenesis/) | on | [tumor](modules/tumor/) | off |

## Configuration

Config is layered TOML, merged at runtime:

```
bdm.core.toml              core tissue params
  + modules/*/config.toml   18 module configs (auto-merged)
  + profiles/*.toml          skin phenotype overlay
  + presets/*.toml           scenario overlay
  = bdm.toml                 runtime config (gitignored)
```

Skin profiles: `normal` (default), `aged`, `diabetic`, `aged_diabetic`. Scenario presets: `wound`, `tumor`, `diabetic_wound`, `tumor_wound`.

## Validation

Validated against published literature across 11 observables, backed by 157 DOI-linked source papers in per-module `SOURCES.yaml` files. Run `python3 literature/validate_all.py` for the RMSE dashboard or `python3 batch/batch.py -n 10 --preset wound --validate` for a 10-run consensus with literature comparison.

## Documentation

See [docs/](docs/README.md) for full documentation including architecture, configuration guide, wound healing model, diabetic impairment, treatments, and validation details. See [studies/](studies/README.md) for packaged experiments and [batch/](batch/README.md) for multi-run consensus and parameter sweeps.

## Third-party

| Library | Version | License |
|---------|---------|---------|
| [toml++](https://github.com/marzer/tomlplusplus) | 3.4.0 | [MIT](third_party/tomlplusplus/LICENSE) |
| [BioDynaMo](https://biodynamo.org) | 1.04+ (master) | Apache 2.0 (linked, not bundled) |

## License

[Apache 2.0](LICENSE)
