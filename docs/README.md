> [Home](../README.md) / Docs

# Skin BioDynaMo (SkiBiDy) Documentation

## Guides

| Document | Description |
|----------|-------------|
| [Guide](guide.md) | Configuration, skin profiles, presets, running, visualization |
| [Architecture](architecture.md) | UWYN hybrid agent-continuum design, CompositeField, field coupling |
| [Parameters](parameters.md) | Parameter index with module links and config layering |
| [Treatment Interventions](treatments.md) | 8 therapeutic interventions for diabetic wounds |
| [Module Development](architecture.md#how-to-add-a-module) | Step-by-step guide to adding new modules |
| [Batch & Sweeps](../batch/README.md) | Multi-run consensus, parameter sensitivity analysis |
| [Studies](../studies/README.md) | Packaged experiments, scenarios, example outputs |

## Module index

Each module is self-contained in `modules/` with its own `config.toml`, source files, validation data, and README.

### Core

| Module | Biology | Fields |
|--------|---------|--------|
| [tissue](../modules/tissue/README.md) | Epidermal homeostasis, cell cycle, calcium-driven differentiation | Calcium, KGF, O2, Water, Stratum |

### Wound healing cascade

| Module | Biology | Fields |
|--------|---------|--------|
| [wound](../modules/wound/README.md) | Punch biopsy event, margin cell spawning | Stratum (zeros) |
| [inflammation](../modules/inflammation/README.md) | Cytokine signaling, Hill-function gating | Inflammation, ImmunePressure |
| [immune](../modules/immune/README.md) | Neutrophils, macrophages, M1/M2 polarization | (agents) |
| [fibroblast](../modules/fibroblast/README.md) | TGF-beta cascade, myofibroblast differentiation, collagen deposition | TGF-beta, Collagen |
| [mmp](../modules/mmp/README.md) | Matrix metalloproteinase ECM remodeling | MMP |
| [fibronectin](../modules/fibronectin/README.md) | Provisional matrix scaffold for migration | Fibronectin |
| [scar](../modules/scar/README.md) | Emergent scar formation from collagen | Scar |
| [perfusion](../modules/perfusion/README.md) | Vascular perfusion, angiogenesis recovery | Vascular |
| [angiogenesis](../modules/angiogenesis/README.md) | VEGF-driven vessel sprouting | VEGF |
| [dermis](../modules/dermis/README.md) | Dermal tissue integrity, sub-layer profile | Dermis |

### Pathology modifiers

| Module | Biology | Fields |
|--------|---------|--------|
| [diabetic](../modules/diabetic/README.md) | Chronic wound impairment across all systems | (modifier) |
| [biofilm](../modules/biofilm/README.md) | Bacterial biofilm colonization, PAMP feedback | Biofilm |

### Neoplasia

| Module | Biology | Fields |
|--------|---------|--------|
| [tumor](../modules/tumor/README.md) | Basal cell carcinoma, soft contact inhibition | Tumor |

### ECM groundwork

| Module | Biology | Fields |
|--------|---------|--------|
| [elastin](../modules/elastin/README.md) | Elastic fiber network, slow turnover | Elastin |
| [hyaluronan](../modules/hyaluronan/README.md) | Hyaluronic acid, water retention, migration scaffold | Hyaluronan |

### Environment

| Module | Biology | Fields |
|--------|---------|--------|
| [ph](../modules/ph/README.md) | Wound bed pH gradient, acid mantle disruption | pH |
| [hemostasis](../modules/hemostasis/README.md) | Fibrin clot scaffold, platelet activation | Fibrin |

## Validation

The validation framework lives in [`literature/`](../literature/) with its own [README](../literature/README.md). Reference data is collocated with each module under `modules/<module>/data/`, with citations in each module's `SOURCES.yaml`.

### Metrics columns

The CSV written to `output/skibidy/metrics.csv` has 35 columns:

| Column | Units | Description | Non-zero when |
|--------|-------|-------------|---------------|
| `step` | - | Simulation step | always |
| `time_h` | hours | Simulated time (step * dt) | always |
| `time_days` | days | Simulated time in days (time_h / 24) | always |
| `n_agents` | count | Active keratinocyte agents | agents spawned |
| `n_stem` | count | Stem cell agents | agents spawned |
| `n_ta` | count | Transit-amplifying agents (cycling, non-stem) | agents spawned |
| `n_g0` | count | Quiescent (G0) agents | agents spawned |
| `n_basal` | count | Agents in stratum basale | agents spawned |
| `n_spinous` | count | Agents in stratum spinosum | agents spawned |
| `n_granular` | count | Agents in stratum granulosum | agents spawned |
| `n_cornified` | count | Agents in stratum corneum | agents spawned |
| `n_neutrophils` | count | Active neutrophil agents | `[skin.wound] enabled` |
| `n_macrophages` | count | Active macrophage agents | `[skin.wound] enabled` |
| `wound_closure_pct` | % | Fraction of wound voxels with Stratum > 0.5 | `[skin.wound] enabled` |
| `mean_o2_wound` | normalized | Mean O2 in wound cylinder | `[skin.wound] enabled` |
| `mean_ca_wound` | mM | Mean Ca2+ in wound cylinder | `[skin.wound] enabled` |
| `mean_infl_wound` | a.u. | Mean inflammation in wound cylinder | `[skin.wound] enabled` |
| `scar_magnitude` | a.u. | Mean scar field in wound cylinder | `[skin.scar] enabled` |
| `mean_anti_infl_wound` | a.u. | Mean anti-inflammatory in wound | `[skin.inflammation] split_enabled` |
| `n_fibroblasts` | count | Active fibroblast agents (all states) | `[skin.fibroblast] enabled` |
| `n_myofibroblasts` | count | Myofibroblast-state fibroblasts | `[skin.fibroblast] enabled` |
| `mean_tgfb_wound` | a.u. | Mean TGF-beta in wound | `[skin.fibroblast] enabled` |
| `mean_collagen_wound` | a.u. | Mean collagen in wound | `[skin.fibroblast] enabled` |
| `mean_perfusion_wound` | normalized | Mean vascular perfusion in wound dermis | `[skin.wound] enabled` |
| `n_tumor_cells` | count | Active tumor cell agents | `[skin.tumor] enabled` |
| `n_tumor_cycling` | count | Tumor agents not in G0 (Ki-67 proxy) | `[skin.tumor] enabled` |
| `tumor_field_cells` | count | Tumor field voxels > 0.5 (handoff) | `[skin.tumor] enabled` |
| `mean_biofilm_wound` | a.u. | Mean biofilm in wound | `[skin.biofilm] enabled` |
| `mean_vegf_wound` | a.u. | Mean VEGF in wound | `[skin.angiogenesis] enabled` |
| `mean_mmp_wound` | a.u. | Mean MMP in wound | `[skin.mmp] enabled` |
| `mean_fibronectin_wound` | a.u. | Mean fibronectin in wound | `[skin.fibronectin] enabled` |
| `mean_elastin_wound` | a.u. | Mean elastin in wound | `[skin.elastin] enabled` |
| `mean_hyaluronan_wound` | a.u. | Mean hyaluronan in wound | `[skin.hyaluronan] enabled` |
| `mean_dermis_papillary` | normalized | Mean dermis integrity in papillary wound voxels | `[skin.dermis] enabled` |
| `mean_dermis_reticular` | normalized | Mean dermis integrity in reticular wound voxels | `[skin.dermis] enabled` |
| `mean_dermis_hypodermis` | normalized | Mean dermis integrity in hypodermis wound voxels | `[skin.dermis] enabled` |

`scripts/analysis/plot_metrics.py` generates figures from this CSV in `output/plots/`.
