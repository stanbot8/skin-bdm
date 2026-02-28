> [Home](../README.md) / [Docs](README.md) / Guide

# Skibidy Guide

Usage guide for the skin tissue simulation. See also: [Architecture](architecture.md) | [Parameters](parameters.md) | [Treatments](treatments.md)

## Configuration

Skibidy uses a layered TOML configuration system. Core tissue parameters live in `bdm.core.toml` (committed to git). Module-specific parameters live in `modules/*/config.toml`. At runtime, `scripts/config/merge_config.py` merges them into `bdm.toml` (gitignored), which BioDynaMo reads.

| Layer | File | Purpose |
|-------|------|---------|
| Core | `bdm.core.toml` | Tissue biology, visualization, simulation bounds |
| Modules | `modules/*/config.toml` | Per-module params in TOML subsections |
| Profiles | `profiles/*.toml` | Biological phenotype overlays (skin type) |
| Presets | `presets/*.toml` | Scenario overlays (what happens to the skin) |

**Config layering order:** `bdm.core.toml` + `modules/*/config.toml` (merged) -> skin profile (`--skin=NAME`) -> preset (`--preset=NAME`)

The merge script auto-discovers module files sorted alphabetically and appends enabled ones to the core config. Modules with `enabled = false` are skipped. Profiles and presets are then applied as sparse overrides using `scripts/config/apply_preset.py`.

### Module files

Each module directory in `modules/` contains a `config.toml` owning one TOML subsection under `[skin]`, plus all related source files (agents, behaviors, events, PDE headers):

| Module | Section | Description |
|--------|---------|-------------|
| `wound/` | `[skin.wound]` | Punch biopsy event |
| `perfusion/` | `[skin.perfusion]` | Vascular perfusion field |
| `inflammation/` | `[skin.inflammation]` | Cytokine diffusion + split inflammation |
| `immune/` | `[skin.immune]` | Neutrophils, macrophages, efferocytosis, chemotaxis |
| `diabetic/` | `[skin.diabetic]` | Diabetic impairment modifiers |
| `fibroblast/` | `[skin.fibroblast]` | Fibroblast lifecycle, TGF-beta, collagen |
| `mmp/` | `[skin.mmp]` | Matrix metalloproteinase dynamics |
| `fibronectin/` | `[skin.fibronectin]` | Provisional matrix scaffold |
| `angiogenesis/` | `[skin.angiogenesis]` | VEGF-driven angiogenesis |
| `biofilm/` | `[skin.biofilm]` | Bacterial biofilm colonization |
| `tumor/` | `[skin.tumor]` | BCC neoplastic growth |
| `scar/` | `[skin.scar]` | Scar field and proportional scarring |
| `dermis/` | `[skin.dermis]` | Dermal tissue integrity (sub-layer profile) |
| `elastin/` | `[skin.elastin]` | Elastic fiber network |
| `hyaluronan/` | `[skin.hyaluronan]` | Hyaluronic acid ground substance |
| `tissue/` | `[skin]` (root) | Core epidermal biology (always on) |

Key names within modules drop their redundant prefix since the section provides context. For example, `wound_enabled` becomes `enabled` under `[skin.wound]`, and `diabetic_m1_duration_factor` becomes `m1_duration_factor` under `[skin.diabetic]`.

For the full parameter reference, see [Parameters](parameters.md).

### Creating a new module

1. Create `modules/mymodule/config.toml` with a `[skin.mymodule]` section
2. Add source files (PDE headers, agents, behaviors) in `modules/mymodule/`
3. Add corresponding C++ members to `sim_param.h`
4. Bind them in `skibidy.cc` with `BDM_ASSIGN_CONFIG_VALUE(member, "skin.mymodule.key")`
5. The merge script picks it up automatically on next run

For the full module development guide (PDE subclass, agents, behaviors, metrics, tests), see [Architecture: How to add a module](architecture.md#how-to-add-a-module).

### Skin profiles

Parameters fall into three categories:

| Category | Examples | Where |
|----------|----------|-------|
| **Skin biology** | Cell cycle, calcium gradient, perfusion baseline, layer geometry | `profiles/*.toml` |
| **Scenario** | Wound timing, tumor seeding, module enable flags, num_steps | `presets/*.toml` |
| **Simulation plumbing** | Metrics interval, visualization export, tissue bounds | `bdm.core.toml` only |

Skin profiles (`profiles/`) define tissue biology, i.e. what kind of skin. Presets (`presets/`) define scenarios, i.e. what happens to that skin. Both are sparse TOML overlays applied to the merged config using `apply_preset.py`. Overlays can target any section (e.g. `[skin]`, `[skin.wound]`, `[skin.diabetic]`).

```bash
./run.sh                              # interactive menu (pick skin + preset)
./run.sh --skin=aged --preset=wound   # aged skin + wound scenario
./run.sh --list-skins                 # show available profiles
```

**Built-in profiles:**

| Profile | Description |
|---------|-------------|
| `normal` | Healthy adult skin (matches defaults, no-op) |
| `aged` | Slower cell cycle, fewer stem cells, thinner corneum, reduced perfusion |
| `diabetic` | Prolonged M1 phase, slow resolution, microangiopathy |

**Creating a custom profile:**

1. Copy `profiles/TEMPLATE.toml` to `profiles/myprofile.toml`
2. Uncomment and change only the keys you want to override
3. Run with `./run.sh --skin=myprofile --preset=wound`

The template contains all biology-relevant keys organized by system (cell cycle, calcium, perfusion, immune kinetics, fibroblast/ECM, diabetic modifiers). Only include keys that differ from the defaults. Profiles are sparse TOML overlays: any key you omit keeps its default value. Missing keys are automatically appended to the correct section.

Custom profiles appear in the interactive menu and in `--list-skins` output.

Example (neonatal skin with faster turnover and thinner epidermis):

```toml
# Profile: neonatal -- fast-cycling thin skin
[skin]
g1_duration = 5.0
s_duration = 4.0
stem_fraction = 0.60
volume_z_cornified = 18.0

[skin.perfusion]
basal = 1.2
```

## Running

```bash
./run.sh                              # interactive menu
./run.sh --preset=wound               # wound healing scenario
./run.sh --preset=tumor               # tumor growth only
./run.sh --preset=tumor_wound         # tumor at t=0, wound at day 8
./run.sh --preset=nothing             # homeostatic epidermis, no events
./run.sh --diabetic                   # shorthand for --preset=diabetic_wound
./run.sh --compare                    # run normal + diabetic back-to-back
./run.sh --no-view                    # skip ParaView viewer
./run.sh --no-validate                # skip validation (literature comparison)
./run.sh --list-presets               # show available presets
./run.sh --list-skins                 # show available skin profiles
./tests/test.sh                       # build and run unit tests
```

`run.sh` handles the full pipeline: merge config from `bdm.core.toml` + modules, apply profile and preset overlays, build (if needed), run the simulation, plot metrics, run validation, and launch the viewer.

The `--compare` mode runs normal and diabetic wound simulations back-to-back and generates a comparison.

## Visualization

The patch script (`scripts/viz/patch_pvsm.py`) reads colormap, glyph, rendering, and camera settings from `bdm.core.toml` and applies them to the ParaView state file:

| Stratum | Value | Color |
|---------|-------|-------|
| Basale | 0 | pinkish-red |
| Spinosum | 1 | salmon pink |
| Granulosum | 2 | pale amber |
| Corneum | 3 | parchment beige |
| Scar Basale | 5 | muted pink |
| Scar Spinosum | 6 | muted salmon |
| Scar Granulosum | 7 | muted amber |
| Scar Corneum | 8 | muted parchment |
| Papillary dermis | 10 | light peach |
| Reticular dermis | 11 | warm brown |
| Hypodermis | 12 | deep brown |

Scar tissue uses stratum+5 encoding (emergent: based on local collagen at dissolution time). Dermal sub-layers start at 10 to avoid collision with scar values.

### Animated GIF

After `./run.sh`, use ParaView **File > Save Animation** to export frames, then:
```bash
./scripts/viz/render_gif.sh        # default 6 fps
./scripts/viz/render_gif.sh 10     # faster playback
```

## Metrics and validation

The simulation writes `output/skibidy/metrics.csv` every `metrics_interval` steps with cell populations, stratum counts, wound closure %, and field means. `scripts/analysis/plot_metrics.py` generates figures in `output/plots/`.

For the full metrics column reference, see [Docs: Metrics columns](README.md#metrics-columns). For the validation framework and literature comparison scripts, see the [validation suite](../literature/README.md).

```bash
python3 literature/validate_all.py          # full dashboard with RMSE
python3 literature/compare_wound.py         # 2x2 wound healing summary
python3 literature/compare_fibroblast.py    # myofibroblast + collagen
python3 literature/compare_tumor.py         # tumor growth + Ki-67
```

## Research modules

| Module | Status | Description |
|--------|--------|-------------|
| Wound healing | Implemented | Punch biopsy, immune response, re-epithelialization, scar formation |
| Fibroblast/scar | Implemented | TGF-beta-driven myofibroblast differentiation, collagen deposition |
| MMP remodeling | Implemented | Matrix metalloproteinase dynamics, collagen/fibronectin degradation |
| Fibronectin | Implemented | Provisional matrix scaffold deposition by fibroblasts |
| VEGF/Angiogenesis | Implemented | Hypoxia-driven VEGF, endothelial sprouting |
| Tumor (BCC) | Implemented | Self-contained neoplastic growth, composable with wound events |
| Vascular perfusion | Implemented | Vessel density field, wound disruption, angiogenesis recovery |
| Biofilm | Implemented | Bacterial colonization, immune clearance, PAMP-driven inflammation |
| Skin cancer (advanced) | Planned | Immune evasion, invasion/metastasis |
| Pathogen response | Planned | Viral/bacterial agents (HSV, HPV, HIV), immune activation |
| Aging | Planned | Stem cell exhaustion, senescence, epidermal thinning |
| Hair follicle | Planned | Bulge stem cells, dermal papilla signaling |
| Skincare | Planned | Topical substance penetration, retinoid effects |
