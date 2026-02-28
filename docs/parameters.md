> [Home](../README.md) / [Docs](README.md) / Parameters

# Parameter Reference

Each module documents its own parameters with defaults, units, and literature sources in its README. This page provides a central index and explains the configuration layering system.

For configuration layering and profiles, see [Guide](guide.md). For architecture details, see [Architecture](architecture.md).

## Configuration layering

Parameters are layered at runtime:

1. `bdm.core.toml` contains core tissue parameters under `[skin]` and visualization/simulation settings
2. `modules/*/config.toml` contains per-module parameters under `[skin.modulename]`
3. `profiles/*.toml` provides skin phenotype overlays (normal, aged, diabetic)
4. `presets/*.toml` provides scenario overlays (wound, tumor, diabetic_wound)

The merge script (`scripts/config/merge_config.py`) combines layers 1 and 2 into `bdm.toml` (gitignored). Profiles and presets are applied as sparse overrides via `scripts/config/apply_preset.py`.

## Module parameter index

| Module | Config file | Parameters |
|--------|-------------|------------|
| [Tissue](../modules/tissue/README.md#parameters) | `bdm.core.toml` `[skin]` | Cell cycle, division, calcium, differentiation, KGF, O2, water, mechanics, shedding |
| [Wound](../modules/wound/README.md#parameters) | `modules/wound/config.toml` | Wound geometry, trigger time, inward bias, resolution |
| [Inflammation](../modules/inflammation/README.md#parameters) | `modules/inflammation/config.toml` | Diffusion/decay, Hill thresholds, wound source, split mode |
| [Perfusion](../modules/perfusion/README.md#parameters) | `modules/perfusion/config.toml` | Per-layer profiles, angiogenesis rate/delay |
| [Immune](../modules/immune/README.md#parameters) | `modules/immune/config.toml` | Neutrophil/macrophage timing, cytokine rates, efferocytosis, chemotaxis |
| [Fibroblast](../modules/fibroblast/README.md#parameters) | `modules/fibroblast/config.toml` | State machine thresholds, TGF-beta/collagen PDE, recruitment waves |
| [MMP](../modules/mmp/README.md#parameters) | `modules/mmp/config.toml` | Diffusion/decay, production rates, degradation rates |
| [Fibronectin](../modules/fibronectin/README.md#parameters) | `modules/fibronectin/config.toml` | Deposition, serum leakage, migration boost |
| [Scar](../modules/scar/README.md#parameters) | `modules/scar/config.toml` | Collagen threshold, proportional mode, accumulation rate |
| [Angiogenesis](../modules/angiogenesis/README.md#parameters) | `modules/angiogenesis/config.toml` | VEGF diffusion/decay/production, hypoxia threshold |
| [Dermis](../modules/dermis/README.md#parameters) | `modules/dermis/config.toml` | Sub-layer densities, collagen recovery, MMP degradation |
| [Diabetic](../modules/diabetic/README.md#parameters) | `modules/diabetic/config.toml` | 16 scaling factors across all systems |
| [Biofilm](../modules/biofilm/README.md#parameters) | `modules/biofilm/config.toml` | Growth, clearance, seeding, M1 block threshold |
| [Tumor](../modules/tumor/README.md#parameters) | `modules/tumor/config.toml` | BCC cycle factor, soft CI, apoptosis, handoff |
| [Elastin](../modules/elastin/README.md#parameters) | `modules/elastin/config.toml` | Sub-layer densities, production, MMP degradation |
| [Hyaluronan](../modules/hyaluronan/README.md#parameters) | `modules/hyaluronan/config.toml` | Sub-layer densities, water retention, migration scaffold |

## Simulation settings

These parameters live in `bdm.core.toml` under `[skin]` and control simulation infrastructure:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `duration_days` | 30 | Simulation duration in days |
| `tissue_min` | 0 | Min x/y for boundary removal |
| `tissue_max` | 30 | Max x/y for boundary removal |
| `metrics_interval_h` | 10 | Hours between CSV rows (0 = disabled) |
| `metrics_autoopen` | true | Open CSV in default app after simulation |
| `subcycle_slow` | 5 | PDE sub-cycling for slow fields (water, hyaluronan, vascular) |
| `subcycle_medium` | 3 | PDE sub-cycling for medium fields (inflammation, TGF-beta, MMP, VEGF) |
| `hot_reload` | false | Watch bdm.toml for runtime parameter changes |

## Debug flags

Under `[skin.debug]` in `bdm.core.toml`:

| Flag | Default | What it enables |
|------|---------|----------------|
| `immune` | false | Immune cell spawning, migration, state transitions |
| `fibroblast` | false | Collagen, TGF-beta, state transitions |
| `wound` | false | Wound creation, resolution, coverage |
| `scaled_grid` | false | ScaledGrid initialization diagnostics |
| `perf` | false | Wall-clock timing of major operations |
