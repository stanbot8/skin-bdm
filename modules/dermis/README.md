> [Home](../../README.md) / [Modules](../README.md) / Dermis

# Dermis

Dermal tissue integrity field representing the structural state of the three dermal sub-layers.

## Biology

The dermis provides structural support and vascular supply to the epidermis. It consists of three functionally distinct sub-layers: the papillary dermis (loose connective tissue with dense capillary loops, directly beneath the epidermis), the reticular dermis (dense collagen and elastic fiber network providing mechanical strength), and the hypodermis (adipose-dominated tissue providing insulation and energy storage).

Wounding destroys dermal tissue within the wound cylinder. Recovery depends on new collagen deposition from myofibroblasts, making dermal repair mechanistically linked to the fibroblast cascade. MMP activity degrades dermal tissue, creating a remodeling feedback loop with collagen deposition.

## Model

**Dermis PDE:** minimal diffusion 0.001 (tissue does not flow), no decay.

**Initialization:** step-function sub-layer profile based on z-coordinate. Papillary dermis (z between 0 and dermal_z_papillary) initializes to `papillary_density` (1.0), reticular dermis to `reticular_density` (0.8), and hypodermis to `hypodermis_density` (0.5). Epidermis = 0.

**Wound response:** zeros all dermal wound voxels.

**Source term:** collagen-gated recovery with per-layer rate factors:
- If local collagen > `collagen_threshold`: `dermis += collagen_recovery_rate * layer_factor * (target - current)` per step
- MMP degradation: `dermis -= mmp_degradation * local_mmp` per step
- Layer factors: papillary 1.5x (recovers fastest), reticular 1.0x (baseline), hypodermis 0.5x (recovers slowest)

## Parameters

From modules/dermis/config.toml:

| Parameter | Default | Units | Description | Source |
|-----------|---------|-------|-------------|--------|
| `enabled` | true | bool | Master switch | Convention |
| `diffusion` | 0.001 | - | Minimal lateral coupling | Calibrated |
| `decay` | 0.0 | - | No spontaneous decay | Convention |
| `papillary_density` | 1.0 | normalized | Healthy papillary dermis | Clark 1996 ([DOI](https://doi.org/10.1007/978-1-4899-0185-9)) |
| `reticular_density` | 0.8 | normalized | Healthy reticular dermis | Clark 1996 |
| `hypodermis_density` | 0.5 | normalized | Healthy hypodermis | Clark 1996 |
| `collagen_threshold` | 0.05 | a.u. | Collagen needed for recovery | Calibrated |
| `collagen_recovery_rate` | 0.002 | per step | Recovery rate from collagen | Calibrated |
| `mmp_degradation` | 0.004 | per step | Tissue breakdown per MMP unit | Schultz & Wysocki 2009 ([DOI](https://doi.org/10.1111/j.1524-475X.2009.00466.x)) |
| `papillary_rate_factor` | 1.5 | multiplier | Papillary recovers fastest | Calibrated |
| `reticular_rate_factor` | 1.0 | multiplier | Reticular baseline rate | Calibrated |
| `hypodermis_rate_factor` | 0.5 | multiplier | Hypodermis recovers slowest | Calibrated |

## Coupling

### Reads
| Field | Source module | How used |
|-------|-------------|----------|
| Collagen | fibroblast | Recovery gate (collagen above threshold drives gain) |
| MMP | mmp | Degradation of dermal tissue |

### Writes
| Field | Consumer modules | What is written |
|-------|-----------------|-----------------|
| Dermis | inflammation (damage-proportional DAMPs) | Dermal tissue integrity per voxel |

## Validation

| Dataset | Observable | Sources | Notes |
|---------|-----------|--------|-------|
| `dermis_parameters` | Sub-layer profiles, recovery dynamics | Clark 1996, Schultz & Wysocki 2009 | Parameter derivation |

## Metrics

| Column | Units | Description |
|--------|-------|-------------|
| `mean_dermis_papillary` | normalized | Mean dermis integrity in papillary wound voxels |
| `mean_dermis_reticular` | normalized | Mean dermis integrity in reticular wound voxels |
| `mean_dermis_hypodermis` | normalized | Mean dermis integrity in hypodermis wound voxels |

## Source files

| File | Purpose |
|------|---------|
| `dermis_pde.h` | Dermis PDE with sub-layer profile and collagen-gated recovery |
| `config.toml` | Module configuration |
