> [Home](../../README.md) / [Modules](../README.md) / Perfusion

# Perfusion

Vascular perfusion field representing dermal microvasculature density and integrity.

## Biology
Dermal microvasculature provides oxygen and nutrients to the epidermis. Three vascular sub-layers have distinct perfusion characteristics: papillary dermis has dense capillary loops (highest perfusion), reticular dermis has arterioles and venules (moderate), and hypodermis has sparse large vessels (lowest). Wounding destroys vessels in the wound cylinder. Recovery occurs through angiogenesis (capillary sprouting from intact vessels), which is gated by tissue presence above the wound and modulated by VEGF when the angiogenesis module is enabled. Vascular perfusion directly controls O2 and water supply to the epidermis, making it a critical upstream field.

## Model
- Vascular PDE with diffusion (capillary sprouting) and no decay (stable vessels)
- Per-layer initialization profile based on Braverman 2000 (papillary=1.0, reticular=0.7, hypodermis=0.3)
- Wound response: zeros vascular field in wound dermis
- Source term: angiogenesis recovery gated by Stratum field (tissue presence required for growth factor signals)
- Per-layer rate factors for recovery (papillary fastest, hypodermis slowest)
- When angiogenesis module is enabled, recovery rate is modulated by local VEGF concentration
- Runs first in CompositeField so O2 and Water read updated perfusion each step

## Parameters
From modules/perfusion/config.toml:
| Parameter | Default | Units | Description | Source |
|-----------|---------|-------|-------------|--------|
| `diffusion` | 0.03 | - | Capillary sprouting / lateral spread | Calibrated |
| `decay` | 0.0 | - | No decay (stable vessels) | Convention |
| `basal` | 1.0 | normalized | Healthy dermis = fully perfused | Convention |
| `angio_rate` | 0.005 | per step | Angiogenesis recovery speed | Schugart et al. 2008 ([DOI](https://doi.org/10.1073/pnas.0711642105)) |
| `angio_delay_h` | 24 | hours | Hours before sprouting begins | Singer & Clark 1999 ([DOI](https://doi.org/10.1056/NEJM199909023411006)) |
| `papillary_fraction` | 1.0 | fraction | Dense capillary loops | Braverman 2000 ([DOI](https://doi.org/10.1046/j.1087-0024.2000.00010.x)) |
| `reticular_fraction` | 0.7 | fraction | Arterioles and venules | Braverman 2000 |
| `hypodermis_fraction` | 0.3 | fraction | Sparse large vessels | Braverman 2000 |
| `angio_papillary_factor` | 1.5 | multiplier | Capillaries regrow fastest | Calibrated |
| `angio_reticular_factor` | 1.0 | multiplier | Baseline rate | Calibrated |
| `angio_hypodermis_factor` | 0.5 | multiplier | Deep vessels recover slowly | Calibrated |

## Coupling
### Reads
| Field | Source module | How used |
|-------|-------------|----------|
| Stratum | tissue | Tissue presence gates angiogenesis recovery |
| VEGF | angiogenesis | Modulates recovery rate (when angiogenesis enabled) |

### Writes
| Field | Consumer modules | What is written |
|-------|-----------------|-----------------|
| Vascular | tissue (O2 pinning, Water pinning), fibronectin (serum leakage) | Perfusion density per voxel |

## Validation
| Dataset | Observable | Sources | Notes |
|---------|-----------|--------|-------|
| `perfusion_parameters` | Per-layer profiles, recovery kinetics | Braverman 2000, Singer & Clark 1999, Johnson & Wilgus 2014, Flegg et al. 2012 | Parameter derivation |

## Metrics
| Column | Units | Description |
|--------|-------|-------------|
| `mean_perfusion_wound` | normalized | Mean vascular perfusion in wound dermis |

## Source files
| File | Purpose |
|------|---------|
| `vascular.h` | Vascular PDE with per-layer profiles and angiogenesis recovery |
| `config.toml` | Module configuration |
