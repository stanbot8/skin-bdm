> [Home](../../README.md) / [Modules](../README.md) / Hyaluronan

# Hyaluronan

Hyaluronic acid ground substance in the dermis providing water retention and migration scaffolding.

## Biology

Hyaluronic acid (HA) is a high-molecular-weight glycosaminoglycan that is a major component of the dermal extracellular matrix. HA has exceptional water-binding capacity, retaining up to 1000 times its weight in water, which maintains tissue hydration and turgor. In the papillary dermis, HA concentration is highest, creating a hydrated environment that facilitates cell migration and nutrient diffusion.

During wound healing, HA in the wound area is disrupted. Fibroblasts produce new HA as part of the provisional matrix, and HA facilitates keratinocyte migration by providing a hydrated scaffold. HA is degraded by hyaluronidases with moderate turnover.

This module is a groundwork module (disabled by default) that adds the hyaluronan PDE field with water retention and migration scaffold coupling.

## Model

**Hyaluronan PDE:** slow diffusion 0.01 (high molecular weight HA is viscous), hyaluronidase decay 0.005.

**Initialization:** sub-layer profile. Papillary dermis has highest concentration (`basal_density` = 0.8, HA-rich). Reticular dermis has lower concentration (`reticular_density` = 0.4). Epidermis = 0.

**Wound response:** zeros hyaluronan in wound voxels.

**Source term:** fibroblast production at `production_rate` per step (when fibroblast module active).

**Water retention coupling:** local water capacity is augmented by `water_retention_factor * local_HA`, increasing tissue hydration where HA is present.

**Migration scaffold coupling:** keratinocyte migration speed is augmented by `migration_scaffold_factor * local_HA`, facilitating cell movement through HA-rich matrix.

## Parameters

From modules/hyaluronan/config.toml:

| Parameter | Default | Units | Description | Source |
|-----------|---------|-------|-------------|--------|
| `enabled` | false | bool | Master switch (groundwork module) | Convention |
| `diffusion` | 0.01 | - | Slow lateral spread (high MW HA) | Calibrated |
| `decay` | 0.005 | per step | Hyaluronidase turnover | Stern et al. 2006 ([DOI](https://doi.org/10.1016/j.ejcb.2006.05.009)) |
| `basal_density` | 0.8 | normalized | Baseline in papillary dermis (HA-rich) | Toole 2004 ([DOI](https://doi.org/10.1038/nrc1391)) |
| `reticular_density` | 0.4 | normalized | Lower in reticular dermis | Toole 2004 |
| `production_rate` | 0.002 | per step | Fibroblast HA synthesis | Calibrated |
| `water_retention_factor` | 0.5 | multiplier | HA contribution to local water capacity | Calibrated |
| `migration_scaffold_factor` | 0.3 | multiplier | HA contribution to migration speed | Calibrated |

## Coupling

### Reads
| Field | Source module | How used |
|-------|-------------|----------|
| (none) | - | HA is produced by fibroblast agents, not read from fields |

### Writes
| Field | Consumer modules | What is written |
|-------|-----------------|-----------------|
| Hyaluronan | tissue (water retention, migration scaffold) | HA concentration per voxel |

## Validation

No validation dataset yet. Parameter values derived from:
- Toole 2004: HA biology and distribution
- Stern et al. 2006: hyaluronidase turnover
- Aya & Bhangal 2014: HA in wound healing
- Chen & Abatangelo 1999: HA functions in wound repair

## Metrics

| Column | Units | Description |
|--------|-------|-------------|
| `mean_hyaluronan_wound` | a.u. | Mean hyaluronan in wound |

## Source files

| File | Purpose |
|------|---------|
| `hyaluronan_pde.h` | Hyaluronan PDE with sub-layer profile and fibroblast production |
| `config.toml` | Module configuration |
