> [Home](../../README.md) / [Modules](../README.md) / Scar

# Scar

Emergent scar formation driven by collagen deposition from the fibroblast cascade.

## Biology

Scar formation is a natural outcome of wound healing where collagen-rich repair tissue replaces the original tissue architecture. In this model, scar is not stamped unconditionally onto healed tissue but instead emerges from the local collagen concentration at the time a keratinocyte dissolves back into the continuum. Wound areas that heal before significant collagen accumulates (such as the wound margin) re-epithelialize normally, while the wound center, where myofibroblast activity was highest, develops scar. This makes scar location and density entirely dependent on the upstream fibroblast/collagen cascade.

## Model

Two scar mechanisms operate depending on configuration:

**Collagen-driven scar (with fibroblasts enabled):** When a cornified keratinocyte dissolves via WriteScarValue, it checks the local collagen concentration. If collagen exceeds `collagen_threshold`, the voxel becomes scar tissue (stratum value + 5); otherwise it becomes normal healthy stratum. Scar magnitude in the Scar PDE field mirrors collagen density directly.

**Inflammation-integral scar (fallback without fibroblasts):** A ScarAccumulationOp runs each step, accumulating scar as `scar += inflammation * accumulation_rate`. This models the clinical observation that prolonged inflammation leads to worse scarring (Ogawa 2017).

The proportional scar mode (`scar_proportional_enabled`) activates the ScarAccumulationOp for continuous scar field updates each step.

## Parameters

From modules/scar/config.toml:

| Parameter | Default | Units | Description | Source |
|-----------|---------|-------|-------------|--------|
| `enabled` | true | bool | Dissolved agents stamp scar into Stratum field | Convention |
| `proportional_enabled` | true | bool | Cumulative scar via collagen or inflammation | Ogawa 2017 ([DOI](https://doi.org/10.3390/ijms18030606)) |
| `accumulation_rate` | 0.001 | per step | Scar from inflammation (fallback) | Calibrated |
| `collagen_threshold` | 0.003 | a.u. | Collagen above this produces scar | Calibrated |

## Coupling

### Reads
| Field | Source module | How used |
|-------|-------------|----------|
| Collagen | fibroblast | Scar threshold check at dissolution |
| Inflammation | inflammation | Fallback accumulation when fibroblasts disabled |

### Writes
| Field | Consumer modules | What is written |
|-------|-----------------|-----------------|
| Scar | (visualization) | Scar magnitude per voxel |
| Stratum | tissue | Scar values (stratum + 5) at keratinocyte dissolution |

## Validation

| Dataset | Observable | Sources | Notes |
|---------|-----------|--------|-------|
| `scar_formation_parameters` | Scar severity and distribution | Ogawa 2017, Gauglitz 2011 | Scar correlates with inflammation integral |
| `scar_remodeling` | Collagen remodeling dynamics | Xue 2015, Bailey 1998 | Long-term scar maturation |

## Metrics

| Column | Units | Description |
|--------|-------|-------------|
| `scar_magnitude` | a.u. | Mean scar field in wound cylinder |

## Source files

| File | Purpose |
|------|---------|
| `scar_pde.h` | Scar accumulation PDE field |
| `scar_op.h` | ScarAccumulationOp for per-step scar updates |
| `config.toml` | Module configuration |
