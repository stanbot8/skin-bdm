> [Home](../../README.md) / [Modules](../README.md) / Elastin

# Elastin

Elastic fiber network in the dermis representing tropoelastin production and MMP-mediated degradation.

## Biology

Elastic fibers provide the dermis with recoil and resilience, allowing skin to return to its original shape after deformation. The elastic fiber network consists of an elastin core surrounded by fibrillin microfibrils. Elastin has an extraordinarily long half-life of approximately 70 years, making it one of the most stable proteins in the human body. New elastic fiber production by fibroblasts (via tropoelastin secretion) is extremely slow, which is why wounds and scars permanently lose their elasticity.

Elastin is degraded by elastase and MMP-9 (produced by M1 macrophages). In the context of wound healing, elastic fibers within the wound cylinder are disrupted and recover very slowly compared to collagen, contributing to the mechanical differences between scar tissue and normal skin.

This module is a groundwork module (disabled by default) that adds the elastin PDE field without affecting wound healing dynamics.

## Model

**Elastin PDE:** no diffusion (structural ECM, immobile), very slow decay 0.0002 (natural turnover, half-life ~70 years).

**Initialization:** sub-layer profile based on z-coordinate. Reticular dermis has the highest density (`basal_density` = 0.5), papillary dermis has thinner fibers (`papillary_density` = 0.3). Epidermis = 0.

**Wound response:** zeros elastin in wound voxels (fiber disruption).

**Source term:** fibroblast production at `production_rate` per step (when fibroblast module active). MMP degradation at `mmp_degradation * local_mmp` per step.

## Parameters

From modules/elastin/config.toml:

| Parameter | Default | Units | Description | Source |
|-----------|---------|-------|-------------|--------|
| `enabled` | false | bool | Master switch (groundwork module) | Convention |
| `diffusion` | 0.0 | - | Structural ECM, immobile | Convention |
| `decay` | 0.0002 | per step | Very slow turnover (t1/2 ~70 years) | Kielty et al. 2002 ([DOI](https://doi.org/10.1002/path.1107)) |
| `basal_density` | 0.5 | normalized | Baseline in reticular dermis | Kielty et al. 2002 |
| `papillary_density` | 0.3 | normalized | Thinner fibers in papillary dermis | Kielty et al. 2002 |
| `production_rate` | 0.0005 | per step | Tropoelastin from fibroblasts | Mithieux & Weiss 2005 ([DOI](https://doi.org/10.1016/S0065-3233(05)70013-9)) |
| `mmp_degradation` | 0.003 | per step | Elastase/MMP-9 degrades elastin | Kielty et al. 2002 |

## Coupling

### Reads
| Field | Source module | How used |
|-------|-------------|----------|
| MMP | mmp | Degradation of elastic fibers |

### Writes
| Field | Consumer modules | What is written |
|-------|-----------------|-----------------|
| Elastin | (visualization, future mechanical coupling) | Elastic fiber density per voxel |

## Validation

No validation dataset yet. Parameter values derived from:
- Kielty et al. 2002: elastic fiber structure and turnover
- Mithieux & Weiss 2005: tropoelastin production rates
- Almine et al. 2012: elastin biology review

## Metrics

| Column | Units | Description |
|--------|-------|-------------|
| `mean_elastin_wound` | a.u. | Mean elastin in wound |

## Source files

| File | Purpose |
|------|---------|
| `elastin_pde.h` | Elastin PDE with sub-layer profile and MMP degradation |
| `config.toml` | Module configuration |
