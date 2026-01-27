> [Home](../../README.md) / [Modules](../README.md) / MMP

# MMP

Matrix metalloproteinase dynamics mediating extracellular matrix remodeling during wound healing.

## Biology

Matrix metalloproteinases (MMPs) are zinc-dependent endopeptidases that degrade extracellular matrix components. During wound healing, M1 macrophages produce MMP-9 (gelatinase B) for debris clearance, while fibroblasts produce MMP-1 (interstitial collagenase) and MMP-3 (stromelysin) for ECM remodeling. MMPs degrade collagen, fibronectin, dermis, and elastin. Tissue inhibitors of metalloproteinases (TIMPs) provide natural negative regulation.

In diabetic wounds, MMP production is elevated 3-fold while TIMP-mediated decay is halved, creating an MMP/TIMP imbalance that prevents collagen accumulation and impairs wound healing (Lobmann et al. 2002).

## Model

MMP PDE with diffusion 0.02 and TIMP-mediated decay 0.02.

Sources: M1 macrophages produce MMP-9 at `m1_rate` per step; fibroblasts produce MMP-1/3 at `fibroblast_rate` per step.

Degradation targets (applied in FusedWoundSourceOp):
- Collagen: `collagen -= collagen_degradation * local_mmp` per step
- Fibronectin: `fibronectin -= fibronectin_degradation * local_mmp` per step
- Dermis: `dermis -= dermis_mmp_degradation * local_mmp` per step (from dermis module)
- Elastin: `elastin -= elastin_mmp_degradation * local_mmp` per step (from elastin module)

## Parameters

From modules/mmp/config.toml:

| Parameter | Default | Units | Description | Source |
|-----------|---------|-------|-------------|--------|
| `enabled` | true | bool | Master switch | Convention |
| `diffusion` | 0.02 | - | MMP diffusion coefficient | Nagase et al. 1999 ([DOI](https://doi.org/10.1074/jbc.274.31.21491)) |
| `decay` | 0.02 | per step | TIMP-mediated inactivation | Calibrated |
| `m1_rate` | 0.003 | per step | MMP-9 per M1 macrophage | Lobmann et al. 2002 ([DOI](https://doi.org/10.1007/s00125-002-0868-8)) |
| `fibroblast_rate` | 0.001 | per step | MMP-1/3 per fibroblast | Nagase et al. 1999 |
| `collagen_degradation` | 0.002 | per step | MMP-mediated collagen lysis | Ladwig et al. 2002 ([DOI](https://doi.org/10.1046/j.1524-475x.2002.10903.x)) |
| `fibronectin_degradation` | 0.002 | per step | MMP-mediated fibronectin lysis | Calibrated |

## Coupling

### Reads
| Field | Source module | How used |
|-------|-------------|----------|
| (none) | - | MMP is written by agents, not read from other fields |

### Writes
| Field | Consumer modules | What is written |
|-------|-----------------|-----------------|
| MMP | collagen (degradation), fibronectin (degradation), dermis (degradation), elastin (degradation) | MMP concentration from immune and fibroblast agents |

## Validation

| Dataset | Observable | Sources | Notes |
|---------|-----------|--------|-------|
| `mmp_ecm_remodeling` | MMP levels, collagen/fibronectin degradation | Lobmann 2002, Nagase 1999, Ladwig 2002 | Parameter derivation and diabetic imbalance |

## Literature data

Reference curves for validation (full citations in [SOURCES.yaml](SOURCES.yaml)):

| Dataset | File | Normalization |
|---------|------|---------------|
| MMP kinetics | [mmp_kinetics.csv](data/mmp_kinetics.csv) | Peak = 1.0 |

## Metrics

| Column | Units | Description |
|--------|-------|-------------|
| `mean_mmp_wound` | a.u. | Mean MMP in wound |

## Source files

| File | Purpose |
|------|---------|
| `mmp_pde.h` | MMP diffusion field with TIMP decay |
| `config.toml` | Module configuration |
