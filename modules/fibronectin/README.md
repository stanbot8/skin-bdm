> [Home](../../README.md) / [Modules](../README.md) / Fibronectin

# Fibronectin

Provisional wound matrix scaffold that facilitates keratinocyte migration during re-epithelialization.

## Biology

Fibronectin forms the provisional matrix that keratinocytes use as a migration scaffold during wound healing. Two sources contribute to the wound fibronectin matrix: plasma fibronectin leaks from damaged blood vessels into the wound bed (early, passive), and cellular fibronectin is deposited by activated fibroblasts (later, active). Keratinocytes bind to fibronectin through alpha5-beta1 integrins, which increases their migration speed across the wound surface (Grinnell 1984).

MMPs degrade fibronectin as part of the ECM remodeling process, eventually replacing the provisional matrix with mature collagen-based tissue.

## Model

Fibronectin PDE with no diffusion (structural ECM) and slow decay 0.0005 (natural turnover).

Sources:
- Activated fibroblasts deposit cellular fibronectin at `deposition_rate` per step
- Vascular leakage deposits plasma fibronectin at `serum_rate` per step (proportional to local perfusion)

Effect on keratinocytes: migration speed multiplied by `(1 + migration_boost * local_fibronectin)`, modeling integrin-FN binding.

MMP degradation applied in FusedWoundSourceOp at `mmp_fibronectin_degradation * local_mmp` per step.

## Parameters

From modules/fibronectin/config.toml:

| Parameter | Default | Units | Description | Source |
|-----------|---------|-------|-------------|--------|
| `enabled` | true | bool | Master switch | Convention |
| `decay` | 0.0005 | per step | Slow turnover (~1400 steps half-life) | Calibrated |
| `deposition_rate` | 0.003 | per step | Cellular fibronectin from activated fibroblasts | Clark 1990 ([DOI](https://doi.org/10.1111/1523-1747.ep12876104)) |
| `serum_rate` | 0.001 | per step | Plasma fibronectin from vascular leakage | Calibrated |
| `migration_boost` | 0.5 | multiplier | Keratinocyte speed from integrin-FN binding | Grinnell 1984 ([DOI](https://doi.org/10.1002/jcb.240260206)) |

## Coupling

### Reads
| Field | Source module | How used |
|-------|-------------|----------|
| Vascular | perfusion | Serum leakage proportional to perfusion |
| MMP | mmp | Fibronectin degradation |

### Writes
| Field | Consumer modules | What is written |
|-------|-----------------|-----------------|
| Fibronectin | tissue (migration speed boost) | Provisional matrix concentration |

## Validation

| Dataset | Observable | Sources | Notes |
|---------|-----------|--------|-------|
| `fibronectin_provisional_matrix` | Fibronectin deposition and turnover | Clark 1990, Grinnell 1984, Hynes 1990 | Parameter derivation |

## Literature data

Reference curves for validation (full citations in [SOURCES.yaml](SOURCES.yaml)):

| Dataset | File | Normalization |
|---------|------|---------------|
| Fibronectin kinetics | [fibronectin_kinetics.csv](data/fibronectin_kinetics.csv) | Peak = 1.0 |

## Metrics

| Column | Units | Description |
|--------|-------|-------------|
| `mean_fibronectin_wound` | a.u. | Mean fibronectin in wound |

## Source files

| File | Purpose |
|------|---------|
| `fibronectin_pde.h` | Fibronectin structural deposit field |
| `config.toml` | Module configuration |
