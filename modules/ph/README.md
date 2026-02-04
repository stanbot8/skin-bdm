> [Home](../../README.md) / [Modules](../README.md) / pH

# pH

Wound surface pH field modeling the acid mantle disruption and recovery after injury.

## Biology

Healthy skin has an acidic surface pH of 5.5 maintained by the acid mantle (lactic acid, free fatty acids, filaggrin degradation products). Acute wounds expose the alkaline interstitial fluid (pH 7.4), creating a pH gradient across the wound bed. Re-acidification depends on keratinocyte re-epithelialization and takes several days. In chronic wounds, the alkaline environment persists and promotes bacterial colonization, impairs oxygen release (Bohr effect), and enhances MMP activity, creating feedback loops that delay healing.

## Model

Continuum PDE field. Semantics: 0 = normal acidic skin (pH 5.5), 1.0 = fresh wound alkaline (pH 7.4). Recovery driven by stratum restoration at configurable rate. Suppresses keratinocyte migration at high values.

## Parameters

From `modules/ph/config.toml`:

| Parameter | Default | Units | Description | Source |
|-----------|---------|-------|-------------|--------|
| `recovery_rate` | 0.003 | per step | Acidification rate, scaled by perfusion | Schneider et al. 2007 ([DOI](https://doi.org/10.1111/j.1524-475X.2007.00230.x)) |
| `migration_suppression` | 0.5 | fraction | Max migration speed reduction at full alkalinity | Gethin 2007 ([DOI](https://doi.org/10.1111/j.1742-481X.2007.00340.x)) |
| `mmp_boost` | 0.5 | fraction | MMP activity amplification at full alkalinity | Leveen et al. 1994 ([DOI](https://doi.org/10.1016/S0140-6736(94)91421-4)) |
| `biofilm_boost` | 0.4 | fraction | Biofilm growth rate boost at full alkalinity | Gethin 2007 |
| `bohr_factor` | 0.3 | fraction | O2 delivery reduction at full alkalinity (Bohr effect) | Leveen et al. 1994 |

## Coupling

### Reads

| Field | Source module | How used |
|-------|-------------|----------|
| Stratum | tissue | Re-epithelialization drives pH recovery |

### Writes

| Field | Consumer modules | What is written |
|-------|-----------------|-----------------|
| pH | tissue (migration suppression) | Wound bed alkalinity level |

## Metrics

| Column | Units | Description |
|--------|-------|-------------|
| `mean_ph_wound` | a.u. | Mean pH in wound |

## Source files

| File | Purpose |
|------|---------|
| `ph_pde.h` | pH continuum field |
| `config.toml` | Module configuration |
