> [Home](../../README.md) / [Modules](../README.md) / Angiogenesis

# Angiogenesis

VEGF-driven vascular recovery replacing flat-rate perfusion repair with hypoxia-dependent angiogenesis.

## Biology

Angiogenesis during wound healing is driven by the HIF-1alpha/VEGF signaling axis. When tissue oxygen drops below a threshold, hypoxia-inducible factor 1-alpha (HIF-1alpha) stabilizes and induces vascular endothelial growth factor (VEGF) production. VEGF stimulates endothelial cell sprouting from intact capillaries at the wound margin, gradually restoring perfusion to the wound bed. M2 macrophages also produce VEGF, coupling immune resolution to vascular repair.

Angiogenesis consumes VEGF proportionally, providing negative feedback that prevents excessive vessel growth. In diabetic wounds, HIF-1alpha is impaired (Thangarajah et al. 2009), reducing VEGF production and delaying vascular recovery.

## Model

**VEGF PDE:** diffusion 0.08, decay 0.01. Two sources:
- Hypoxia-driven: where O2 < `vegf_hypoxia_threshold`, VEGF is produced proportional to hypoxia severity (`vegf_production_rate * (threshold - local_O2)`)
- M2 macrophages: `m2_vegf_rate` per M2 macrophage per step

**Perfusion modulation:** replaces the flat `angio_rate` in the perfusion module with `angio_vegf_rate * local_vegf * demand`, where demand is the deficit from healthy perfusion.

**VEGF consumption:** angiogenesis consumes VEGF at `vegf_consumption_rate * recovery`, providing negative feedback.

**Feedback loop:**
```
Wound --> Destroyed vessels --> Low perfusion --> Low O2
                                                    |
                                HIF-1alpha --> VEGF production
                           (diabetic: impaired)     |
                                VEGF --> Endothelial sprouting --> Perfusion recovery
                                                    |
                                VEGF consumed (negative feedback)
                                                    |
                      M2 Macrophages --> VEGF (immune-vascular coupling)
```

## Parameters

From modules/angiogenesis/config.toml:

| Parameter | Default | Units | Description | Source |
|-----------|---------|-------|-------------|--------|
| `enabled` | true | bool | Master switch | Convention |
| `vegf_diffusion` | 0.08 | - | VEGF diffusion coefficient | Calibrated |
| `vegf_decay` | 0.01 | per step | VEGF natural decay | Calibrated |
| `vegf_production_rate` | 0.002 | per step | VEGF per hypoxic voxel | Johnson & Wilgus 2014 ([DOI](https://doi.org/10.1089/wound.2013.0517)) |
| `vegf_hypoxia_threshold` | 0.5 | normalized | O2 for HIF-1alpha induction | Schugart et al. 2008 ([DOI](https://doi.org/10.1073/pnas.0711642105)) |
| `vegf_consumption_rate` | 0.1 | fraction | VEGF consumed per unit recovery | Calibrated |
| `angio_vegf_rate` | 0.01 | per step | Perfusion recovery per unit VEGF | Flegg et al. 2012 ([DOI](https://doi.org/10.1016/j.jtbi.2012.01.043)) |
| `m2_vegf_rate` | 0.001 | per step | VEGF per M2 macrophage | Jetten et al. 2014 ([DOI](https://doi.org/10.1007/s10456-013-9381-6)) |

## Coupling

### Reads
| Field | Source module | How used |
|-------|-------------|----------|
| O2 | tissue | Hypoxia detection for VEGF production |

### Writes
| Field | Consumer modules | What is written |
|-------|-----------------|-----------------|
| VEGF | perfusion (modulates recovery rate) | VEGF concentration from hypoxia and M2 macrophages |

## Validation

| Dataset | Observable | Sources | Notes |
|---------|-----------|--------|-------|
| `angiogenesis_parameters` | VEGF dynamics, perfusion recovery coupling | Catrina 2004, Thangarajah 2009, Jetten 2014 | Parameter derivation |

## Literature data

Reference curves for validation (full citations in [SOURCES.yaml](SOURCES.yaml)):

| Dataset | File | Normalization |
|---------|------|---------------|
| VEGF kinetics | [vegf_kinetics.csv](data/vegf_kinetics.csv) | Peak = 1.0 |

## Metrics

| Column | Units | Description |
|--------|-------|-------------|
| `mean_vegf_wound` | a.u. | Mean VEGF in wound |

## Source files

| File | Purpose |
|------|---------|
| `vegf_pde.h` | VEGF diffusion field with hypoxia source |
| `vegf_op.h` | VEGF-modulated perfusion recovery operation |
| `config.toml` | Module configuration |
