> [Home](../../README.md) / [Modules](../README.md) / Diabetic

# Diabetic

Chronic wound modifiers that impair healing across immune, keratinocyte, fibroblast, vascular, and ECM remodeling systems.

## Biology

Diabetic wound healing is impaired through compounding dysfunctions across multiple systems. Advanced glycation end-products (AGEs) bind RAGE receptors on immune cells, prolonging the M1 pro-inflammatory macrophage phenotype and delaying the M1 to M2 transition critical for wound resolution. Excess neutrophils with extended lifespans sustain tissue damage. Keratinocyte proliferation and migration are directly impaired. Fibroblasts activate slowly, deposit less collagen, and have shortened lifespans. MMP production is elevated while TIMP-mediated inhibition is reduced, creating an ECM degradation imbalance. Vascular HIF-1alpha signaling is impaired, reducing VEGF production and delaying angiogenesis.

These dysfunctions compound: prolonged inflammation damages tissue, impaired fibroblasts fail to rebuild ECM, excess MMPs degrade what little collagen is deposited, and poor vascularization limits oxygen delivery. The result is a chronic non-healing wound phenotype.

This module has no PDE field of its own. It operates by scaling parameters of other modules at runtime, plus a BaselineInflammationOp that adds constant low-grade AGE/RAGE-driven inflammation to wound voxels.

## Model

**Multiplicative scaling factors** applied to other modules:

| System | Factor | Effect |
|--------|--------|--------|
| Immune: M1 duration | `m1_duration_factor` (3.0x) | Prolonged pro-inflammatory phase |
| Immune: M2 resolution | `resolution_factor` (0.3x) | Impaired anti-inflammatory clearance |
| Immune: efferocytosis | `efferocytosis_factor` (0.5x) | Reduced neutrophil clearance |
| Immune: neutrophil count | `neutrophil_factor` (1.8x) | Excess neutrophil infiltration |
| Immune: neutrophil lifespan | `neutrophil_lifespan_factor` (1.5x) | Delayed neutrophil apoptosis |
| Tissue: proliferation | `prolif_factor` (0.5x) | Impaired keratinocyte division |
| Tissue: migration | `migration_factor` (0.6x) | Slower re-epithelialization |
| Fibroblast: activation | `fibroblast_activation_factor` (2.0x) | Delayed fibroblast activation |
| Fibroblast: collagen | `collagen_factor` (0.4x) | Reduced collagen deposition |
| Fibroblast: lifespan | `fibroblast_lifespan_factor` (0.6x) | Shortened fibroblast survival |
| MMP: production | `mmp_factor` (3.0x) | Elevated MMP levels |
| MMP: TIMP decay | `timp_factor` (0.5x) | Reduced MMP inactivation |
| Angiogenesis: VEGF | `vegf_factor` (0.4x) | Impaired HIF-1alpha response |
| Biofilm: clearance | `biofilm_clearance_factor` (0.5x) | Impaired neutrophil phagocytosis |

**BaselineInflammationOp:** adds `baseline_inflammation` per step to wound voxels, modeling chronic low-grade AGE/RAGE-driven NF-kB activation independent of immune cell activity.

**Inflammation sensitivity:** keratinocyte Hill-function thresholds are divided by `inflammation_sensitivity` (2.0x), making keratinocytes more responsive to lower inflammation levels.

## Parameters

From modules/diabetic/config.toml:

| Parameter | Default | Units | Description | Source |
|-----------|---------|-------|-------------|--------|
| `mode` | false | bool | Master switch for diabetic impairments | Convention |
| `m1_duration_factor` | 3.0 | multiplier | M1 duration scaling | Mirza & Koh 2011 ([DOI](https://doi.org/10.1016/j.cyto.2011.06.016)) |
| `resolution_factor` | 0.3 | multiplier | M2 resolution rate scaling | Louiselle et al. 2021 ([DOI](https://doi.org/10.1016/j.trsl.2021.05.006)) |
| `efferocytosis_factor` | 0.5 | multiplier | Efferocytosis probability scaling | Khanna et al. 2010 ([DOI](https://doi.org/10.1371/journal.pone.0009539)) |
| `neutrophil_factor` | 1.8 | multiplier | Neutrophil spawn count scaling | Wong et al. 2015 ([DOI](https://doi.org/10.1038/nm.3887)) |
| `neutrophil_lifespan_factor` | 1.5 | multiplier | Neutrophil lifespan scaling | Brem & Tomic-Canic 2007 ([DOI](https://doi.org/10.1172/JCI32169)) |
| `prolif_factor` | 0.5 | multiplier | Keratinocyte G1 to S rate scaling | Wang & Graves 2020 ([DOI](https://doi.org/10.1155/2020/3714704)) |
| `migration_factor` | 0.6 | multiplier | Keratinocyte migration speed scaling | Wang & Graves 2020 |
| `fibroblast_activation_factor` | 2.0 | multiplier | Activation delay scaling (>1 = slower) | Lerman et al. 2003 |
| `collagen_factor` | 0.4 | multiplier | Collagen deposition scaling | Goldberg et al. 2007 ([DOI](https://doi.org/10.1038/sj.jid.5700890)) |
| `fibroblast_lifespan_factor` | 0.6 | multiplier | Fibroblast lifespan scaling | Lerman et al. 2003 |
| `baseline_inflammation` | 0.001 | per step | AGE/RAGE chronic inflammation | Bierhaus et al. 2005 ([DOI](https://doi.org/10.1007/s00109-005-0688-7)) |
| `inflammation_sensitivity` | 2.0 | multiplier | Keratinocyte inflammation sensitivity | Rasik & Shukla 2000 ([DOI](https://doi.org/10.1046/j.1365-2613.2000.00158.x)) |
| `mmp_factor` | 3.0 | multiplier | MMP production scaling | Lobmann et al. 2002 ([DOI](https://doi.org/10.1007/s00125-002-0868-8)) |
| `timp_factor` | 0.5 | multiplier | TIMP decay scaling | Lobmann et al. 2002 |
| `biofilm_clearance_factor` | 0.5 | multiplier | Immune clearance scaling | Bjarnsholt et al. 2008 ([DOI](https://doi.org/10.1111/j.1524-475X.2007.00283.x)) |
| `vegf_factor` | 0.4 | multiplier | HIF-1alpha impairment | Thangarajah et al. 2009 ([DOI](https://doi.org/10.1073/pnas.0906670106)) |

## Coupling

### Modifies

| Target module | Parameters affected | Effect |
|--------------|-------------------|--------|
| immune | M1 duration, resolution rate, efferocytosis, neutrophil count/lifespan | Prolonged inflammation, excess neutrophils |
| tissue | Proliferation rate, migration speed, inflammation sensitivity | Impaired re-epithelialization |
| fibroblast | Activation delay, collagen rate, lifespan | Impaired ECM repair |
| mmp | MMP production, TIMP decay | ECM degradation imbalance |
| angiogenesis | VEGF production factor | Impaired vascular recovery |
| biofilm | Immune clearance rates | Increased infection susceptibility |
| inflammation | Baseline source rate | Chronic low-grade inflammation |

### Writes

| Field | Consumer modules | What is written |
|-------|-----------------|-----------------|
| Inflammation | (via BaselineInflammationOp) | Constant AGE/RAGE-driven inflammation in wound voxels |

## Validation

| Dataset | Observable | Sources | Notes |
|---------|-----------|--------|-------|
| `diabetic_wound_healing` | Parameter derivation for all factors | 28 sources | Comprehensive literature review |
| `diabetic_closure_kinetics` | Delayed wound closure | Multiple | ~42 days vs ~22 days normal |
| `diabetic_immune_cell_kinetics` | Prolonged neutrophil/macrophage counts | Multiple | Elevated and sustained |
| `diabetic_inflammation_timecourse` | Chronic inflammation plateau | Multiple | Failed resolution |

## Literature data

Reference curves for validation (full citations in [SOURCES.yaml](SOURCES.yaml)):

| Dataset | File | Normalization |
|---------|------|---------------|
| Diabetic closure kinetics | [diabetic_closure_kinetics.csv](data/diabetic_closure_kinetics.csv) | Absolute 0 to 100% |
| Diabetic inflammation timecourse | [diabetic_inflammation_timecourse.csv](data/diabetic_inflammation_timecourse.csv) | Peak = 1.0 |
| Diabetic immune cell kinetics | [diabetic_immune_cell_kinetics.csv](data/diabetic_immune_cell_kinetics.csv) | Peak = 1.0 per type |

## Metrics

This module does not add its own metrics columns. Its effects are visible through existing columns: elevated `mean_infl_wound`, increased `n_neutrophils`/`n_macrophages`, reduced `mean_collagen_wound`, delayed `wound_closure_pct`.

## Source files

| File | Purpose |
|------|---------|
| `baseline_inflammation.h` | BaselineInflammationOp for AGE/RAGE chronic inflammation |
| `config.toml` | Module configuration with all scaling factors |
