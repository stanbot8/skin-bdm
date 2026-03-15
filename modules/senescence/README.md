> [Home](../../README.md) / [Modules](../README.md) / Senescence

# Senescence

Cellular senescence accumulation, SASP (Senescence-Associated Secretory Phenotype) output, and immune-mediated clearance. Transient senescence promotes healing via PDGF-AA; excessive accumulation in diabetic or aged wounds impairs repair through chronic inflammatory SASP.

## Biology

Cellular senescence is a permanent growth arrest triggered by DNA damage, oxidative stress, or oncogene activation. In acute wounds, senescent fibroblasts and endothelial cells accumulate transiently and secrete PDGF-AA, which promotes myofibroblast differentiation and optimal wound closure (Demaria et al. 2014). This beneficial transient senescence is cleared by natural killer (NK) cells that recognize NKG2D ligands (MICA/MICB) upregulated on senescent cell surfaces (Sagiv et al. 2016), along with macrophage-mediated efferocytosis (Kang et al. 2011).

In diabetic or aged wounds, senescent cells accumulate excessively due to hyperglycemia-driven AGE (Advanced Glycation End-product) stress, chronic inflammation via ROS and NF-kB, and impaired immune clearance. The resulting persistent SASP creates a pro-inflammatory milieu rich in IL-6, IL-8, MMPs, and TGF-beta that sustains chronic inflammation and drives fibrosis (Coppe et al. 2008; He and Sharpless 2017). Senescence also limits fibrosis by arresting activated myofibroblasts (Krizhanovsky et al. 2008), creating a double-edged role in tissue repair.

## Model

**Senescence field:** Non-diffusing structural field representing local senescent cell density. Accumulates from three sources and is cleared by two pathways:

```
Accumulation sources:
  1. Wound-induced DNA damage (gamma-H2AX)        [wound_accumulation_rate * wound_gate]
  2. Inflammation-driven (ROS, NF-kB activation)   [inflammation_accumulation * local_inflammation]
  3. AGE-driven glycative stress (diabetic mode)    [age_accumulation * local_AGE]

Clearance pathways:
  1. Immune surveillance (NK + macrophages)         [immune_clearance_rate * senescence * immune_pressure]
  2. Senolytic therapy (treatment parameter)        [senolytic_clearance_rate * senescence]
```

**SASP output:** When senescent cell density exceeds threshold, three secretory outputs are produced proportional to local senescence concentration:
- Pro-inflammatory cytokines (IL-6, IL-8) deposited to the Inflammation field
- MMP-3/MMP-9 deposited to the ProMMP field as zymogens
- TGF-beta1 deposited to the TGF-beta field (fibrotic component)

**Diabetic mode:** A multiplicative acceleration factor (default 2.5x) represents hyperglycemia-accelerated senescence via Brownlee pathways (polyol, AGE, PKC, hexosamine).

## Parameters

From modules/senescence/config.toml:

| Parameter | Default | Units | Description | Source |
|-----------|---------|-------|-------------|--------|
| `enabled` | true | bool | Master switch | Convention |
| `diffusion` | 0.0 | - | Immobile (senescent cells are growth-arrested) | Convention |
| `decay` | 0.00005 | per step | Very slow clearance (immune-mediated, days to weeks) | Calibrated |
| `wound_accumulation_rate` | 0.00002 | per step | Basal wound-induced DNA damage (gamma-H2AX) | Demaria et al. 2014 ([DOI](https://doi.org/10.1016/j.devcel.2014.11.012)) |
| `inflammation_accumulation` | 0.00005 | per step | Inflammation-driven senescence (ROS, NF-kB) | Campisi 2007 ([DOI](https://doi.org/10.1038/nrm2233)) |
| `age_accumulation` | 0.0002 | per step | AGE-driven senescence in diabetic mode | He and Sharpless 2017 ([DOI](https://doi.org/10.1016/j.cell.2017.05.015)) |
| `sasp_inflammation_rate` | 0.00002 | per step | SASP pro-inflammatory cytokine output (IL-6, IL-8) | Coppe et al. 2008 ([DOI](https://doi.org/10.1371/journal.pbio.0060301)) |
| `sasp_mmp_rate` | 0.000005 | per step | SASP MMP-3/MMP-9 production (zymogen) | Coppe et al. 2008 ([DOI](https://doi.org/10.1371/journal.pbio.0060301)) |
| `sasp_tgfb_rate` | 0.000005 | per step | SASP TGF-beta1 (fibrotic component) | Coppe et al. 2008 ([DOI](https://doi.org/10.1371/journal.pbio.0060301)) |
| `immune_clearance_rate` | 0.002 | per step | NK/macrophage-mediated clearance via NKG2D ligands | Sagiv et al. 2016 ([DOI](https://doi.org/10.18632/aging.100897)) |
| `senolytic_clearance_rate` | 0.0 | per step | ABT-263/dasatinib+quercetin clearance (0 = no treatment) | Xu et al. 2018 ([DOI](https://doi.org/10.1038/s41591-018-0092-9)) |
| `diabetic_accumulation_factor` | 2.5 | multiplier | Hyperglycemia accelerates senescence | He and Sharpless 2017 ([DOI](https://doi.org/10.1016/j.cell.2017.05.015)) |

## Coupling

### Reads
| Field | Source module | How used |
|-------|-------------|----------|
| Inflammation | inflammation | Drives inflammation-dependent senescence accumulation |
| AGE | glucose | AGE-driven glycative stress accumulation (diabetic mode) |
| Immune Pressure | immune | Scales immune-mediated senescence clearance (NK + macrophage density proxy) |

### Writes
| Field | Consumer modules | What is written |
|-------|-----------------|-----------------|
| Senescence | (self) | Net senescence density (accumulation minus clearance) |
| Inflammation | immune, fibroblast | SASP pro-inflammatory cytokines (IL-6, IL-8) |
| ProMMP | mmp | SASP MMP-3/MMP-9 as zymogens |
| TGF-beta | fibroblast | SASP TGF-beta1 (fibrotic) |

## Source files

| File | Purpose |
|------|---------|
| `senescence_pde.h` | Non-diffusing structural field for senescent cell density |
| `post_hook.h` | Accumulation (DNA damage, inflammation, AGE), SASP output, immune clearance |
| `senescence_op.h` | Organizational header (logic lives in fused_post.h) |
| `params.h` | Parameter struct (SenescenceParams) |
| `config.toml` | Module configuration |
