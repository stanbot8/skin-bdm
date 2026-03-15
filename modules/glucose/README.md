> [Home](../../README.md) / [Modules](../README.md) / Glucose

# Glucose

Interstitial glucose transport and advanced glycation end-product (AGE) formation in wound tissue. Healthy tissue maintains a normalized glucose level of 1.0 via vascular perfusion. In diabetic tissue, hyperglycemia elevates the basal level, driving non-enzymatic glycation of extracellular matrix proteins and sustained RAGE-mediated inflammation.

## Biology

Glucose is the primary metabolic substrate for all wound cells. It reaches the interstitial space via capillary perfusion and is consumed by resident and infiltrating cells. In healthy acute wounds, glucose supply temporarily drops in the wound bed due to vascular disruption, then recovers as angiogenesis restores perfusion.

In diabetic patients, chronic hyperglycemia elevates interstitial glucose to approximately 1.5 to 2 times the normal fasting level (Brem and Tomic-Canic 2007). The excess glucose reacts non-enzymatically with lysine and arginine residues on extracellular matrix proteins via the Maillard reaction (Schiff base formation followed by Amadori rearrangement), producing irreversible AGEs (Brownlee 2005). AGEs accumulate in the ECM because cross-linked glycated proteins resist normal enzymatic turnover.

AGEs exert their pathological effects by binding the Receptor for Advanced Glycation End-products (RAGE) on macrophages and endothelial cells, activating NF-kB and sustaining pro-inflammatory cytokine production (Bierhaus et al. 2005). This creates a feed-forward loop: hyperglycemia produces AGEs, AGEs drive chronic inflammation via RAGE, and inflammation impairs the healing process. AGE cross-linking of collagen also renders it resistant to MMP-mediated remodeling (Mott et al. 1997), while glycation of fibronectin disrupts integrin binding and impairs fibroblast migration (Rana et al. 2007).

Bacteria in wound biofilm consume glucose as a metabolic substrate, so hyperglycemia also fuels accelerated biofilm growth (Bowler et al. 2001).

## Model

Two coupled continuum fields: Glucose (diffusing metabolic substrate) and AGE (non-diffusing accumulator).

**Glucose field:**
- Initialized everywhere to `basal_conc` (normalized healthy fasting level)
- Wounding drops wound bed glucose to 50% of basal (vascular disruption)
- Perfusion supply: dermal voxels receive glucose proportional to local vascular density, driving recovery toward `basal_conc * vascular_fraction`
- Background decay represents baseline cellular metabolic consumption
- Bacterial consumption removes glucose where biofilm is present

**AGE field:**
- Non-diffusing structural accumulator (diffusion = 0, prescribed grid)
- AGE formation rate is proportional to glucose squared (second-order kinetics: doubling glucose quadruples AGE formation, reflecting the non-enzymatic nature of the Maillard reaction)
- Very slow decay (`age_decay = 0.0001`) representing the near-irreversible nature of cross-linked glycated proteins
- Only active in diabetic mode (`diabetic.mode = true`)

**AGE downstream effects (via post hook):**
- RAGE-mediated NF-kB inflammation: deposits inflammation proportional to local AGE concentration, gated by wound openness (1 minus stratum)
- Collagen cross-linking: blocks a fraction of MMP-mediated collagen degradation proportional to local AGE
- M1 prolongation: extends effective M1 macrophage duration via RAGE signaling
- Fibroblast migration impairment: glycated ECM reduces fibroblast migration speed

## Parameters

From modules/glucose/config.toml and params.h:

| Parameter | Default | Units | Description | Source |
|-----------|---------|-------|-------------|--------|
| `enabled` | true | bool | Master switch | Convention |
| `diffusion` | 0.005 | - | Interstitial glucose diffusion coefficient | Brem and Tomic-Canic 2007 ([DOI](https://doi.org/10.1172/JCI32169)) |
| `decay` | 0.0002 | per step | Baseline cellular uptake (metabolic consumption) | Calibrated |
| `basal_conc` | 1.0 | normalized | Healthy fasting glucose level (1.0 = 5 mM) | Convention |
| `perfusion_supply` | 0.008 | per step | Perfusion-driven glucose delivery rate | Calibrated |
| `age_rate` | 0.0001 | per step | AGE formation rate from glucose (Maillard reaction) | Brownlee 2005 ([DOI](https://doi.org/10.2337/diabetes.54.6.1615)) |
| `age_inflammation` | 0.002 | per step | AGE-driven inflammation per step (legacy path) | Bierhaus et al. 2005 ([DOI](https://doi.org/10.1007/s00109-005-0688-7)) |
| `bacterial_consumption` | 0.01 | per step | Biofilm glucose consumption rate | Calibrated |
| `prolif_threshold` | 0.05 | normalized | Minimum glucose for full proliferation rate | Calibrated |
| `age_decay` | 0.0001 | per step | AGE turnover (very slow; cross-linked proteins persist months) | Brownlee 2005 ([DOI](https://doi.org/10.2337/diabetes.54.6.1615)) |
| `age_rage_inflammation` | 0.001 | per step | RAGE-mediated NF-kB inflammation per AGE unit per step | Bierhaus et al. 2005 ([DOI](https://doi.org/10.1007/s00109-005-0688-7)) |
| `age_collagen_crosslink` | 0.15 | fraction | Fraction of MMP collagen degradation blocked by local AGE (0 to 1) | Mott et al. 1997 ([DOI](https://doi.org/10.1038/ki.1997.455)) |
| `age_m1_prolongation` | 0.05 | multiplier | AGE-RAGE extends effective M1 duration (multiplier on AGE level) | Chavakis et al. 2004 ([DOI](https://doi.org/10.1016/j.micinf.2004.08.004)) |
| `age_fibroblast_migration_impair` | 0.1 | multiplier | Glycated ECM reduces fibroblast migration speed | Rana et al. 2007 ([DOI](https://doi.org/10.1109/nebc.2007.4413357)) |

## Coupling

### Reads
| Field | Source module | How used |
|-------|-------------|----------|
| Vascular | angiogenesis | Perfusion-driven glucose supply rate |
| Stratum | tissue | Gates AGE inflammation to open wound (1 minus stratum) |

### Writes
| Field | Consumer modules | What is written |
|-------|-----------------|-----------------|
| Glucose | lactate (anaerobic glycolysis substrate), biofilm (bacterial metabolism), diabetic (inflammation scaling), fibroblast (proliferation gating) | Perfusion supply, wound depletion |
| AGE | immune (M1 prolongation), mmp (collagen crosslink resistance), senescence (glycative stress), fibroblast (migration impairment) | Non-enzymatic glycation from glucose (diabetic mode only) |
| Inflammation | immune, scar | RAGE-mediated NF-kB inflammation from AGE accumulation |

## Source files

| File | Purpose |
|------|---------|
| `glucose_pde.h` | Glucose diffusion field: initialization, wound depletion |
| `age_pde.h` | AGE non-diffusing accumulator field (prescribed grid) |
| `source_hook.h` | Perfusion-driven glucose supply in dermal voxels |
| `post_hook.h` | AGE formation from glucose (Maillard reaction) and RAGE-mediated inflammation |
| `params.h` | GlucoseParams struct |
| `config.toml` | Module configuration |
