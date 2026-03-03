> [Home](../../README.md) / [Modules](../README.md) / ROS

# ROS (Reactive Oxygen Species)

Oxidative stress field mediating tissue damage from neutrophil respiratory burst, macrophage oxidative burst, and mitochondrial electron transport chain leak. Provides the mechanistic link between hyperglycemia and impaired healing in diabetic wounds.

## Biology

Reactive oxygen species (superoxide, hydrogen peroxide, hydroxyl radical) are produced during wound healing from three primary sources: (1) neutrophil NADPH oxidase (NOX2) respiratory burst for antimicrobial defense (Babior 1999), (2) M1 macrophage oxidative burst during phagocytosis, and (3) constitutive mitochondrial electron transport chain leak that increases under hypoxia (reverse electron transport at Complex I) and hyperglycemia (electron donor overload at Complex III; Nishikawa et al. 2000).

At physiological levels, ROS functions as a signaling molecule promoting immune defense. At pathological levels (chronic wounds, diabetes), excessive ROS overwhelms the enzymatic antioxidant defense (SOD, catalase, glutathione peroxidase) and causes oxidative tissue damage through multiple pathways:

1. **DNA damage**: 8-oxoguanine formation drives cellular senescence accumulation
2. **Pro-MMP activation**: ROS oxidizes the cysteine switch domain of latent MMPs, converting zymogens to active proteases (Rajagopalan et al. 1996)
3. **NF-kB activation**: ROS amplifies pro-inflammatory signaling
4. **Collagen damage**: oxidative cross-link disruption degrades ECM
5. **Endothelial dysfunction**: ROS impairs VEGFR2 signaling, reducing angiogenic response

In diabetic wounds, Brownlee's unifying hypothesis (2001) identifies mitochondrial superoxide overproduction as the central mediator: hyperglycemia overloads the electron transport chain, generating excess superoxide that activates all major pathways of diabetic tissue damage (polyol, AGE, PKC, hexosamine). This module makes that pathway explicit rather than using parametric multipliers.

## Model

**ROS PDE:** diffusion 0.03, decay 0.04 (fast enzymatic clearance). Three sources:
- Neutrophil respiratory burst: `ros_neutrophil_burst` per neutrophil per step (agent-deposited)
- M1 macrophage oxidative burst: `ros_m1_burst` per M1 macrophage per step (agent-deposited)
- Mitochondrial ETC leak: `ros_mitochondrial_rate` per wound voxel (tissue-wide, in fused_source)

**Hypoxia amplification:** when local O2 < `ros_hypoxia_threshold`, mitochondrial ROS increases up to `ros_hypoxia_amplification` fold (reverse electron transport at Complex I).

**Diabetic amplification:** in diabetic mode, mitochondrial rate multiplied by `diabetic_mito_ros_factor` (hyperglycemia superoxide overproduction) and antioxidant clearance reduced by `diabetic_antioxidant_factor` (impaired enzymatic defense).

**Clearance:**
- Enzymatic (SOD, catalase, GPx): `ros_tissue_antioxidant * tissue_density * ros_concentration` (scales with tissue cellularity; open wound = low clearance)
- Vascular washout: `ros_perfusion_clearance * perfusion * ros_concentration`
- PDE decay: intrinsic scavenging at `decay = 0.04`

**Feedback loop:**
```
Wound --> Neutrophils/M1 --> NADPH oxidase --> ROS
                                                |
                   Mitochondria (hypoxia, glucose) --> ROS
                                                |
              +--- DNA damage --> Senescence --> SASP --> Inflammation
              |                                           |
              +--- Cysteine switch --> Pro-MMP activation --> ECM damage
              |                                           |
              +--- NF-kB --> Inflammation amplification --+
              |
              +--- VEGFR2 dysfunction --> Impaired angiogenesis
              |
              +--- Collagen cross-link damage --> ECM quality loss
```

In diabetic wounds, the hyperglycemia -> mitochondrial superoxide pathway creates a vicious cycle: ROS amplifies inflammation, inflammation recruits more neutrophils/M1 macrophages, which produce more ROS. The impaired antioxidant capacity means this cycle does not resolve normally.

## Parameters

From modules/ros/config.toml:

| Parameter | Default | Units | Description | Source |
|-----------|---------|-------|-------------|--------|
| `enabled` | true | bool | Master switch | Convention |
| `diffusion` | 0.03 | - | ROS diffusion coefficient (small molecule) | Calibrated |
| `decay` | 0.04 | per step | Enzymatic scavenging (SOD, catalase, GPx) | Schafer & Werner 2008 ([DOI](https://doi.org/10.1016/j.phrs.2008.06.004)) |
| `neutrophil_burst` | 0.015 | per step | NADPH oxidase respiratory burst per neutrophil | Babior 1999 ([DOI](https://doi.org/10.1182/blood.v93.5.1464)) |
| `m1_burst` | 0.008 | per step | M1 macrophage oxidative burst | Babior 1999 |
| `mitochondrial_rate` | 0.001 | per step | Basal mitochondrial ETC leak (tissue-wide) | Nishikawa et al. 2000 ([DOI](https://doi.org/10.1038/35008121)) |
| `hypoxia_threshold` | 0.4 | normalized | O2 below this amplifies mitochondrial ROS | Calibrated |
| `hypoxia_amplification` | 2.0 | fold | Max mitochondrial ROS increase under hypoxia | Calibrated |
| `diabetic_mito_factor` | 3.0 | fold | Hyperglycemia superoxide overproduction | Nishikawa et al. 2000 |
| `diabetic_antioxidant_factor` | 0.5 | fraction | Impaired SOD/catalase/GPx capacity | Schafer & Werner 2008 |
| `senescence_ros_rate` | 0.0004 | per step | DNA damage driving senescence accumulation | Calibrated |
| `mmp_activation_rate` | 0.02 | per step | Oxidative pro-MMP cysteine switch activation | Rajagopalan et al. 1996 ([DOI](https://doi.org/10.1172/jci119076)) |
| `angiogenesis_impairment` | 0.15 | fraction | Endothelial VEGFR2 dysfunction | Calibrated |
| `inflammation_amplification` | 0.002 | per step | NF-kB activation amplifies inflammation | Dunnill et al. 2017 ([DOI](https://doi.org/10.1111/iwj.12557)) |
| `collagen_damage_rate` | 0.005 | per step | Oxidative collagen cross-link damage | Calibrated |
| `tissue_antioxidant_rate` | 0.02 | per step | Enzymatic clearance by resident cells | Calibrated |
| `perfusion_clearance` | 0.005 | per step | Vascular washout of ROS byproducts | Calibrated |

## Coupling

### Reads
| Field | Source module | How used |
|-------|-------------|----------|
| O2 | tissue | Hypoxia amplifies mitochondrial ROS production |
| Vascular | perfusion | Tissue density for antioxidant clearance, vascular washout |

### Writes (via fused_post effects)
| Field | Consumer modules | What is written |
|-------|-----------------|-----------------|
| Senescence | senescence | DNA damage accumulation from ROS |
| MMP | mmp | Pro-MMP activation via cysteine switch oxidation |
| Inflammation | inflammation | NF-kB mediated amplification |
| Collagen | fibroblast | Oxidative cross-link damage |
| Vascular recovery | perfusion | Endothelial dysfunction (rate reduction) |

## Validation

| Dataset | Observable | Sources | Notes |
|---------|-----------|--------|-------|
| `ros_parameters` | ROS dynamics, antioxidant kinetics | Schafer 2008, Dunnill 2017, Nishikawa 2000 | Parameter derivation |

## Metrics

| Column | Units | Description |
|--------|-------|-------------|
| `mean_ros_wound` | a.u. | Mean ROS concentration in wound |

## Source files

| File | Purpose |
|------|---------|
| `ros_pde.h` | ROS diffusion field |
| `config.toml` | Module configuration |
