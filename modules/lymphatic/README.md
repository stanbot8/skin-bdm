> [Home](../../README.md) / [Modules](../README.md) / Lymphatic

# Lymphatic

Lymphangiogenesis and interstitial fluid dynamics, modeling lymphatic vessel regeneration and edema resolution during wound healing.

## Biology
Lymphatic vessels in the dermis drain interstitial fluid, clear cytokines, and transport antigen-presenting cells. Wounding disrupts the local lymphatic network, causing interstitial edema as inflammatory vascular leak (governed by Starling forces) exceeds drainage capacity. Lymphangiogenesis, the regrowth of lymphatic vessels from intact margins, is driven by macrophage-derived VEGF-C/D signaling through VEGFR-3. This process lags angiogenesis by several days and proceeds more slowly than blood vessel regrowth.

Edema has two functional consequences in the model: it reduces oxygen delivery to the wound bed by acting as a diffusion barrier, and it impairs cell migration through hydrostatic resistance. Lymphatic drainage resolves edema by removing interstitial fluid proportionally to local lymphatic vessel density. In diabetic tissue, both lymphatic regrowth rate and target density are reduced, prolonging edema and delaying healing.

The lymphatic module also provides a post-wound TGF-beta drainage pathway. Local lymphatic density determines how much TGF-beta is cleared from the wound, preventing excessive fibrosis signals from accumulating.

## Model
Continuum PDE fields with source/sink dynamics:

**Lymphatic field:** Represents lymphatic vessel density (normalized, 0 to 1). Uses slow diffusion to model sprouting from wound margins. Regenerates toward `basal_density` at `regen_rate`, boosted by local VEGF concentration. In diabetic mode, both the target density and regeneration rate are scaled by `diabetic_lymphatic_factor`.

**Edema field:** Represents interstitial fluid accumulation. Non-diffusing (D=0) since edema is local. Source: inflammatory vascular leak proportional to `edema_leak_rate * inflammation * perfusion`. Sink: lymphatic drainage proportional to `edema_drainage_rate * lymphatic_density * edema`.

**Post hook (TGF-beta drainage):** Local lymphatic density determines TGF-beta clearance rate in epidermal wound voxels. Drainage is proportional to `edema_drainage_rate * lymphatic_density * tgfb_concentration`.

## Parameters
From modules/lymphatic/config.toml:
| Parameter | Default | Units | Description | Source |
|-----------|---------|-------|-------------|--------|
| `enabled` | true | bool | Master switch | Convention |
| `diffusion` | 0.005 | - | Slow lymphatic sprouting diffusion | Paavonen et al. 2000 ([DOI](https://doi.org/10.1016/s0002-9440(10)65021-3)) |
| `basal_density` | 0.8 | normalized | Healthy dermis lymphatic density | Kataru et al. 2009 ([DOI](https://doi.org/10.1182/blood-2008-09-176776)) |
| `regen_rate` | 0.0003 | per step | Regeneration rate (slower than vascular) | Kataru et al. 2009 ([DOI](https://doi.org/10.1182/blood-2008-09-176776)) |
| `vegf_boost` | 0.3 | - | VEGF-C promotes lymphangiogenesis | Makinen et al. 2001 ([DOI](https://doi.org/10.1093/emboj/20.17.4762)) |
| `edema_leak_rate` | 0.003 | per step | Vascular leak from inflammation (Starling forces) | Michel and Curry 1999 ([DOI](https://doi.org/10.1152/physrev.1999.79.3.703)) |
| `edema_drainage_rate` | 0.01 | per step | Lymphatic drainage of interstitial fluid | Calibrated |
| `edema_o2_impairment` | 0.3 | fraction | Edema reduces O2 delivery (diffusion barrier) | Calibrated |
| `edema_migration_impairment` | 0.2 | fraction | Edema slows cell crawling (hydrostatic resistance) | Rutkowski and Swartz 2007 ([DOI](https://doi.org/10.1016/j.tcb.2006.11.007)) |
| `diabetic_lymphatic_factor` | 0.4 | fraction | Impaired lymphangiogenesis in diabetes | Asai et al. 2012 ([DOI](https://doi.org/10.1016/j.ajpath.2012.08.023)) |

## Coupling
### Reads
| Field | Source module | How used |
|-------|-------------|----------|
| VEGF | angiogenesis | Boosts lymphatic regeneration rate |
| Inflammation | inflammation | Drives edema source (vascular leak) |
| Vascular (perfusion) | perfusion | Scales vascular leak rate |
| TGF-beta | immune | Drained by lymphatic post hook |
| Oxygen | oxygen | Reduced by edema O2 impairment |

### Writes
| Field | Consumer modules | What is written |
|-------|-----------------|-----------------|
| Lymphatic | lymphatic (self) | Vessel density regeneration |
| Edema | tissue (migration gating) | Interstitial fluid accumulation and drainage |
| Oxygen | all O2-dependent modules | Reduced by edema diffusion barrier |
| TGF-beta | fibroblast, scar | Drained proportional to lymphatic density |

## Source files
| File | Purpose |
|------|---------|
| `lymphatic_pde.h` | LymphaticPDE: lymphatic vessel density field with slow diffusion |
| `edema_pde.h` | EdemaPDE: interstitial fluid field (non-diffusing, mechanistic source/sink) |
| `source_hook.h` | LymphaticSourceHook: lymphatic regen, edema source/sink, O2 impairment |
| `post_hook.h` | LymphaticPostHook: TGF-beta drainage proportional to lymphatic density |
| `params.h` | LymphaticParams struct |
| `config.toml` | Module configuration |
