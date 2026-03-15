> [Home](../../README.md) / [Modules](../README.md) / Lactate

# Lactate

Anaerobic glycolysis metabolite that accumulates in hypoxic wound tissue. Wound lactate levels reach 4 to 12 mM, approximately 3 to 5 times higher than arterial blood (Trabold et al. 2003). Lactate acts as an angiogenesis signal by stabilizing HIF-1alpha (boosting VEGF production) and directly stimulates collagen synthesis by fibroblasts.

## Biology

When tissue oxygen drops below a critical threshold, cells shift from aerobic oxidative phosphorylation to anaerobic glycolysis, producing lactate as an end product. In wounds, the avascular wound bed creates a steep oxygen gradient: the wound center is profoundly hypoxic, and lactate accumulates there to concentrations several-fold higher than in perfused tissue.

Rather than being merely a metabolic waste product, wound lactate serves important signaling functions. It stabilizes hypoxia-inducible factor 1-alpha (HIF-1alpha) even under normoxic conditions, which boosts vascular endothelial growth factor (VEGF) production in wound fibroblasts and macrophages (Constant et al. 2000). This lactate-HIF-1alpha-VEGF axis provides a metabolic signal for angiogenesis that complements the direct hypoxia-driven VEGF pathway.

Lactate also directly stimulates collagen synthesis by fibroblasts, approximately doubling the deposition rate at wound-relevant concentrations, through a redox-dependent mechanism independent of hypoxia itself (Hunt et al. 2007). This ensures that the proliferative and remodeling phases proceed robustly even in the metabolically stressed wound environment.

Lactate production couples to glucose availability: anaerobic glycolysis consumes glucose as substrate. When the glucose module is active, the lactate production rate scales with local glucose concentration. This means diabetic hyperglycemia produces more lactate under equal hypoxia, reflecting the increased glycolytic substrate available.

As vascular perfusion is restored through angiogenesis, lactate is cleared via venous washout, providing a natural feedback loop that resolves the metabolic signal as the wound heals.

## Model

Single diffusing continuum field representing normalized tissue lactate concentration.

**Lactate production:**
- Produced in wound voxels where O2 falls below `o2_threshold`
- Production rate scales linearly with hypoxia fraction: `(threshold - O2) / threshold`
- When glucose module is enabled, production also scales with local glucose concentration (anaerobic glycolysis requires glucose substrate)
- Active in both epidermal wound and dermal wound voxels

**Lactate clearance:**
- Background PDE decay (systemic clearance, maintains equilibrium)
- Perfusion-driven washout proportional to local vascular density (only where perfusion exceeds 0.1)

**Downstream signaling (via source hook):**
- HIF-1alpha stabilization: boosts VEGF production proportional to `vegf_boost * lactate * base_vegf_rate`, applied where O2 is below the VEGF production threshold
- Collagen boost: read by fibroblast behavior to scale collagen deposition by `(1 + collagen_boost * lactate)`

## Parameters

From modules/lactate/config.toml and params.h:

| Parameter | Default | Units | Description | Source |
|-----------|---------|-------|-------------|--------|
| `enabled` | true | bool | Master switch | Convention |
| `diffusion` | 0.01 | - | Tissue diffusion coefficient | Calibrated |
| `decay` | 0.02 | per step | Systemic clearance (maintains equilibrium) | Calibrated |
| `production_rate` | 0.003 | per step | Anaerobic glycolysis rate under hypoxia | Trabold et al. 2003 ([DOI](https://doi.org/10.1046/j.1524-475x.2003.11621.x)) |
| `o2_threshold` | 0.3 | normalized | O2 below this triggers lactate production | Trabold et al. 2003 ([DOI](https://doi.org/10.1046/j.1524-475x.2003.11621.x)) |
| `vegf_boost` | 0.03 | multiplier | HIF-1alpha stabilization boost to VEGF production | Constant et al. 2000 ([DOI](https://doi.org/10.1046/j.1524-475x.2000.00353.x)) |
| `collagen_boost` | 0.03 | multiplier | Collagen synthesis enhancement factor for fibroblasts | Hunt et al. 2007 ([DOI](https://doi.org/10.1089/ars.2007.1674)) |
| `perfusion_clearance` | 0.005 | per step | Perfusion-driven lactate venous washout rate | Calibrated |

## Coupling

### Reads
| Field | Source module | How used |
|-------|-------------|----------|
| O2 | oxygen | Hypoxia fraction drives lactate production (below `o2_threshold`) |
| Glucose | glucose | Scales lactate production by glucose availability (anaerobic glycolysis substrate) |
| Vascular | angiogenesis | Perfusion-driven lactate clearance (venous washout) |
| VEGF | angiogenesis | Base VEGF production rate read from SignalBoard for HIF-1alpha boost calculation |

### Writes
| Field | Consumer modules | What is written |
|-------|-----------------|-----------------|
| Lactate | fibroblast (collagen boost) | Hypoxia-driven production, perfusion clearance |
| VEGF | angiogenesis | HIF-1alpha stabilization boosts VEGF production in wound voxels |

## Source files

| File | Purpose |
|------|---------|
| `lactate_pde.h` | Lactate diffusion field: initialization (starts at 0, builds from hypoxia) |
| `source_hook.h` | Hypoxia-driven lactate production, perfusion clearance, HIF-1alpha VEGF boost |
| `params.h` | LactateParams struct |
| `config.toml` | Module configuration |
