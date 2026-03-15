> [Home](../../README.md) / [Modules](../README.md) / Nitric Oxide

# Nitric Oxide

Nitric oxide (NO) produced by inducible nitric oxide synthase (iNOS) in M1 macrophages and neutrophils during the inflammatory phase of wound healing. NO has a very short half-life (seconds to minutes in tissue), modeled as a high decay rate. It serves three functional roles: vasodilation (increases local perfusion), antimicrobial activity (suppresses biofilm growth via reactive nitrogen species), and anti-fibrotic signaling (reduces excessive collagen deposition).

## Biology

Nitric oxide is a small gaseous signaling molecule produced by the inducible isoform of nitric oxide synthase (iNOS) in activated immune cells. During the inflammatory phase, M1 macrophages are the primary source, with neutrophils contributing at a lower rate (Witte and Barbul 2002). M2 macrophages do not express iNOS and therefore do not produce NO.

NO mediates vasodilation by activating soluble guanylyl cyclase in vascular smooth muscle cells, increasing cyclic GMP and causing relaxation. In wounds, this increases local perfusion and oxygen delivery to the healing tissue (Luo and Chen 2005). NO also serves as a direct antimicrobial effector: reactive nitrogen species (peroxynitrite, nitrogen dioxide) generated from NO are toxic to bacteria and suppress biofilm growth (Fang 1997).

In the later phases of healing, NO exerts an anti-fibrotic effect by reducing collagen deposition. NO donors decrease collagen synthesis and scar formation in wound models (Schaffer et al. 1996), providing a natural brake on excessive fibrosis during the inflammatory phase.

In diabetic wounds, iNOS expression is reduced by approximately 40%, contributing to impaired healing through diminished vasodilation, weakened antimicrobial defense, and loss of anti-fibrotic regulation (Witte et al. 2002).

## Model

Single diffusing continuum field with rapid decay representing the short-lived nature of NO in tissue.

**NO production:**
- M1 macrophages produce NO at `m1_production` rate, tapered by M1 state age (declining iNOS expression over time)
- Neutrophils produce NO at `neutrophil_production` rate, tapered by cell age
- In diabetic mode, production rates are multiplied by `diabetic.no_factor` (approximately 0.6, representing 40% reduction in iNOS expression)
- Production uses agent-level deposition via ScaledGrid

**NO decay:**
- High PDE decay rate (0.05) models rapid autoxidation (seconds to minutes half-life)
- Fast diffusion (0.02) reflects the small gaseous molecule crossing cell membranes freely

**Downstream effects (consumed by other modules):**
- Vasodilation: angiogenesis source hook boosts effective perfusion by `angio_rate * vasodilation_factor * NO` in wound voxels
- Antimicrobial: biofilm post hook scales effective growth rate by `max(0, 1 - antimicrobial_factor * NO)`
- Anti-fibrotic: fibroblast behavior scales collagen deposition by `max(0, 1 - collagen_suppression * NO)`

## Parameters

From modules/nitric_oxide/config.toml and sim_param.h:

| Parameter | Default | Units | Description | Source |
|-----------|---------|-------|-------------|--------|
| `enabled` | true | bool | Master switch | Convention |
| `diffusion` | 0.02 | - | Fast diffusion (small gaseous molecule) | Calibrated |
| `decay` | 0.05 | per step | Rapid autoxidation (half-life seconds to minutes) | Calibrated |
| `m1_production` | 0.01 | per step | iNOS-derived NO from M1 macrophages | Witte and Barbul 2002 ([DOI](https://doi.org/10.1016/s0002-9610(02)00815-2)) |
| `neutrophil_production` | 0.005 | per step | iNOS-derived NO from neutrophils | Witte and Barbul 2002 ([DOI](https://doi.org/10.1016/s0002-9610(02)00815-2)) |
| `vasodilation_factor` | 0.05 | multiplier | Perfusion boost per unit NO | Luo and Chen 2005 ([DOI](https://doi.org/10.1111/j.1745-7254.2005.00058.x)) |
| `antimicrobial_factor` | 0.3 | multiplier | Biofilm growth suppression per unit NO | Calibrated |
| `collagen_suppression` | 0.03 | multiplier | Anti-fibrotic collagen deposition reduction | Schaffer et al. 1996 ([DOI](https://doi.org/10.1006/jsre.1996.0254)) |

## Coupling

### Reads
| Field | Source module | How used |
|-------|-------------|----------|
| (Agent positions) | immune | M1 macrophages and neutrophils deposit NO at their locations via iNOS |

### Writes
| Field | Consumer modules | What is written |
|-------|-----------------|-----------------|
| NitricOxide | angiogenesis (vasodilation boost), biofilm (antimicrobial suppression), fibroblast (collagen suppression) | iNOS production from M1 macrophages and neutrophils |

## Source files

| File | Purpose |
|------|---------|
| `nitric_oxide_pde.h` | NO diffusion field: fast diffusion, high decay (starts at 0, produced by immune cells) |
| `params.h` | NitricOxideParams struct (enabled flag; remaining params in sim_param.h) |
| `config.toml` | Module configuration |
