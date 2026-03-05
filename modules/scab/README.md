> [Home](../../README.md) / [Modules](../README.md) / Scab

# Scab

Eschar and wound crust formation, modeling the dried fibrin/platelet/exudate layer that protects the wound surface and its gradual dissolution during healing.

## Biology
The scab (eschar) forms rapidly at the wound surface from coagulation cascade products: fibrin, platelets, dried exudate, and trapped red blood cells polymerize into a protective crust. This layer serves as a mechanical barrier against infection and desiccation, reducing evaporative water loss from the wound bed. However, the scab also impedes keratinocyte migration during re-epithelialization, as migrating cells must tunnel underneath the crust rather than gliding freely over a moist surface. This is the basis of moist wound healing: keeping the wound surface hydrated prevents scab formation and accelerates epithelial closure.

Scab dissolution occurs through multiple mechanisms. Baseline shedding removes the outermost layers over days. Matrix metalloproteinases (MMPs) enzymatically degrade the fibrin and collagen matrix of the scab. As keratinocytes re-epithelialize underneath, they physically undermine the scab from below, causing it to lift. Moisture softens the scab and accelerates its breakdown, which is why occlusive dressings speed healing.

## Model
Structural PDE field (non-diffusing, D=0) with multi-pathway decay:

**Wound seeding:** When a wound is created, epidermal wound voxels (from the surface down to the cornified layer) are seeded with scab at `wound_seed` density, representing immediate coagulation cascade activation.

**Decay mechanisms (all multiplicative with current scab density):**
1. Baseline turnover at `decay` rate (shedding, drying, mechanical loss; half-life approximately 3.6 days)
2. Re-epithelialization undermining at `reepith_rate`, proportional to local stratum corneum recovery (barrier fraction)
3. MMP proteolysis at `mmp_degradation * local_mmp_concentration`
4. Moisture softening at `moisture_softening * local_water_concentration`

**Functional effects (consumed by other modules):**
- `evaporation_shield`: scab reduces surface water loss (0 = no shield, 1 = full)
- `migration_penalty`: scab impedes keratinocyte migration speed (0 = no penalty, 1 = full block)

## Parameters
From modules/scab/config.toml:
| Parameter | Default | Units | Description | Source |
|-----------|---------|-------|-------------|--------|
| `enabled` | false | bool | Master switch | Convention |
| `wound_seed` | 1.0 | normalized | Initial scab density in wound | Singer and Clark 1999 ([DOI](https://doi.org/10.1056/NEJM199909023411006)) |
| `decay` | 0.0008 | per step | Baseline turnover (half-life approximately 3.6 days) | Singer and Clark 1999 ([DOI](https://doi.org/10.1056/NEJM199909023411006)) |
| `mmp_degradation` | 0.001 | per step | MMP-mediated scab breakdown | Martin 1997 ([DOI](https://doi.org/10.1126/science.276.5309.75)) |
| `reepith_rate` | 0.005 | per step | Degradation proportional to barrier recovery | Singer and Clark 1999 ([DOI](https://doi.org/10.1056/NEJM199909023411006)) |
| `moisture_softening` | 0.002 | per step | Moisture accelerates scab softening | Winter 1962 ([DOI](https://doi.org/10.1038/193293a0)) |
| `evaporation_shield` | 0.5 | fraction | Scab reduces surface water loss | Martin 1997 ([DOI](https://doi.org/10.1126/science.276.5309.75)) |
| `migration_penalty` | 0.15 | fraction | Scab impedes keratinocyte migration | Winter 1962 ([DOI](https://doi.org/10.1038/193293a0)) |

## Coupling
### Reads
| Field | Source module | How used |
|-------|-------------|----------|
| Stratum (corneum) | tissue | Re-epithelialization undermining (barrier fraction) |
| MMP | mmp | Enzymatic scab degradation |
| Water | hydration | Moisture softening of scab |

### Writes
| Field | Consumer modules | What is written |
|-------|-----------------|-----------------|
| Scab | tissue (migration penalty), hydration (evaporation shield) | Scab density (decay only, no production after wound seeding) |

## Source files
| File | Purpose |
|------|---------|
| `scab_pde.h` | ScabPDE: structural field with wound seeding in epidermal voxels |
| `source_hook.h` | ScabDecayHook: multi-pathway decay (baseline, re-epithelialization, MMP, moisture) |
| `params.h` | ScabParams struct |
| `config.toml` | Module configuration |
