> [Home](../../README.md) / [Modules](../README.md) / Body Site

# Body Site

Anatomical site-specific parameter scaling, adjusting perfusion, nerve density, and temperature for different body regions.

## Biology
Skin properties vary substantially across the body. The scalp has rich blood supply from external carotid branches and dense hair follicles whose bulge stem cells contribute to re-epithelialization, leading to faster wound healing. The plantar foot has very thick glabrous (hairless) epidermis, reduced perfusion due to its distal location and weight-bearing pressure, and cooler surface temperature. The face has thin skin with high vascularity and nerve density. These regional differences affect wound healing speed, infection risk, and scar formation.

Hair follicle density is particularly important because follicular bulge stem cells migrate radially into wounded epidermis and contribute to re-epithelialization. Sites with high follicle density (scalp, face) heal faster than follicle-sparse sites (shin, foot plantar), independent of other factors.

The module uses a site TOML overlay approach: each body site has a TOML file that overrides specific parameters from the baseline (forearm), which is the reference site for most wound healing studies. Only parameters that differ from baseline need to be specified.

## Model
Parameter overlay system (no PDE fields or source hooks):

**Site selection:** Enabled via `--site <name>` on the batch command line. The site TOML file from `modules/body_site/sites/<name>.toml` is merged into the simulation configuration, overriding the relevant baseline values.

**Available sites:**
| Site | Key characteristics |
|------|-------------------|
| `forearm` | Reference baseline (volar forearm, most study data) |
| `face` | Thin skin, high vascularity and nerve density |
| `scalp` | Very high perfusion (1.5x), dense follicles, fast healing |
| `torso` | Average properties |
| `shin` | Lower perfusion, slower healing |
| `hand_dorsum` | Thin skin, moderate vascularity |
| `foot_dorsum` | Distal, reduced perfusion and temperature |
| `foot_plantar` | Thick glabrous skin, no follicles, low perfusion (0.6x), cool (31C) |

**Overridden parameters by site:**
- `skin.perfusion.basal`: regional blood supply (0.6 for plantar foot to 1.5 for scalp)
- `skin.neuropathy.basal_nerve_density`: sensory innervation density
- `skin.temperature.wound_surface`: local surface temperature (normalized)

## Parameters
From modules/body_site/config.toml:
| Parameter | Default | Units | Description | Source |
|-----------|---------|-------|-------------|--------|
| `enabled` | false | bool | Master switch (opt-in) | Convention |
| `site` | "forearm" | string | Anatomical site identifier | Convention |

Site-specific overrides are defined in individual TOML files under `sites/`. Skin thickness data from Lee and Hwang 2002 ([DOI](https://doi.org/10.1007/s00276-002-0034-5)). Follicular contribution to healing from Ito et al. 2005 ([DOI](https://doi.org/10.1038/nm1328)) and Rittie 2016 ([DOI](https://doi.org/10.1007/s12079-016-0330-1)).

## Coupling
### Reads
No fields read directly. The module operates as a configuration overlay.

### Writes
No fields written directly. The module modifies configuration parameters that are consumed by other modules:
| Parameter | Consumer modules | Effect |
|-----------|-----------------|--------|
| `perfusion.basal` | perfusion, angiogenesis | Regional blood supply scaling |
| `neuropathy.basal_nerve_density` | neuropathy | Sensory innervation density |
| `temperature.wound_surface` | temperature | Local surface temperature |

## Source files
| File | Purpose |
|------|---------|
| `config.toml` | Module configuration (enabled flag and site selector) |
| `sites/forearm.toml` | Forearm baseline (reference site, no overrides) |
| `sites/face.toml` | Face site overrides |
| `sites/scalp.toml` | Scalp site overrides (high perfusion, dense follicles) |
| `sites/torso.toml` | Torso site overrides |
| `sites/shin.toml` | Shin site overrides |
| `sites/hand_dorsum.toml` | Hand dorsum site overrides |
| `sites/foot_dorsum.toml` | Foot dorsum site overrides |
| `sites/foot_plantar.toml` | Foot plantar site overrides (diabetic foot ulcer target) |
