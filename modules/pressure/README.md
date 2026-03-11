> [Home](../../README.md) / [Modules](../README.md) / Pressure

# Pressure

Pressure ulcer mechanics modeling ischemia-reperfusion injury from sustained mechanical load, with NPUAP staging and repositioning protocols.

## Biology
Pressure ulcers develop when sustained mechanical load exceeds capillary closing pressure (approximately 32 mmHg), occluding blood flow and causing progressive tissue ischemia. Skeletal muscle is more susceptible to compression damage than skin, so deep tissue injury often initiates at the bone-muscle interface before becoming visible at the surface. The relationship between pressure magnitude and time to tissue necrosis is inverse: high pressures cause damage rapidly while lower pressures require prolonged exposure.

When pressure is relieved, reperfusion of ischemic tissue generates a burst of reactive oxygen species (ROS) through ischemia-reperfusion injury. This ROS burst can paradoxically damage surviving cells and amplify the inflammatory response, extending the injury beyond the original compression zone. Shear stress at the skin surface (from sliding against support surfaces) amplifies the effective compression damage.

Moisture from incontinence or perspiration macerates the skin and increases susceptibility to pressure damage. Regular repositioning (pressure relief cycling) is the primary clinical intervention, allowing intermittent perfusion recovery between compression periods.

The NPUAP staging system classifies pressure injuries by depth: Stage I (non-blanchable erythema, reversible), Stage II (partial thickness skin loss), Stage III (full thickness skin loss into subcutaneous tissue), and Stage IV (full thickness tissue loss exposing muscle or bone).

## Model
Global pressure state with per-voxel field effects:

**Compression state:** A global flag tracks whether the tissue is currently under compression. When `reposition_interval_h > 0`, the compression state alternates on and off at that interval, modeling clinical repositioning protocols. When set to 0, pressure is constant (worst case).

**Tissue damage accumulator:** A global scalar (0 to 1) that increases at `tissue_damage_rate` during compression periods (pressure-time relationship). The accumulator maps directly to NPUAP staging via configurable thresholds.

**Per-voxel effects during compression:**
- Perfusion is reduced by `ischemia_rate * shear_factor` (capillary occlusion)
- Inflammation is generated proportionally to cumulative tissue damage once Stage I threshold is exceeded

**Per-voxel effects during offloading (pressure release):**
- ROS burst proportional to `reperfusion_ros_burst * tissue_damage` (ischemia-reperfusion injury)
- Perfusion recovery toward `1.0 - tissue_damage` (damaged tissue has permanently reduced capacity)

## Parameters
From modules/pressure/config.toml:
| Parameter | Default | Units | Description | Source |
|-----------|---------|-------|-------------|--------|
| `enabled` | false | bool | Master switch | Convention |
| `compression_threshold` | 0.3 | normalized | Compression above this occludes capillaries | Gefen 2007 ([DOI](https://doi.org/10.12968/jowc.2007.16.8.27854)) |
| `ischemia_rate` | 0.02 | per step | Tissue oxygen drop rate under compression | Gefen 2007 ([DOI](https://doi.org/10.12968/jowc.2007.16.8.27854)) |
| `reperfusion_ros_burst` | 0.15 | a.u. | ROS burst magnitude on pressure release | Stekelenburg et al. 2008 ([DOI](https://doi.org/10.1016/j.apmr.2008.01.012)) |
| `shear_factor` | 1.5 | multiplier | Shear stress amplifies compression damage | Gefen 2007 ([DOI](https://doi.org/10.12968/jowc.2007.16.8.27854)) |
| `tissue_damage_rate` | 0.003 | per hour | Tissue damage accumulation under compression | Kosiak 1959, Stekelenburg et al. 2006 ([DOI](https://doi.org/10.1152/japplphysiol.00889.2005)) |
| `reposition_interval_h` | 0 | hours | Hours between pressure relief (0 = constant) | NPUAP/EPUAP 2019 |
| `moisture_damage_rate` | 0.001 | per step | Moisture-associated skin damage | NPUAP/EPUAP 2019 |
| `stage_1` | 0.15 | fraction | Stage I: non-blanchable erythema | NPUAP/EPUAP 2019 |
| `stage_2` | 0.35 | fraction | Stage II: partial thickness skin loss | NPUAP/EPUAP 2019 |
| `stage_3` | 0.60 | fraction | Stage III: full thickness skin loss | NPUAP/EPUAP 2019 |
| `stage_4` | 0.85 | fraction | Stage IV: full thickness tissue loss | NPUAP/EPUAP 2019 |

## Coupling
### Reads
| Field | Source module | How used |
|-------|-------------|----------|
| Water | hydration | Moisture damage assessment |

### Writes
| Field | Consumer modules | What is written |
|-------|-----------------|-----------------|
| Vascular (perfusion) | perfusion, angiogenesis | Reduced during compression, recovers during offloading |
| ROS | ros | Reperfusion burst on pressure release |
| Inflammation | immune, scar | Generated proportional to tissue damage above Stage I |

## Source files
| File | Purpose |
|------|---------|
| `source_hook.h` | PressureSourceHook: ischemia-reperfusion cycle, tissue damage, staging, repositioning |
| `params.h` | PressureParams struct with staging thresholds |
| `config.toml` | Module configuration |
