> [Home](../../README.md) / [Modules](../README.md) / Temperature

# Temperature

Wound bed thermal regulation modeling the interplay between vascular perfusion warming and surface evaporative cooling. Healthy skin is maintained at core body temperature (37C, normalized to 1.0) by vascular perfusion. Wounding exposes tissue to the environment, cooling the wound surface. Temperature affects cell migration speed, enzyme kinetics, and bacterial growth via Q10 coefficients.

## Biology

Skin temperature is maintained near core body temperature (approximately 37C) by blood flow through the dermal vascular plexus. When a wound disrupts this vasculature, the exposed wound surface cools via evaporation and radiative heat loss, dropping to approximately 33 to 35C depending on wound depth and ambient conditions (Fierheller and Sibbald 2010).

Temperature recovery tracks vascular restoration: as angiogenesis rebuilds the capillary network, perfusion-driven warming gradually returns the wound bed toward body temperature (McGuiness et al. 2004). The wound surface remains cooler than deeper tissue because evaporative cooling acts at the air-tissue interface while perfusion warming comes from below.

Temperature affects biological rates via the Q10 temperature coefficient, which describes how a 10C change alters reaction rates (Mundim et al. 2020). In the wound context, three Q10-sensitive processes are modeled: MMP enzymatic activity (cooler wounds have slower matrix remodeling), bacterial growth (warmer wounds promote faster biofilm expansion), and optionally cell migration and proliferation (disabled by default for normal wounds, active in diabetic and biofilm scenarios).

The photon module also couples to temperature: absorbed photon energy heats tissue locally, with the thermal coupling scaled by material conductivity.

## Model

Single diffusing continuum field representing normalized tissue temperature (37C = 1.0).

**Initialization and wounding:**
- All voxels initialized to 1.0 (body temperature)
- Wounding sets epidermal wound voxels to `wound_surface` temperature
- Dermal wound voxels receive depth-dependent partial cooling: deeper tissue retains more warmth via `wound_surface + (1.0 - wound_surface) * depth_factor`

**Temperature sources and sinks (via source hook):**
- Dermal perfusion warming: wound voxels warm toward 1.0 at `perfusion_warming * vascular_fraction`
- Surface cooling: epidermal wound voxels lose heat proportional to exposed surface area (1 minus barrier), capped at `wound_surface` floor
- Perfusion warming from below: epidermal voxels receive warming from the vascular plexus beneath them

**Q10 coupling (consumed by other modules):**
- Temperature is denormalized to Celsius as `temp_val * 37.0`
- Rate scaling: `Q10 ^ ((temp_C - 37.0) / 10.0)`
- At wound surface (35C): MMP activity scales to `Q10_mmp ^ (-0.2)`, biofilm growth scales to `Q10_biofilm ^ (-0.2)`
- As perfusion restores temperature toward 37C, Q10 scaling returns toward 1.0

## Parameters

From modules/temperature/config.toml and params.h:

| Parameter | Default | Units | Description | Source |
|-----------|---------|-------|-------------|--------|
| `enabled` | true | bool | Master switch | Convention |
| `diffusion` | 0.001 | - | Thermal conduction through tissue | Calibrated |
| `decay` | 0.0 | per step | No intrinsic decay (energy conservation) | Convention |
| `wound_surface` | 0.946 | normalized | Wound surface temperature (approximately 35C / 37C) | Fierheller and Sibbald 2010 ([DOI](https://doi.org/10.1097/01.asw.0000383197.28192.98)) |
| `perfusion_warming_rate` | 0.025 | per step | Perfusion-driven rewarming rate | McGuiness et al. 2004 ([DOI](https://doi.org/10.12968/jowc.2004.13.9.26702)) |
| `surface_cooling_rate` | 0.002 | per step | Evaporative and radiative heat loss at exposed surface | Calibrated |
| `q10_migration` | 1.0 | coefficient | Cell migration Q10 (disabled for normal wound, active in biofilm and diabetic) | Mundim et al. 2020 ([DOI](https://doi.org/10.1016/j.ecolmodel.2020.109127)) |
| `q10_proliferation` | 1.0 | coefficient | Cell cycle Q10 (disabled for normal wound, active in diabetic) | McGuiness et al. 2004 ([DOI](https://doi.org/10.12968/jowc.2004.13.9.26702)) |
| `q10_mmp` | 1.3 | coefficient | MMP enzymatic activity Q10 | Calibrated |
| `q10_biofilm` | 2.0 | coefficient | Bacterial growth Q10 | Calibrated |

## Coupling

### Reads
| Field | Source module | How used |
|-------|-------------|----------|
| Vascular | angiogenesis | Perfusion-driven warming toward body temperature |
| Barrier | tissue | Surface cooling gated by exposed surface area (1 minus barrier) |

### Writes
| Field | Consumer modules | What is written |
|-------|-----------------|-----------------|
| Temperature | mmp (enzyme kinetics Q10), biofilm (bacterial growth Q10), photon (thermal absorption target) | Perfusion warming, surface cooling |

## Source files

| File | Purpose |
|------|---------|
| `temperature_pde.h` | Temperature diffusion field: initialization to body temp, wound cooling (depth-dependent) |
| `source_hook.h` | Perfusion warming (dermal), surface evaporative cooling (epidermal), sub-surface perfusion warming |
| `params.h` | TemperatureParams struct |
| `config.toml` | Module configuration |
