> [Home](../../README.md) / [Modules](../README.md) / Photon

# Photon

Light transport for photobiomodulation (PBM) and optogenetic stimulation, modeling photon diffusion in tissue with opsin kinetics and phototoxicity.

## Biology
Photobiomodulation uses low-level light (typically red or near-infrared wavelengths) to stimulate cellular processes including mitochondrial cytochrome c oxidase activation, ATP production, and wound healing acceleration. Light propagation through tissue is governed by absorption and scattering: melanin, hemoglobin, and water absorb photons while collagen fibers and cell membranes scatter them. The diffusion approximation provides an effective framework for computing the fluence rate (photon density) at depth when scattering dominates over absorption.

Channelrhodopsins are light-gated ion channels that open in response to photon absorption. The two-state model describes the channel as switching between a closed resting state and an open conducting state, with light driving the closed-to-open transition and thermal relaxation driving recovery. The open-state fraction determines the photocurrent and downstream cellular effects.

Excessive light exposure generates reactive oxygen species (ROS) through photosensitization reactions and can cause local tissue heating from absorbed photon energy. These phototoxicity effects set an upper bound on therapeutic light dose.

## Model
Two continuum PDE fields coupled through a source hook:

**Fluence field:** Represents photon density (light intensity) throughout the tissue volume. Uses diffusion approximation transport where D = 1 / (3 * (mu_a + mu_s')), with absorption acting as a decay sink. A Gaussian beam profile injects photons at the wound surface, centered at `(beam_center_x, beam_center_y)` with radius `beam_radius`. Per-voxel material absorption properties enable tissue-appropriate attenuation.

**Opsin field:** Represents channelrhodopsin open-state fraction (0 to `opsin_saturation`). Non-diffusing structural field (opsins are membrane-bound proteins). Kinetics follow the two-state model: dO/dt = k_on * L * (sat - O) - k_off * O, where L is local fluence, O is open fraction, k_on is the light-driven activation rate, and k_off is thermal dark recovery.

**Phototoxicity coupling:** When local fluence exceeds `phototoxicity_threshold`, ROS is generated at `phototoxicity_rate * excess`. Absorbed photon energy also drives local temperature rise via `thermal_coupling`, scaled inversely by material conductivity.

## Parameters
From modules/photon/config.toml:
| Parameter | Default | Units | Description | Source |
|-----------|---------|-------|-------------|--------|
| `enabled` | false | bool | Master switch | Convention |
| `absorption_coeff` | 0.1 | mm^-1 | Tissue absorption coefficient (mu_a) | Jacques 2013 ([DOI](https://doi.org/10.1088/0031-9155/58/11/R37)) |
| `scattering_coeff` | 1.0 | mm^-1 | Reduced scattering coefficient (mu_s') | Jacques 2013 ([DOI](https://doi.org/10.1088/0031-9155/58/11/R37)) |
| `anisotropy` | 0.9 | - | Scattering anisotropy factor (g) | Jacques 2013 ([DOI](https://doi.org/10.1088/0031-9155/58/11/R37)) |
| `irradiance` | 1.0 | normalized | Surface irradiance (0 to 1) | Calibrated |
| `beam_radius` | 0.3 | sim units | Beam spot radius | Mardinly et al. 2018 ([DOI](https://doi.org/10.1038/s41593-018-0139-8)) |
| `beam_center_x` | 0.0 | sim units | Beam center x (0 = wound center) | Convention |
| `beam_center_y` | 0.0 | sim units | Beam center y (0 = wound center) | Convention |
| `diffusion` | 0.3 | - | Effective photon diffusion coefficient | Jacques 2013 ([DOI](https://doi.org/10.1088/0031-9155/58/11/R37)) |
| `decay` | 0.1 | - | Absorption sink (mu_a contribution) | Jacques 2013 ([DOI](https://doi.org/10.1088/0031-9155/58/11/R37)) |
| `opsin_activation_rate` | 0.5 | per step | Light-driven closed-to-open rate | Fenno et al. 2011 ([DOI](https://doi.org/10.1146/annurev-neuro-061010-113817)) |
| `opsin_deactivation_rate` | 0.1 | per step | Thermal dark-state recovery rate | Fenno et al. 2011 ([DOI](https://doi.org/10.1146/annurev-neuro-061010-113817)) |
| `opsin_saturation` | 0.8 | fraction | Max opsin open fraction (photocurrent plateau) | Fenno et al. 2011 ([DOI](https://doi.org/10.1146/annurev-neuro-061010-113817)) |
| `phototoxicity_threshold` | 0.5 | a.u. | Fluence threshold for ROS generation | Mardinly et al. 2018 ([DOI](https://doi.org/10.1038/s41593-018-0139-8)) |
| `phototoxicity_rate` | 0.002 | per step | ROS generation rate above threshold | Mardinly et al. 2018 ([DOI](https://doi.org/10.1038/s41593-018-0139-8)) |
| `thermal_coupling` | 0.001 | - | Absorbed light to temperature rise rate | Packer et al. 2015 ([DOI](https://doi.org/10.1038/nmeth.3217)) |

## Coupling
### Reads
| Field | Source module | How used |
|-------|-------------|----------|
| Fluence | photon (self) | Drives opsin kinetics and phototoxicity |
| Material properties | tissue | Per-voxel absorption and conductivity |

### Writes
| Field | Consumer modules | What is written |
|-------|-----------------|-----------------|
| Fluence | photon (self) | Beam source injection and absorption sink |
| Opsin | (downstream consumers) | Channelrhodopsin open-state fraction |
| ROS | ros | Phototoxicity from excessive fluence |
| Temperature | temperature | Absorbed photon energy converted to heat |

## Source files
| File | Purpose |
|------|---------|
| `fluence_pde.h` | FluencePDE: photon density field with diffusion transport and absorption decay |
| `opsin_pde.h` | OpsinPDE: channelrhodopsin open-state fraction (non-diffusing structural field) |
| `source_hook.h` | PhotonSourceHook: beam injection, absorption, opsin kinetics, phototoxicity, thermal coupling |
| `params.h` | PhotonParams struct |
| `config.toml` | Module configuration |
