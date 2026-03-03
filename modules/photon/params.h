#ifndef PHOTON_PARAMS_H_
#define PHOTON_PARAMS_H_

#include "biodynamo.h"

namespace bdm {
namespace skibidy {

struct PhotonParams {

  // Photon transport: diffusion approximation for light in tissue
  // Jacques 2013 (doi:10.1088/0031-9155/58/11/R37)
  bool enabled = false;
  real_t absorption_coeff = 0.1;  // mu_a (mm^-1), brain grey matter at 920 nm
  real_t scattering_coeff = 1.0;  // mu_s' reduced scattering (mm^-1)
  real_t anisotropy = 0.9;  // g, scattering anisotropy
  real_t irradiance = 1.0;  // surface irradiance (normalized)
  real_t beam_radius = 0.3;  // beam spot radius (sim units)
  real_t beam_center_x = 0.0;  // beam center x
  real_t beam_center_y = 0.0;  // beam center y
  real_t diffusion = 0.3;  // effective photon diffusion coefficient
  real_t decay = 0.1;  // absorption sink
  // Opsin kinetics (ChR2 two-state model, Fenno et al. 2011)
  real_t opsin_activation_rate = 0.5;  // light-driven closed-to-open rate
  real_t opsin_deactivation_rate = 0.1;  // thermal dark-state recovery
  real_t opsin_saturation = 0.8;  // max open fraction
  // Phototoxicity (Mardinly et al. 2018)
  real_t phototoxicity_threshold = 0.5;  // fluence threshold for ROS
  real_t phototoxicity_rate = 0.002;  // ROS generation rate above threshold
  real_t thermal_coupling = 0.001;  // absorbed light to temperature rise
};

}  // namespace skibidy
}  // namespace bdm

#endif  // PHOTON_PARAMS_H_
