#ifndef DERMIS_PARAMS_H_
#define DERMIS_PARAMS_H_

#include "biodynamo.h"

namespace bdm {
namespace skibidy {

struct DermisParams {

  // Dermis (dermal tissue integrity -- ECM structural state)
  bool enabled = true;
  real_t diffusion = 0.001;
  real_t decay = 0.0;
  real_t papillary_density = 1.0;
  real_t reticular_density = 0.8;
  real_t hypodermis_density = 0.5;
  real_t collagen_threshold = 0.05;
  real_t collagen_recovery_rate = 0.002;
  real_t mmp_degradation = 0.004;
  real_t papillary_rate_factor = 1.5;
  real_t reticular_rate_factor = 1.0;
  real_t hypodermis_rate_factor = 0.5;
};

}  // namespace skibidy
}  // namespace bdm

#endif  // DERMIS_PARAMS_H_
