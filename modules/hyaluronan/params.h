#ifndef HYALURONAN_PARAMS_H_
#define HYALURONAN_PARAMS_H_

#include "biodynamo.h"

namespace bdm {
namespace skibidy {

struct HyaluronanParams {

  // Hyaluronan (HA / GAG ground substance)
  bool enabled = false;
  real_t diffusion = 0.01;
  real_t decay = 0.005;
  real_t basal_density = 0.8;
  real_t reticular_density = 0.4;
  real_t production_rate = 0.002;
  real_t water_retention_factor = 0.5;
  real_t migration_scaffold_factor = 0.3;
};

}  // namespace skibidy
}  // namespace bdm

#endif  // HYALURONAN_PARAMS_H_
