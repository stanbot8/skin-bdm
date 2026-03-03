#ifndef SENESCENCE_PARAMS_H_
#define SENESCENCE_PARAMS_H_

#include "biodynamo.h"

namespace bdm {
namespace skibidy {

struct SenescenceParams {

  // Senescence module (cellular senescence and SASP)
  // Demaria et al. 2014 (doi:10.1016/j.devcel.2014.11.012)
  bool enabled = true;
  real_t diffusion = 0.0;  // immobile
  real_t decay = 0.00005;  // very slow clearance
  real_t wound_rate = 0.00002;  // basal wound-induced DNA damage (gamma-H2AX)
  real_t infl_rate = 0.00005;  // inflammation-driven accumulation (ROS/NF-kB; Campisi 2007)
  real_t age_rate = 0.0002;  // AGE-driven accumulation (diabetic glycative stress)
  real_t immune_clearance_rate = 0.002;  // NK/macrophage-mediated clearance
                                          // (Sagiv et al. 2016, doi:10.18632/aging.100897)
};

}  // namespace skibidy
}  // namespace bdm

#endif  // SENESCENCE_PARAMS_H_
