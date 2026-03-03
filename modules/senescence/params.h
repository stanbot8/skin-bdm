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
  real_t wound_rate = 0.0001;  // basal wound-induced accumulation
  real_t infl_rate = 0.0003;  // inflammation-driven accumulation
  real_t age_rate = 0.0005;  // AGE-driven accumulation (diabetic)
};

}  // namespace skibidy
}  // namespace bdm

#endif  // SENESCENCE_PARAMS_H_
