#ifndef SCAB_PARAMS_H_
#define SCAB_PARAMS_H_

#include "biodynamo.h"

namespace bdm {
namespace skibidy {

struct ScabParams {

  // Scab (protective wound crust)
  // Eming et al. 2007 (doi:10.1038/sj.jid.5700701)
  bool enabled = false;
  real_t wound_seed = 1.0;  // initial scab density in wound
  real_t decay = 0.0008;  // baseline turnover (shedding)
  real_t mmp_degradation = 0.001;  // MMP-mediated breakdown
  real_t reepith_rate = 0.005;  // degradation proportional to barrier
  real_t moisture_softening = 0.002;  // moisture accelerates scab softening
  real_t evaporation_shield = 0.5;  // scab reduces surface water loss
  real_t migration_penalty = 0.15;  // scab impedes keratinocyte migration
};

}  // namespace skibidy
}  // namespace bdm

#endif  // SCAB_PARAMS_H_
