#ifndef PH_PARAMS_H_
#define PH_PARAMS_H_

#include "biodynamo.h"

namespace bdm {
namespace skibidy {

struct PhParams {

  // pH (wound alkalinity -- acid mantle disruption)
  // Schneider et al. 2007 (doi:10.1111/j.1524-475X.2007.00230.x): acidic wound
  // bed (pH 5.5-6.5) promotes keratinocyte migration and re-epithelialization.
  // Field stores alkalinity excess: 0 = normal skin, 1.0 = fresh wound.
  real_t recovery_rate = 0.003;  // acidification rate per step (perfusion-scaled)
  real_t migration_suppression = 0.5;  // max speed reduction at full alkalinity
  real_t mmp_boost = 0.5;  // MMP activity amplification at full alkalinity
  real_t biofilm_boost = 0.4;  // biofilm growth boost at full alkalinity
  real_t bohr_factor = 0.3;  // O2 delivery reduction at full alkalinity (Bohr effect)
};

}  // namespace skibidy
}  // namespace bdm

#endif  // PH_PARAMS_H_
