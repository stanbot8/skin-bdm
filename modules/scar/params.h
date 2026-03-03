#ifndef SCAR_PARAMS_H_
#define SCAR_PARAMS_H_

#include "biodynamo.h"

namespace bdm {
namespace skibidy {

struct ScarParams {

  // Inflammation-proportional scarring (Extension 1)
  bool proportional_enabled = false;  // cumulative inflammation integral
  real_t accumulation_rate = 0.001;  // scar += inflammation * rate per step
  real_t collagen_threshold = 0.003;  // local collagen above this -> scar, below -> normal
};

}  // namespace skibidy
}  // namespace bdm

#endif  // SCAR_PARAMS_H_
