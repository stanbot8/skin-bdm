#ifndef ELASTIN_PARAMS_H_
#define ELASTIN_PARAMS_H_

#include "biodynamo.h"

namespace bdm {
namespace skibidy {

struct ElastinParams {

  // Elastin (elastic fiber network)
  bool enabled = false;
  real_t diffusion = 0.0;
  real_t decay = 0.0002;
  real_t basal_density = 0.5;
  real_t papillary_density = 0.3;
  real_t production_rate = 0.0005;
  real_t mmp_degradation = 0.003;
};

}  // namespace skibidy
}  // namespace bdm

#endif  // ELASTIN_PARAMS_H_
