#ifndef LACTATE_PARAMS_H_
#define LACTATE_PARAMS_H_

#include "biodynamo.h"

namespace bdm {
namespace skibidy {

struct LactateParams {

  // Lactate (hypoxia metabolite and angiogenesis signal)
  // Hunt et al. 2007 (doi:10.1089/ten.2007.0115)
  bool enabled = true;
  real_t diffusion = 0.01;  // tissue diffusion
  real_t decay = 0.02;  // clearance rate
  real_t production_rate = 0.003;  // hypoxia-driven production
  real_t o2_threshold = 0.3;  // O2 below this triggers production
  real_t vegf_boost = 0.03;  // VEGF production boost factor
  real_t collagen_boost = 0.03;  // collagen synthesis boost factor
  real_t perfusion_clearance = 0.005;  // perfusion-driven lactate clearance
};

}  // namespace skibidy
}  // namespace bdm

#endif  // LACTATE_PARAMS_H_
