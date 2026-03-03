#ifndef LYMPHATIC_PARAMS_H_
#define LYMPHATIC_PARAMS_H_

#include "biodynamo.h"

namespace bdm {
namespace skibidy {

struct LymphaticParams {

  // Lymphatic system (lymphangiogenesis and interstitial fluid)
  // Kataru et al. 2009 (doi:10.1182/blood-2008-09-176776)
  bool enabled = true;
  real_t diffusion = 0.005;  // slow lymphatic sprouting
  real_t basal_density = 0.8;  // healthy dermis lymphatic density
  real_t regen_rate = 0.0003;  // regeneration rate (slower than vascular)
  real_t vegf_boost = 0.3;  // VEGF promotes lymphangiogenesis
};

}  // namespace skibidy
}  // namespace bdm

#endif  // LYMPHATIC_PARAMS_H_
