#ifndef NEUROPATHY_PARAMS_H_
#define NEUROPATHY_PARAMS_H_

#include "biodynamo.h"

namespace bdm {
namespace skibidy {

struct NeuropathyParams {

  // Neuropathy module (nerve density and neuropeptide signaling)
  // Suvas 2017 (doi:10.4049/jimmunol.1601751)
  bool enabled = true;
  real_t diffusion = 0.01;  // slow neurite extension
  real_t decay = 0.0;  // structural (no spontaneous loss)
  real_t basal_density = 1.0;  // normalized healthy innervation
};

}  // namespace skibidy
}  // namespace bdm

#endif  // NEUROPATHY_PARAMS_H_
