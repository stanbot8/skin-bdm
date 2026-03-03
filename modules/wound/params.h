#ifndef WOUND_PARAMS_H_
#define WOUND_PARAMS_H_

#include "biodynamo.h"

namespace bdm {
namespace skibidy {

struct WoundParams {

  // Wound event
  real_t center_x = 0.0;  // x-center of punch biopsy (origin)
  real_t center_y = 0.0;  // y-center of punch biopsy (origin)
  real_t radius = 12.0;  // radius of circular wound (calibrated: 96 wound voxels)
  int trigger_step = 0;  // TOML: trigger_h; 0h = start wounded
  bool enabled = true;  // master switch
  real_t inward_bias = 0.3;  // 0=pure random, 1=pure center-directed
  real_t inflammation_source_rate = 0.003;  // DAMPs/IL-1a from open wound (per step, scaled by 1-stratum)
  real_t inflammation_source_taper = 0.0;  // exp decay rate for DAMP clearance (per step; 0=no taper)
};

}  // namespace skibidy
}  // namespace bdm

#endif  // WOUND_PARAMS_H_
