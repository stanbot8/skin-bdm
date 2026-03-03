#ifndef INFLAMMATION_PARAMS_H_
#define INFLAMMATION_PARAMS_H_

#include "biodynamo.h"

namespace bdm {
namespace skibidy {

struct InflammationParams {

  // Inflammation field
  real_t diffusion = 0.05;  // lateral cytokine spread
  real_t decay = 0.003;  // background clearance
  real_t migration_threshold = 1.0;  // above this, migration suppressed
  real_t prolif_threshold = 1.5;  // above this, proliferation suppressed

  // Split pro/anti-inflammatory fields (Extension 4)
  bool split_inflammation_enabled = false;  // use separate pro/anti fields
};

}  // namespace skibidy
}  // namespace bdm

#endif  // INFLAMMATION_PARAMS_H_
