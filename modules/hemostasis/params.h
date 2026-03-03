#ifndef HEMOSTASIS_PARAMS_H_
#define HEMOSTASIS_PARAMS_H_

#include "biodynamo.h"

namespace bdm {
namespace skibidy {

struct HemostasisParams {

  // Hemostasis / fibrin clot (provisional wound matrix)
  bool enabled = false;
  real_t decay = 0.0003;  // slow turnover
  real_t wound_seed = 1.0;  // initial clot density
  real_t tgfb_coupling = 0.0005;  // fibrin -> TGF-beta rate
  real_t fibronectin_coupling = 0.001;  // fibrin -> fibronectin rate
  real_t mmp_degradation = 0.002;  // MMP fibrinolysis rate
  real_t migration_boost = 0.3;  // keratinocyte migration boost
};

}  // namespace skibidy
}  // namespace bdm

#endif  // HEMOSTASIS_PARAMS_H_
