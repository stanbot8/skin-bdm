#ifndef FIBRONECTIN_PARAMS_H_
#define FIBRONECTIN_PARAMS_H_

#include "biodynamo.h"

namespace bdm {
namespace skibidy {

struct FibronectinParams {
  real_t haptotaxis_weight = 0.3;  // relative weight of FN haptotaxis (vs O2=1.0)

  // Fibronectin (provisional wound matrix -- migration scaffold)
  // Grinnell 1984 (doi:10.1002/jcb.240260206): fibronectin essential for wound migration
  bool enabled = true;
  real_t decay = 0.005;  // ECM turnover (half-life ~5.8 days)
  real_t deposition_rate = 0.003;  // deposited by activated fibroblasts
  real_t serum_rate = 0.001;  // plasma fibronectin from wound serum
  real_t migration_boost = 0.5;  // max migration speed boost on fibronectin
  real_t wound_seed = 0.3;  // plasma fibronectin co-deposited with fibrin clot (Clark 1996)
  real_t carrying_capacity = 1.0;  // local saturation limit for fibroblast deposition
};

}  // namespace skibidy
}  // namespace bdm

#endif  // FIBRONECTIN_PARAMS_H_
