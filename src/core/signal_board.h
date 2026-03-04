#ifndef SIGNAL_BOARD_H_
#define SIGNAL_BOARD_H_

#include "biodynamo.h"

namespace bdm {
namespace skibidy {

// ---------------------------------------------------------------------------
// SignalBoard -- cross-hook coupling channel (Unreal Delegate/Event pattern).
// Step-level signals set during Init(), per-voxel accumulators reset each
// iteration. Replaces ad hoc member access between hook structs.
// ---------------------------------------------------------------------------
struct SignalBoard {
  // --- Step-level signals (set in Init, stable for entire step) ---

  // AngioSourceHook -> LactateSourceHook
  real_t vegf_prod_rate = 0;
  real_t vegf_threshold = 0;
  bool do_vegf = false;
  bool vasc_eligible = false;

  // --- Per-voxel accumulators (reset each iteration) ---

  // MMPPostHook -> matrikine feedback
  real_t total_ecm_degraded = 0;

  void ResetPerVoxel() {
    total_ecm_degraded = 0;
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // SIGNAL_BOARD_H_
