#ifndef FIBRONECTIN_POST_HOOK_H_
#define FIBRONECTIN_POST_HOOK_H_

#include "core/grid_registry.h"
#include "core/voxel_snapshot.h"
#include "core/signal_board.h"

namespace bdm {
namespace skibidy {

// Post hook for fibronectin serum source.
// Plasma fibronectin seeps into open wound bed from serum, providing
// provisional matrix for cell migration. Gated by stratum coverage
// and subject to carrying capacity.
struct FibronectinPostHook {
  DiffusionGrid* fn_grid = nullptr;
  const SimParam* sp_ = nullptr;
  bool active = false;

  inline void Init(const GridRegistry& reg, SignalBoard& sig) {
    sp_ = reg.Params();
    active = sp_->fibronectin.enabled && sp_->fibronectin.serum_rate > 0;
    if (!active) return;
    fn_grid = reg.Get(fields::kFibronectinId);
    if (!fn_grid) active = false;
  }

  inline void ApplyEpiWound(const VoxelSnapshot& snap, SignalBoard& sig) {
    if (snap.stratum >= 1.0) return;
    size_t fn_si = snap.coarse_si;
    real_t gate = std::max(static_cast<real_t>(0),
                           static_cast<real_t>(1) - snap.stratum);
    // Carrying capacity limits serum FN
    real_t fn_cap = sp_->fibronectin.carrying_capacity;
    if (fn_cap > 0) {
      real_t fn_cur = fn_grid->GetConcentration(fn_si);
      gate *= std::max(static_cast<real_t>(0), 1.0 - fn_cur / fn_cap);
    }
    if (gate > 1e-10) {
      fn_grid->ChangeConcentrationBy(fn_si,
          sp_->fibronectin.serum_rate * gate * snap.coarse_w);
    }
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // FIBRONECTIN_POST_HOOK_H_
