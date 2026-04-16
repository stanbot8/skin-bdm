#ifndef FIBRONECTIN_SOURCE_HOOK_H_
#define FIBRONECTIN_SOURCE_HOOK_H_

#include "core/grid_registry.h"
#include "core/voxel_snapshot.h"
#include "core/signal_board.h"

namespace bdm {
namespace skibidy {

// Source hook for fibronectin manual decay (D=0 prescribed field).
// Fibronectin is a structural ECM grid with manual first-order decay.
struct FibronectinDecayHook {
  DiffusionGrid* fn_grid = nullptr;
  const SimParam* sp_ = nullptr;
  const GridRegistry* reg_ = nullptr;
  bool active = false;
  inline void Init(const GridRegistry& reg, SignalBoard& sig) {
    sp_ = reg.Params();
    reg_ = &reg;
    active = sp_->fibronectin.enabled && sp_->fibronectin.decay > 0;
    if (!active) return;
    fn_grid = reg.Get(fields::kFibronectinId);
    if (!fn_grid) { active = false; return; }
  }

  // Self-contained iteration: fibronectin decay
  inline void ApplyDecay(const GridContext::WoundMaskData& mask, real_t dt) {
    real_t fn_mu_dt = sp_->fibronectin.decay * dt;
    real_t cw = reg_->CoarseWeight();
    for (size_t idx : mask.epi_wound) {
      size_t si = reg_->CoarseIndex(idx);
      real_t v = fn_grid->GetConcentration(si);
      if (v > 1e-10)
        fn_grid->ChangeConcentrationBy(si, -v * fn_mu_dt * cw);
    }
    for (size_t idx : mask.dermal_wound) {
      size_t si = reg_->CoarseIndex(idx);
      real_t v = fn_grid->GetConcentration(si);
      if (v > 1e-10)
        fn_grid->ChangeConcentrationBy(si, -v * fn_mu_dt * cw);
    }
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // FIBRONECTIN_SOURCE_HOOK_H_
