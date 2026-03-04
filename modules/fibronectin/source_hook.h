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
  bool active = false;
  bool coarse_ = false;
  real_t coarse_w_ = 1.0;
  std::vector<size_t> coarse_map_;

  inline void Init(const GridRegistry& reg, SignalBoard& sig) {
    sp_ = reg.Params();
    active = sp_->fibronectin.enabled && sp_->fibronectin.decay > 0;
    if (!active) return;
    fn_grid = reg.Get(fields::kFibronectinId);
    if (!fn_grid) { active = false; return; }

    coarse_ = sp_->grid_resolution_structural > 0 &&
              sp_->grid_resolution_structural != sp_->grid_resolution;
    if (coarse_) {
      real_t rf = static_cast<real_t>(sp_->grid_resolution);
      real_t rc = static_cast<real_t>(sp_->grid_resolution_structural);
      coarse_w_ = (rc * rc * rc) / (rf * rf * rf);
      auto* strat = reg.Get(fields::kStratumId);
      if (strat) {
        GridContext ctx(strat, sp_);
        coarse_map_.resize(ctx.n);
        for (size_t i = 0; i < ctx.n; i++) {
          Real3 pos = {ctx.X(i), ctx.Y(i), ctx.Z(i)};
          coarse_map_[i] = fn_grid->GetBoxIndex(pos);
        }
      }
    }
  }

  size_t CSI(size_t fine) const {
    return coarse_ ? coarse_map_[fine] : fine;
  }

  // Self-contained iteration: fibronectin decay
  inline void ApplyDecay(const GridContext::WoundMaskData& mask, real_t dt) {
    real_t fn_mu_dt = sp_->fibronectin.decay * dt;
    real_t cw = coarse_ ? coarse_w_ : 1.0;
    for (size_t idx : mask.epi_wound) {
      size_t si = CSI(idx);
      real_t v = fn_grid->GetConcentration(si);
      if (v > 1e-10)
        fn_grid->ChangeConcentrationBy(si, -v * fn_mu_dt * cw);
    }
    for (size_t idx : mask.dermal_wound) {
      size_t si = CSI(idx);
      real_t v = fn_grid->GetConcentration(si);
      if (v > 1e-10)
        fn_grid->ChangeConcentrationBy(si, -v * fn_mu_dt * cw);
    }
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // FIBRONECTIN_SOURCE_HOOK_H_
