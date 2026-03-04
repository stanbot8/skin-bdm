#ifndef ELASTIN_SOURCE_HOOK_H_
#define ELASTIN_SOURCE_HOOK_H_

#include "core/grid_registry.h"
#include "core/voxel_snapshot.h"
#include "core/signal_board.h"

namespace bdm {
namespace skibidy {

// Source hook for elastin manual decay (D=0 prescribed field).
// Elastin is a structural ECM grid with manual first-order decay.
struct ElastinDecayHook {
  DiffusionGrid* el_grid = nullptr;
  const SimParam* sp_ = nullptr;
  bool active = false;
  bool coarse_ = false;
  real_t coarse_w_ = 1.0;
  std::vector<size_t> coarse_map_;

  inline void Init(const GridRegistry& reg, SignalBoard& sig) {
    sp_ = reg.Params();
    active = sp_->elastin.enabled && sp_->elastin.decay > 0;
    if (!active) return;
    el_grid = reg.Get(fields::kElastinId);
    if (!el_grid) { active = false; return; }

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
          coarse_map_[i] = el_grid->GetBoxIndex(pos);
        }
      }
    }
  }

  size_t CSI(size_t fine) const {
    return coarse_ ? coarse_map_[fine] : fine;
  }

  // Self-contained iteration: elastin decay
  inline void ApplyDecay(const GridContext::WoundMaskData& mask, real_t dt) {
    real_t el_mu_dt = sp_->elastin.decay * dt;
    real_t cw = coarse_ ? coarse_w_ : 1.0;
    for (size_t idx : mask.dermal_all) {
      size_t si = CSI(idx);
      real_t v = el_grid->GetConcentration(si);
      if (v > 1e-10)
        el_grid->ChangeConcentrationBy(si, -v * el_mu_dt * cw);
    }
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // ELASTIN_SOURCE_HOOK_H_
