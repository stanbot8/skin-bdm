#ifndef FIBRONECTIN_SOURCE_HOOK_H_
#define FIBRONECTIN_SOURCE_HOOK_H_

#include "core/hook_api.h"

namespace bdm {
namespace skibidy {

// Source hook for fibronectin manual decay (D=0 prescribed field).
// Fibronectin is a structural ECM grid with manual first-order decay.
struct FibronectinDecayHook {
  DiffusionGrid* fn_grid = nullptr;
  DiffusionGrid* col_grid = nullptr;
  const SimParam* sp_ = nullptr;
  const GridRegistry* reg_ = nullptr;
  bool active = false;
  bool do_replacement = false;

  inline void Init(const GridRegistry& reg, SignalBoard& sig) {
    sp_ = reg.Params();
    reg_ = &reg;
    active = sp_->fibronectin.enabled && sp_->fibronectin.decay > 0;
    if (!active) return;
    fn_grid = reg.Get(fields::kFibronectinId);
    if (!fn_grid) { active = false; return; }
    do_replacement = sp_->fibronectin.collagen_replacement > 0 &&
                     sp_->fibroblast.enabled;
    if (do_replacement) col_grid = reg.Get(fields::kCollagenId);
  }

  // Fibronectin decay with collagen-driven replacement (Clark 1990):
  // as mature collagen accumulates, it physically displaces provisional FN.
  inline void ApplyDecay(const GridContext::WoundMaskData& mask, real_t dt) {
    real_t base_mu_dt = sp_->fibronectin.decay * dt;
    real_t k_replace = sp_->fibronectin.collagen_replacement * dt;
    real_t cw = reg_->CoarseWeight();
    auto decay_voxel = [&](size_t idx) {
      size_t si = reg_->CoarseIndex(idx);
      real_t v = fn_grid->GetConcentration(si);
      if (v <= 1e-10) return;
      real_t eff_mu = base_mu_dt;
      if (do_replacement && col_grid) {
        real_t col = col_grid->GetConcentration(si);
        eff_mu += k_replace * col;
      }
      fn_grid->ChangeConcentrationBy(si, -v * eff_mu * cw);
    };
    for (size_t idx : mask.epi_wound) decay_voxel(idx);
    for (size_t idx : mask.dermal_wound) decay_voxel(idx);
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // FIBRONECTIN_SOURCE_HOOK_H_
