#ifndef ELASTIN_SOURCE_HOOK_H_
#define ELASTIN_SOURCE_HOOK_H_

#include "core/hook_api.h"

namespace bdm {
namespace skibidy {

// Source hook for elastin manual decay (D=0 prescribed field).
// Elastin is a structural ECM grid with manual first-order decay.
struct ElastinDecayHook {
  DiffusionGrid* el_grid = nullptr;
  const SimParam* sp_ = nullptr;
  const GridRegistry* reg_ = nullptr;
  bool active = false;
  inline void Init(const GridRegistry& reg, SignalBoard& sig) {
    sp_ = reg.Params();
    reg_ = &reg;
    active = sp_->elastin.enabled && sp_->elastin.decay > 0;
    if (!active) return;
    el_grid = reg.Get(fields::kElastinId);
    if (!el_grid) { active = false; return; }
  }

  // Self-contained iteration: elastin decay
  inline void ApplyDecay(const GridContext::WoundMaskData& mask, real_t dt) {
    real_t el_mu_dt = sp_->elastin.decay * dt;
    real_t cw = reg_->CoarseWeight();
    for (size_t idx : mask.dermal_all) {
      size_t si = reg_->CoarseIndex(idx);
      real_t v = el_grid->GetConcentration(si);
      if (v > 1e-10)
        el_grid->ChangeConcentrationBy(si, -v * el_mu_dt * cw);
    }
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // ELASTIN_SOURCE_HOOK_H_
