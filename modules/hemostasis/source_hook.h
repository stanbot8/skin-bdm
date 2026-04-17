#ifndef HEMOSTASIS_SOURCE_HOOK_H_
#define HEMOSTASIS_SOURCE_HOOK_H_

#include "core/hook_api.h"

namespace bdm {
namespace skibidy {

// Source hook for hemostasis module.
// Fibrin clot decay + coupling to downstream signals (TGF-beta release
// from platelet degranulation, fibronectin provisional matrix).
// Clark 1996 "The Molecular and Cellular Biology of Wound Repair"
struct HemostasisSourceHook {
  DiffusionGrid* fb_grid = nullptr;
  DiffusionGrid* tgfb_grid = nullptr;
  DiffusionGrid* fn_grid = nullptr;
  const SimParam* sp_ = nullptr;
  const GridRegistry* reg_ = nullptr;
  bool active = false;

  inline void Init(const GridRegistry& reg, SignalBoard& sig) {
    sp_ = reg.Params();
    reg_ = &reg;
    active = sp_->hemostasis.enabled;
    if (!active) return;
    fb_grid = reg.Get(fields::kFibrinId);
    if (!fb_grid) { active = false; return; }
    if (sp_->fibroblast.enabled)
      tgfb_grid = reg.Get(fields::kTGFBetaId);
    if (sp_->fibronectin.enabled)
      fn_grid = reg.Get(fields::kFibronectinId);
  }

  inline void ApplyDecay(const GridContext::WoundMaskData& mask, real_t dt) {
    real_t cw = reg_->CoarseWeight();

    if (sp_->hemostasis.decay > 0) {
      real_t fb_mu_dt = sp_->hemostasis.decay * dt;
      for (size_t idx : mask.epi_wound) {
        size_t si = reg_->CoarseIndex(idx);
        real_t v = fb_grid->GetConcentration(si);
        if (v > 1e-10)
          fb_grid->ChangeConcentrationBy(si, -v * fb_mu_dt * cw);
      }
    }
    if (reg_->PostWound()) {
      real_t tgfb_rate = sp_->hemostasis.tgfb_coupling;
      real_t fn_rate = sp_->hemostasis.fibronectin_coupling;
      if ((tgfb_grid && tgfb_rate > 0) || (fn_grid && fn_rate > 0)) {
        for (size_t idx : mask.epi_wound) {
          size_t fb_si = reg_->CoarseIndex(idx);
          real_t fb = fb_grid->GetConcentration(fb_si);
          if (fb < 1e-10) continue;
          if (tgfb_grid && tgfb_rate > 0) {
            tgfb_grid->ChangeConcentrationBy(idx, fb * tgfb_rate);
          }
          if (fn_grid && fn_rate > 0) {
            fn_grid->ChangeConcentrationBy(
                reg_->CoarseIndex(idx), fb * fn_rate * cw);
          }
        }
      }
    }
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // HEMOSTASIS_SOURCE_HOOK_H_
