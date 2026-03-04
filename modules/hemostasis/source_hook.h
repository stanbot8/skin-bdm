#ifndef HEMOSTASIS_SOURCE_HOOK_H_
#define HEMOSTASIS_SOURCE_HOOK_H_

#include "core/grid_registry.h"
#include "core/voxel_snapshot.h"
#include "core/signal_board.h"

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
  bool active = false;
  bool post_wound_ = false;
  bool coarse_ = false;
  real_t coarse_w_ = 1.0;
  std::vector<size_t> coarse_map_;
  // Separate coarse map for fn_grid (may differ from fb_grid)
  std::vector<size_t> fn_coarse_map_;

  inline void Init(const GridRegistry& reg, SignalBoard& sig) {
    sp_ = reg.Params();
    active = sp_->hemostasis.enabled;
    if (!active) return;
    fb_grid = reg.Get(fields::kFibrinId);
    if (!fb_grid) { active = false; return; }
    if (sp_->fibroblast.enabled)
      tgfb_grid = reg.Get(fields::kTGFBetaId);
    if (sp_->fibronectin.enabled)
      fn_grid = reg.Get(fields::kFibronectinId);

    uint64_t wound_step = static_cast<uint64_t>(sp_->wound.trigger_step);
    post_wound_ = (reg.Step() > wound_step);

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
          coarse_map_[i] = fb_grid->GetBoxIndex(pos);
        }
        if (fn_grid) {
          fn_coarse_map_.resize(ctx.n);
          for (size_t i = 0; i < ctx.n; i++) {
            Real3 pos = {ctx.X(i), ctx.Y(i), ctx.Z(i)};
            fn_coarse_map_[i] = fn_grid->GetBoxIndex(pos);
          }
        }
      }
    }
  }

  size_t CSI(size_t fine) const {
    return coarse_ ? coarse_map_[fine] : fine;
  }

  size_t FnCSI(size_t fine) const {
    return coarse_ ? fn_coarse_map_[fine] : fine;
  }

  // Self-contained iteration: fibrin decay + coupling loops
  inline void ApplyDecay(const GridContext::WoundMaskData& mask, real_t dt) {
    real_t cw = coarse_ ? coarse_w_ : 1.0;

    // Manual decay
    if (sp_->hemostasis.decay > 0) {
      real_t fb_mu_dt = sp_->hemostasis.decay * dt;
      for (size_t idx : mask.epi_wound) {
        size_t si = CSI(idx);
        real_t v = fb_grid->GetConcentration(si);
        if (v > 1e-10)
          fb_grid->ChangeConcentrationBy(si, -v * fb_mu_dt * cw);
      }
    }
    // Fibrin -> TGF-beta and fibronectin coupling
    if (post_wound_) {
      real_t tgfb_rate = sp_->hemostasis.tgfb_coupling;
      real_t fn_rate = sp_->hemostasis.fibronectin_coupling;
      if ((tgfb_grid && tgfb_rate > 0) || (fn_grid && fn_rate > 0)) {
        for (size_t idx : mask.epi_wound) {
          size_t fb_si = CSI(idx);
          real_t fb = fb_grid->GetConcentration(fb_si);
          if (fb < 1e-10) continue;
          // TGF-beta is fine grid: write at fine idx
          if (tgfb_grid && tgfb_rate > 0) {
            tgfb_grid->ChangeConcentrationBy(idx, fb * tgfb_rate);
          }
          // Fibronectin is structural: write at coarse idx, scaled
          if (fn_grid && fn_rate > 0) {
            fn_grid->ChangeConcentrationBy(
                FnCSI(idx), fb * fn_rate * cw);
          }
        }
      }
    }
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // HEMOSTASIS_SOURCE_HOOK_H_
