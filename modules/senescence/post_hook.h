#ifndef SENESCENCE_POST_HOOK_H_
#define SENESCENCE_POST_HOOK_H_

#include "core/grid_registry.h"
#include "core/voxel_snapshot.h"
#include "core/signal_board.h"

namespace bdm {
namespace skibidy {

// Post hook for senescence module.
// Senescent cells accumulate from DNA damage, inflammation, and AGE glycation.
// SASP output: pro-inflammatory cytokines, MMPs, TGF-beta.
// Demaria et al. 2014 (doi:10.1016/j.devcel.2014.11.012)
// Coppe et al. 2008 (doi:10.1371/journal.pbio.0060301)
struct SenescencePostHook {
  DiffusionGrid* sen_grid = nullptr;
  DiffusionGrid* infl_grid = nullptr;
  DiffusionGrid* age_grid = nullptr;
  DiffusionGrid* prommp_grid = nullptr;
  DiffusionGrid* tgfb_grid = nullptr;
  const SimParam* sp_ = nullptr;
  bool active = false;
  bool do_age = false;

  inline void Init(const GridRegistry& reg, SignalBoard& sig) {
    sp_ = reg.Params();
    active = sp_->senescence.enabled;
    if (!active) return;
    sen_grid = reg.Get(fields::kSenescenceId);
    if (!sen_grid) { active = false; return; }
    infl_grid = reg.InflammationGrid();
    do_age = sp_->glucose_mod.enabled && sp_->diabetic.mode;
    if (do_age)
      age_grid = reg.Get(fields::kAGEId);
    if (sp_->mmp.enabled)
      prommp_grid = reg.Get(fields::kProMMPId);
    if (sp_->fibroblast.enabled)
      tgfb_grid = reg.Get(fields::kTGFBetaId);
  }

  // Called per epidermal wound voxel.
  inline void ApplyEpiWound(const VoxelSnapshot& snap, SignalBoard& sig) {
    size_t sen_si = snap.coarse_si;

    // Accumulation: wound-induced DNA damage
    real_t wound_gate = std::max(static_cast<real_t>(0),
                                 static_cast<real_t>(1) - snap.stratum);
    real_t accum = sp_->senescence.wound_rate * wound_gate;

    // Accumulation: inflammation-driven (ROS, NF-kB)
    if (infl_grid) {
      real_t infl_val = infl_grid->GetConcentration(snap.idx);
      if (infl_val > 1e-10) {
        accum += sp_->senescence.infl_rate * infl_val;
      }
    }

    // Accumulation: AGE-driven (diabetic glycative stress)
    if (do_age && age_grid) {
      real_t age_val = age_grid->GetConcentration(snap.idx);
      if (age_val > 1e-10) {
        accum += sp_->senescence.age_rate * age_val;
      }
    }

    // Diabetic mode: accelerated senescence
    if (sp_->diabetic.mode) {
      accum *= sp_->diabetic.senescence_factor;
    }

    // Senolytic clearance (treatment pathway)
    real_t sen_val = sen_grid->GetConcentration(sen_si);
    if (sp_->senolytic_clearance_rate > 0 && sen_val > 1e-10) {
      accum -= sp_->senolytic_clearance_rate * sen_val;
    }

    if (std::abs(accum) > 1e-10) {
      real_t new_val = sen_val + accum;
      if (new_val < 0) accum = -sen_val;
      sen_grid->ChangeConcentrationBy(sen_si, accum * snap.coarse_w);
    }

    // SASP output (only when senescent cells are present)
    sen_val = sen_grid->GetConcentration(sen_si);
    if (sen_val > 1e-10) {
      // Pro-inflammatory cytokines (IL-6, IL-8)
      if (infl_grid && sp_->sasp_inflammation_rate > 0) {
        infl_grid->ChangeConcentrationBy(snap.idx,
            sp_->sasp_inflammation_rate * sen_val);
      }
      // MMP-3/MMP-9 (as pro-MMP zymogen)
      if (prommp_grid && sp_->sasp_mmp_rate > 0) {
        prommp_grid->ChangeConcentrationBy(snap.idx,
            sp_->sasp_mmp_rate * sen_val);
      }
      // TGF-beta1 (fibrotic)
      if (tgfb_grid && sp_->sasp_tgfb_rate > 0) {
        tgfb_grid->ChangeConcentrationBy(snap.idx,
            sp_->sasp_tgfb_rate * sen_val);
      }
    }
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // SENESCENCE_POST_HOOK_H_
