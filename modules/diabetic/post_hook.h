#ifndef DIABETIC_POST_HOOK_H_
#define DIABETIC_POST_HOOK_H_

#include "core/hook_api.h"

namespace bdm {
namespace skibidy {

// Post hook for diabetic baseline inflammation.
// AGE/RAGE axis sustains NF-kB activation in open wound bed.
// Gated by stratum: re-epithelialized tissue seals the source.
// Brem & Tomic-Canic 2007 (doi:10.1172/JCI32169)
struct DiabeticPostHook {
  DiffusionGrid* infl_grid = nullptr;
  DiffusionGrid* glucose_grid = nullptr;
  const SimParam* sp_ = nullptr;
  bool active = false;
  bool do_glucose_age = false;

  inline void Init(const GridRegistry& reg, SignalBoard& sig) {
    sp_ = reg.Params();
    active = sp_->diabetic.mode &&
             sp_->diabetic.baseline_inflammation > 0;
    if (!active) return;
    infl_grid = reg.InflammationGrid();
    if (!infl_grid) { active = false; return; }
    do_glucose_age = sp_->glucose_mod.enabled;
    if (do_glucose_age)
      glucose_grid = reg.Get(fields::kGlucoseId);
  }

  inline void ApplyEpiWound(const VoxelSnapshot& snap, SignalBoard& sig) {
    real_t gate = snap.WoundGate();
    if (gate <= 1e-10) return;
    real_t eff_baseline = sp_->diabetic.baseline_inflammation;
    if (do_glucose_age && glucose_grid) {
      real_t gluc = glucose_grid->GetConcentration(snap.idx);
      eff_baseline += sp_->glucose_mod.age_inflammation * gluc;
    }
    infl_grid->ChangeConcentrationBy(snap.idx, eff_baseline * gate);
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // DIABETIC_POST_HOOK_H_
