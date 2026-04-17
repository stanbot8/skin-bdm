#ifndef GLUCOSE_POST_HOOK_H_
#define GLUCOSE_POST_HOOK_H_

#include "core/hook_api.h"

namespace bdm {
namespace skibidy {

// Post hook for AGE (advanced glycation end-product) formation.
// Non-enzymatic glycation: glucose reacts with lysine/arginine residues
// on ECM proteins via Amadori rearrangement to form irreversible AGEs.
// RAGE receptor on macrophages and endothelial cells activates NF-kB.
// Brownlee 2005 (doi:10.2337/diabetes.54.6.1615)
// Chavakis et al. 2004 (doi:10.1016/j.micinf.2004.08.004)
struct AGEPostHook {
  DiffusionGrid* age_grid = nullptr;
  DiffusionGrid* glucose_grid = nullptr;
  DiffusionGrid* infl_grid = nullptr;
  const SimParam* sp_ = nullptr;
  bool active = false;

  inline void Init(const GridRegistry& reg, SignalBoard& sig) {
    sp_ = reg.Params();
    active = sp_->glucose_mod.enabled && sp_->diabetic.mode;
    if (!active) return;
    age_grid = reg.Get(fields::kAGEId);
    glucose_grid = reg.Get(fields::kGlucoseId);
    if (!age_grid || !glucose_grid) { active = false; return; }
    infl_grid = reg.InflammationGrid();
  }

  // Called per epidermal wound voxel. Fully self-contained.
  inline void ApplyEpiWound(const VoxelSnapshot& snap, SignalBoard& sig) {
    // AGE formation from glucose via non-enzymatic Maillard reaction.
    // Rate is proportional to glucose^2 (second-order: glucose + protein
    // lysine/arginine produces Schiff base, then Amadori rearrangement).
    // Linear model underestimates hyperglycemia severity: doubling glucose
    // should quadruple AGE formation rate (Brownlee 2001, 2005).
    real_t gluc = glucose_grid->GetConcentration(snap.idx);
    if (gluc > 1e-10) {
      real_t age_formation = sp_->glucose_mod.age_rate * gluc * gluc;
      age_grid->ChangeConcentrationBy(snap.idx, age_formation);
    }
    // RAGE-mediated NF-kB inflammation
    real_t age_val = age_grid->GetConcentration(snap.idx);
    if (age_val > 1e-10 && infl_grid) {
      real_t gate = snap.WoundGate();
      if (gate > 1e-10) {
        infl_grid->ChangeConcentrationBy(snap.idx,
            sp_->age_rage_inflammation * age_val * gate);
      }
    }
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // GLUCOSE_POST_HOOK_H_
