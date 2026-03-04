#ifndef GLUCOSE_SOURCE_HOOK_H_
#define GLUCOSE_SOURCE_HOOK_H_

#include "core/grid_registry.h"
#include "core/voxel_snapshot.h"
#include "core/signal_board.h"

namespace bdm {
namespace skibidy {

// Source hook for glucose module.
// Perfusion delivers glucose to dermal voxels proportional to vascular
// density. In diabetic mode, hyperglycemia drives AGE accumulation and
// bacterial metabolism.
// Brownlee 2005 (doi:10.2337/diabetes.54.6.1615)
struct GlucoseSourceHook {
  DiffusionGrid* glucose_grid = nullptr;
  const SimParam* sp_ = nullptr;
  bool active = false;

  inline void Init(const GridRegistry& reg, SignalBoard& sig) {
    sp_ = reg.Params();
    active = sp_->glucose_mod.enabled;
    if (!active) return;
    glucose_grid = reg.Get(fields::kGlucoseId);
    if (!glucose_grid) { active = false; return; }
  }

  // Dermal: perfusion supplies glucose
  inline void ApplyDermal(const VoxelSnapshot& snap, SignalBoard& sig) {
    real_t gluc_cur = glucose_grid->GetConcentration(snap.idx);
    real_t gluc_target = sp_->glucose_mod.basal_conc * snap.vasc;
    if (gluc_cur < gluc_target - 1e-10) {
      real_t gain = std::min(gluc_target - gluc_cur,
                             sp_->glucose_mod.perfusion_supply);
      if (gain > 1e-10) {
        glucose_grid->ChangeConcentrationBy(snap.idx, gain);
      }
    }
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // GLUCOSE_SOURCE_HOOK_H_
