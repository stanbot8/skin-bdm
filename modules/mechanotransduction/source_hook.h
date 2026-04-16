#ifndef MECHANOTRANSDUCTION_SOURCE_HOOK_H_
#define MECHANOTRANSDUCTION_SOURCE_HOOK_H_

#include "core/grid_registry.h"
#include "core/voxel_snapshot.h"
#include "core/signal_board.h"

namespace bdm {
namespace skibidy {

// Source hook for mechanotransduction module.
// Tissue stiffness = weighted sum of collagen (dominant) and elastin
// (compliance). Healthy dermis ~0.3; scar tissue > 0.8.
// YAP/TAZ nuclear translocation on stiff substrates drives
// myofibroblast differentiation and scar amplification.
// Tomasek et al. 2002 (doi:10.1038/nrm809)
struct MechanoSourceHook {
  DiffusionGrid* stiffness_grid = nullptr;
  DiffusionGrid* col_grid = nullptr;
  DiffusionGrid* elastin_grid = nullptr;
  const SimParam* sp_ = nullptr;
  bool active = false;

  inline void Init(const GridRegistry& reg, SignalBoard& sig) {
    sp_ = reg.Params();
    active = sp_->mechanotransduction.enabled;
    if (!active) return;
    stiffness_grid = reg.Get(fields::kStiffnessId);
    if (!stiffness_grid) { active = false; return; }
    if (sp_->fibroblast.enabled)
      col_grid = reg.Get(fields::kCollagenId);
    if (sp_->elastin.enabled)
      elastin_grid = reg.Get(fields::kElastinId);
  }

  // Dermal: derive stiffness from collagen + elastin composition
  inline void ApplyDermal(const VoxelSnapshot& snap, SignalBoard& sig) {
    real_t col_val = col_grid ?
        col_grid->GetConcentration(snap.coarse_si) : 0;
    real_t el_val = elastin_grid ?
        elastin_grid->GetConcentration(snap.coarse_si) : 0;
    // Collagen dominates stiffness; elastin provides compliance
    real_t stiff = std::min(static_cast<real_t>(1.0),
        col_val * 0.8 - el_val * 0.2);
    stiff = std::max(real_t{0}, stiff);
    real_t stiff_cur = stiffness_grid->GetConcentration(snap.idx);
    real_t stiff_delta = stiff - stiff_cur;
    if (std::abs(stiff_delta) > 1e-10) {
      stiffness_grid->ChangeConcentrationBy(snap.idx, stiff_delta);
    }
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // MECHANOTRANSDUCTION_SOURCE_HOOK_H_
