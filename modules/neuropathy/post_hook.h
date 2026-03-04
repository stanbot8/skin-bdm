#ifndef NEUROPATHY_POST_HOOK_H_
#define NEUROPATHY_POST_HOOK_H_

#include "core/grid_registry.h"
#include "core/voxel_snapshot.h"
#include "core/signal_board.h"

namespace bdm {
namespace skibidy {

// Post hook for neuropathy module (neuropeptide effects).
// Sensory nerve terminals release substance P and CGRP. Local nerve
// density determines neuropeptide flux: denervated wound = reduced
// neurogenic signaling = impaired healing.
// Suvas 2017 (doi:10.4049/jimmunol.1601751)
struct NeuropathyPostHook {
  DiffusionGrid* nerve_grid = nullptr;
  DiffusionGrid* infl_grid = nullptr;
  DiffusionGrid* vasc_grid = nullptr;
  const SimParam* sp_ = nullptr;
  bool active = false;

  inline void Init(const GridRegistry& reg, SignalBoard& sig) {
    sp_ = reg.Params();
    active = sp_->neuropathy.enabled;
    if (!active) return;
    nerve_grid = reg.Get(fields::kNerveId);
    if (!nerve_grid) { active = false; return; }
    infl_grid = reg.InflammationGrid();
    vasc_grid = reg.Get(fields::kVascularId);
  }

  // Called per epidermal wound voxel.
  inline void ApplyEpiWound(const VoxelSnapshot& snap, SignalBoard& sig) {
    size_t nerve_si = snap.coarse_si;
    real_t nerve_val = nerve_grid->GetConcentration(nerve_si);
    if (nerve_val <= 1e-10) return;

    // Substance P: neurogenic inflammation (mast cell degranulation)
    if (infl_grid && sp_->substance_p_inflammation > 0) {
      infl_grid->ChangeConcentrationBy(snap.idx,
          sp_->substance_p_inflammation * nerve_val);
    }

    // CGRP: vasodilation boosts local perfusion
    if (vasc_grid && sp_->cgrp_vasodilation > 0) {
      real_t vasc_val = vasc_grid->GetConcentration(snap.idx);
      real_t vasc_target = sp_->perfusion.basal;
      if (vasc_val < vasc_target) {
        real_t boost = sp_->cgrp_vasodilation * nerve_val;
        real_t headroom = vasc_target - vasc_val;
        real_t gain = std::min(headroom, boost);
        if (gain > 1e-10) {
          vasc_grid->ChangeConcentrationBy(snap.idx, gain);
        }
      }
    }
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // NEUROPATHY_POST_HOOK_H_
