#ifndef NEUROPATHY_SOURCE_HOOK_H_
#define NEUROPATHY_SOURCE_HOOK_H_

#include "core/grid_registry.h"
#include "core/voxel_snapshot.h"
#include "core/signal_board.h"

namespace bdm {
namespace skibidy {

// Source hook for neuropathy module (nerve regeneration).
// Schwann cell-guided neurite extension from wound margins. Rate is
// slow (weeks to months) and modulated by VEGF (neurovascular coupling)
// and inhibited by TGF-beta (scar barrier to neurite growth).
// Boulton et al. 2005 (doi:10.2337/diacare.28.4.956)
struct NeuropathySourceHook {
  DiffusionGrid* nerve_grid = nullptr;
  DiffusionGrid* vegf_grid = nullptr;
  DiffusionGrid* tgfb_grid = nullptr;
  const SimParam* sp_ = nullptr;
  bool active = false;

  inline void Init(const GridRegistry& reg, SignalBoard& sig) {
    sp_ = reg.Params();
    active = sp_->neuropathy.enabled;
    if (!active) return;
    nerve_grid = reg.Get(fields::kNerveId);
    if (!nerve_grid) { active = false; return; }
    if (sp_->angiogenesis.enabled)
      vegf_grid = reg.Get(fields::kVEGFId);
    if (sp_->fibroblast.enabled)
      tgfb_grid = reg.Get(fields::kTGFBetaId);
  }

  // Called per dermal voxel. Uses snap.coarse_si/coarse_w for structural nerve grid.
  inline void ApplyDermal(const VoxelSnapshot& snap, SignalBoard& sig) {
    if (!snap.in_wound || !snap.post_wound) return;

    size_t n_si = snap.coarse_si;
    real_t nerve_val = nerve_grid->GetConcentration(n_si);
    real_t nerve_target = sp_->neuropathy.basal_density;
    if (sp_->diabetic.mode) {
      nerve_target *= sp_->diabetic.nerve_factor;
    }
    if (nerve_val >= nerve_target - 1e-10) return;

    real_t rate = sp_->nerve_regeneration_rate;
    if (sp_->diabetic.mode) {
      rate *= sp_->diabetic.nerve_regen_factor;
    }
    // VEGF neurovascular coupling: promotes nerve sprouting
    if (vegf_grid) {
      real_t local_vegf = vegf_grid->GetConcentration(snap.idx);
      rate *= (1.0 + sp_->nerve_vegf_boost * local_vegf);
    }
    // TGF-beta scar barrier: inhibits neurite extension
    if (tgfb_grid) {
      real_t local_tgfb = tgfb_grid->GetConcentration(snap.idx);
      rate *= std::max(real_t{0},
                       1.0 - sp_->nerve_tgfb_inhibit * local_tgfb);
    }
    real_t gain = std::min(nerve_target - nerve_val, rate);
    if (gain > 1e-10) {
      nerve_grid->ChangeConcentrationBy(n_si, gain * snap.coarse_w);
    }
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // NEUROPATHY_SOURCE_HOOK_H_
