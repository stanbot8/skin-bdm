#ifndef FIBROBLAST_POST_HOOK_H_
#define FIBROBLAST_POST_HOOK_H_

#include "core/grid_registry.h"
#include "core/voxel_snapshot.h"
#include "core/signal_board.h"

namespace bdm {
namespace skibidy {

// Post hook for fibroblast TGF-beta clearance pathways.
// Decorin-mediated sequestration: small leucine-rich proteoglycan that
// co-deposits with collagen, directly binds active TGF-beta1.
// Yamaguchi et al. 1990 (doi:10.1038/346281a0)
// Tissue receptor clearance: TbRII/TbRI internalization via
// clathrin-mediated endocytosis scales with tissue density.
struct FibroblastPostHook {
  DiffusionGrid* tgfb_grid = nullptr;
  DiffusionGrid* col_grid = nullptr;
  const SimParam* sp_ = nullptr;
  bool active = false;
  bool do_decorin = false;
  bool do_tissue_tgfb = false;

  inline void Init(const GridRegistry& reg, SignalBoard& sig) {
    sp_ = reg.Params();
    do_decorin = sp_->fibroblast.enabled &&
                 sp_->decorin_sequestration_rate > 0;
    do_tissue_tgfb = sp_->fibroblast.enabled &&
                     sp_->fibroblast.tgfb_tissue_clearance > 0;
    active = do_decorin || do_tissue_tgfb;
    if (!active) return;
    tgfb_grid = reg.Get(fields::kTGFBetaId);
    if (!tgfb_grid) { active = false; return; }
    if (do_decorin)
      col_grid = reg.Get(fields::kCollagenId);
  }

  inline void ApplyEpiWound(const VoxelSnapshot& snap, SignalBoard& sig) {
    // Decorin-mediated TGF-beta sequestration
    if (do_decorin && col_grid) {
      real_t col_val = col_grid->GetConcentration(snap.coarse_si);
      if (col_val > 1e-10) {
        real_t tgfb_val = tgfb_grid->GetConcentration(snap.idx);
        if (tgfb_val > 1e-10) {
          real_t sink = sp_->decorin_sequestration_rate *
              col_val * tgfb_val;
          tgfb_grid->ChangeConcentrationBy(snap.idx,
              -std::min(sink, tgfb_val));
        }
      }
    }

    // Tissue receptor-mediated clearance
    if (do_tissue_tgfb) {
      real_t tgfb_val = tgfb_grid->GetConcentration(snap.idx);
      if (tgfb_val > 1e-10) {
        real_t tissue_density = std::max(
            std::max(static_cast<real_t>(0), snap.stratum),
            std::max(static_cast<real_t>(0), snap.vasc));
        tissue_density = std::min(static_cast<real_t>(1), tissue_density);
        if (tissue_density > 1e-10) {
          real_t sink = sp_->fibroblast.tgfb_tissue_clearance *
              tissue_density * tgfb_val;
          tgfb_grid->ChangeConcentrationBy(snap.idx,
              -std::min(sink, tgfb_val));
        }
      }
    }
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // FIBROBLAST_POST_HOOK_H_
