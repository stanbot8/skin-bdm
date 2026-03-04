#ifndef LACTATE_SOURCE_HOOK_H_
#define LACTATE_SOURCE_HOOK_H_

#include "core/grid_registry.h"
#include "core/voxel_snapshot.h"
#include "core/signal_board.h"

namespace bdm {
namespace skibidy {

// Source hook for lactate module.
// Hypoxic wound tissue shifts to anaerobic glycolysis, producing lactate.
// Lactate stabilizes HIF-1alpha, boosting VEGF production (Warburg effect
// in wound healing). Perfusion clears lactate via venous washout.
// Hunt et al. 2007 (doi:10.1089/wound.2007.0402)
struct LactateSourceHook {
  DiffusionGrid* lactate_grid = nullptr;
  DiffusionGrid* vegf_grid = nullptr;
  const SimParam* sp_ = nullptr;
  bool active = false;

  inline void Init(const GridRegistry& reg, SignalBoard& sig) {
    sp_ = reg.Params();
    active = sp_->lactate.enabled;
    if (!active) return;
    lactate_grid = reg.Get(fields::kLactateId);
    if (!lactate_grid) { active = false; return; }
    if (sp_->angiogenesis.enabled)
      vegf_grid = reg.Get(fields::kVEGFId);
  }

  // Dermal: hypoxia-driven production + perfusion clearance
  inline void ApplyDermal(const VoxelSnapshot& snap, SignalBoard& sig) {
    if (!snap.in_wound || !snap.post_wound) return;
    if (snap.o2 < sp_->lactate.o2_threshold) {
      real_t prod = sp_->lactate.production_rate *
                    (sp_->lactate.o2_threshold - snap.o2) /
                    sp_->lactate.o2_threshold;
      if (prod > 1e-10) {
        lactate_grid->ChangeConcentrationBy(snap.idx, prod);
      }
    }
    // Perfusion clears lactate (washout)
    real_t lac_cur = lactate_grid->GetConcentration(snap.idx);
    if (lac_cur > 1e-10 && snap.vasc > 0.1) {
      real_t clear = std::min(lac_cur,
                               sp_->lactate.perfusion_clearance * snap.vasc);
      lactate_grid->ChangeConcentrationBy(snap.idx, -clear);
    }
  }

  // Epidermal wound: hypoxia production + lactate-HIF-1a VEGF boost
  inline void ApplyEpiWound(const VoxelSnapshot& snap, SignalBoard& sig) {
    if (!snap.post_wound) return;
    if (snap.o2 < sp_->lactate.o2_threshold) {
      real_t prod = sp_->lactate.production_rate *
                    (sp_->lactate.o2_threshold - snap.o2) /
                    sp_->lactate.o2_threshold;
      if (prod > 1e-10) {
        lactate_grid->ChangeConcentrationBy(snap.idx, prod);
      }
    }

    // Lactate-HIF-1a VEGF boost (merged from ApplyVEGFBoost)
    if (sig.do_vegf && vegf_grid) {
      real_t lac_val = lactate_grid->GetConcentration(snap.idx);
      if (lac_val > 1e-10) {
        if (snap.o2 < sig.vegf_threshold) {
          real_t base_vegf = sig.vegf_prod_rate *
                             (sig.vegf_threshold - snap.o2) / sig.vegf_threshold;
          real_t boost = base_vegf * sp_->lactate.vegf_boost * lac_val;
          if (boost > 1e-10) {
            vegf_grid->ChangeConcentrationBy(snap.idx, boost);
          }
        }
      }
    }
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // LACTATE_SOURCE_HOOK_H_
