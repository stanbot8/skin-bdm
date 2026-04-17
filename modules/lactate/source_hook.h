#ifndef LACTATE_SOURCE_HOOK_H_
#define LACTATE_SOURCE_HOOK_H_

#include "core/hook_api.h"

namespace bdm {
namespace skibidy {

// Source hook for lactate module.
// Hypoxic wound tissue shifts to anaerobic glycolysis, producing lactate.
// Lactate stabilizes HIF-1alpha, boosting VEGF production (Warburg effect
// in wound healing). Perfusion clears lactate via venous washout.
// Production rate now couples to both O2 deficit AND glucose availability:
//   lactate_rate = production_rate * (1 - O2/threshold) * glucose
// This ensures lactate accumulation requires glucose substrate (anaerobic
// glycolysis consumes glucose), making the diabetic hyperglycemia phenotype
// produce more lactate under equal hypoxia.
// Hunt et al. 2007 (doi:10.1089/wound.2007.0402)
// Trabold et al. 2003 (doi:10.1046/j.1524-475x.2003.11621.x)
struct LactateSourceHook {
  DiffusionGrid* lactate_grid = nullptr;
  DiffusionGrid* vegf_grid = nullptr;
  DiffusionGrid* glucose_grid = nullptr;
  const SimParam* sp_ = nullptr;
  bool active = false;
  bool has_glucose = false;

  inline void Init(const GridRegistry& reg, SignalBoard& sig) {
    sp_ = reg.Params();
    active = sp_->lactate.enabled;
    if (!active) return;
    lactate_grid = reg.Get(fields::kLactateId);
    if (!lactate_grid) { active = false; return; }
    if (sp_->angiogenesis.enabled)
      vegf_grid = reg.Get(fields::kVEGFId);
    if (sp_->glucose_mod.enabled) {
      glucose_grid = reg.Get(fields::kGlucoseId);
      has_glucose = (glucose_grid != nullptr);
    }
  }

  // Dermal: hypoxia-driven production + perfusion clearance
  inline void ApplyDermal(const VoxelSnapshot& snap, SignalBoard& sig) {
    if (!snap.in_wound || !snap.post_wound) return;
    if (snap.o2 < sp_->lactate.o2_threshold) {
      real_t hypoxia_frac = (sp_->lactate.o2_threshold - snap.o2) /
                             sp_->lactate.o2_threshold;
      // Scale by glucose availability (default 1.0 if glucose module disabled)
      real_t gluc_factor = 1.0;
      if (has_glucose) {
        gluc_factor = glucose_grid->GetConcentration(snap.idx);
      }
      real_t prod = sp_->lactate.production_rate * hypoxia_frac * gluc_factor;
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
      real_t hypoxia_frac = (sp_->lactate.o2_threshold - snap.o2) /
                             sp_->lactate.o2_threshold;
      real_t gluc_factor = 1.0;
      if (has_glucose) {
        gluc_factor = glucose_grid->GetConcentration(snap.idx);
      }
      real_t prod = sp_->lactate.production_rate * hypoxia_frac * gluc_factor;
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
