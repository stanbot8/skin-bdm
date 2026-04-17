#ifndef PH_SOURCE_HOOK_H_
#define PH_SOURCE_HOOK_H_

#include "core/hook_api.h"

namespace bdm {
namespace skibidy {

// Source hook for pH module.
// Wound acidification: tissue metabolism produces lactic acid. Recovery
// scales with perfusion (O2 delivery) and tissue coverage (barrier).
// Schneider et al. 2007 (doi:10.1007/s00403-006-0713-x)
struct PHSourceHook {
  DiffusionGrid* ph_grid = nullptr;
  DiffusionGrid* vasc_grid = nullptr;
  const SimParam* sp_ = nullptr;
  bool active = false;

  inline void Init(const GridRegistry& reg, SignalBoard& sig) {
    sp_ = reg.Params();
    ph_grid = reg.Get(fields::kPHId);
    active = (ph_grid != nullptr);
    if (active)
      vasc_grid = reg.Get(fields::kVascularId);
  }

  // Epidermal wound: pH recovery (wound acidification)
  inline void ApplyEpiWound(const VoxelSnapshot& snap, SignalBoard& sig) {
    size_t ph_idx = snap.coarse_si;
    real_t ph_val = ph_grid->GetConcentration(ph_idx);
    if (ph_val <= 1e-6) return;

    Real3 below = {snap.x, snap.y, -1.0};
    real_t local_perf = vasc_grid->GetValue(below);
    real_t eff_rate = sp_->ph.recovery_rate *
        std::max(static_cast<real_t>(0.2), local_perf) * snap.barrier;
    if (sp_->diabetic.mode) {
      eff_rate *= sp_->diabetic.ph_recovery_factor;
    }
    real_t drop = std::min(ph_val, eff_rate);
    ph_grid->ChangeConcentrationBy(ph_idx, -drop * snap.coarse_w);
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // PH_SOURCE_HOOK_H_
