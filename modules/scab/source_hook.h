#ifndef SCAB_SOURCE_HOOK_H_
#define SCAB_SOURCE_HOOK_H_

#include "core/grid_registry.h"
#include "core/voxel_snapshot.h"
#include "core/signal_board.h"

namespace bdm {
namespace skibidy {

// Decay hook for scab module.
// Scab degrades via three mechanisms:
//   1. Baseline turnover (shedding, drying, mechanical loss)
//   2. Re-epithelialization undermining (barrier-proportional: as keratinocytes
//      migrate underneath, scab lifts from below)
//   3. MMP proteolysis (enzymatic degradation of fibrin/collagen scab matrix)
//   4. Moisture softening (wet wound healing dissolves/softens scab)
// Singer & Clark 1999 (doi:10.1056/NEJM199909023411006)
struct ScabDecayHook {
  DiffusionGrid* scab_grid = nullptr;
  DiffusionGrid* stratum_grid = nullptr;
  DiffusionGrid* mmp_grid = nullptr;
  DiffusionGrid* water_grid = nullptr;
  const SimParam* sp_ = nullptr;
  const GridRegistry* reg_ = nullptr;
  bool active = false;
  inline void Init(const GridRegistry& reg, SignalBoard& sig) {
    sp_ = reg.Params();
    reg_ = &reg;
    active = sp_->scab.enabled;
    if (!active) return;
    scab_grid = reg.Get(fields::kScabId);
    if (!scab_grid) { active = false; return; }
    stratum_grid = reg.Get(fields::kStratumId);
    if (sp_->mmp.enabled) mmp_grid = reg.Get(fields::kMMPId);
    water_grid = reg.Get(fields::kWaterId);
  }

  // Self-contained iteration: scab decay over epidermal wound voxels.
  inline void ApplyDecay(const GridContext::WoundMaskData& mask, real_t dt) {
    if (!active) return;
    real_t cw = reg_->CoarseWeight();
    real_t base_decay = sp_->scab.decay * dt;
    real_t reepith = sp_->scab.reepith_rate * dt;
    real_t mmp_deg = sp_->scab.mmp_degradation * dt;
    real_t moist = sp_->scab.moisture_softening * dt;

    for (size_t idx : mask.epi_wound) {
      size_t si = reg_->CoarseIndex(idx);
      real_t v = scab_grid->GetConcentration(si);
      if (v <= 1e-10) continue;

      real_t loss = base_decay;

      // Re-epithelialization undermining: barrier goes 0 (void) to 1 (healed)
      if (stratum_grid) {
        real_t stratum = stratum_grid->GetConcentration(idx);
        real_t barrier = std::clamp((stratum + 1) / 4, real_t{0}, real_t{1});
        loss += reepith * barrier;
      }

      // MMP proteolysis
      if (mmp_grid && mmp_deg > 0) {
        real_t mmp_val = mmp_grid->GetConcentration(idx);
        loss += mmp_deg * mmp_val;
      }

      // Moisture softening
      if (water_grid && moist > 0) {
        real_t water_val = water_grid->GetConcentration(idx);
        loss += moist * water_val;
      }

      real_t delta = std::min(v, v * loss);
      if (delta > 1e-10) {
        scab_grid->ChangeConcentrationBy(si, -delta * cw);
      }
    }
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // SCAB_SOURCE_HOOK_H_
