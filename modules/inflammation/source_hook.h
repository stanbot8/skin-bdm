#ifndef INFLAMMATION_SOURCE_HOOK_H_
#define INFLAMMATION_SOURCE_HOOK_H_

#include "core/grid_registry.h"
#include "core/voxel_snapshot.h"
#include "core/signal_board.h"

namespace bdm {
namespace skibidy {

// Source hook for wound inflammation (DAMPs).
// Surface DAMPs: keratinocyte IL-1alpha from open epidermis drives
// early innate immune response. Rate tapers as macrophage phagocytosis
// clears damage signals.
// Chen & Nunez 2010 (doi:10.1126/science.1183530)
struct WoundInflSourceHook {
  DiffusionGrid* infl_grid = nullptr;
  const SimParam* sp_ = nullptr;
  bool active = false;
  real_t wound_infl_rate = 0;

  inline void Init(const GridRegistry& reg, SignalBoard& sig) {
    sp_ = reg.Params();
    uint64_t step = reg.Step();
    uint64_t wound_step = static_cast<uint64_t>(sp_->wound.trigger_step);
    bool post_wound = (step > wound_step);

    wound_infl_rate = sp_->wound.inflammation_source_rate;
    active = post_wound && wound_infl_rate > 0;
    if (!active) return;

    // Temporal taper: DAMPs cleared by macrophage phagocytosis
    real_t taper = sp_->wound.inflammation_source_taper;
    if (taper > 0) {
      uint64_t wound_age = step - wound_step;
      wound_infl_rate *= std::exp(-taper * static_cast<real_t>(wound_age));
      if (wound_infl_rate < 1e-10) { active = false; return; }
    }

    infl_grid = reg.InflammationGrid();
    if (!infl_grid) active = false;
  }

  // Epidermal wound: surface DAMPs (keratinocyte IL-1a from open epidermis)
  inline void ApplyEpiWound(const VoxelSnapshot& snap, SignalBoard& sig) {
    if (snap.stratum >= 1.0) return;
    // Clamp at 0: wound void (-1) has same openness as kBasal (0)
    real_t s_clamped = std::max(static_cast<real_t>(0.0), snap.stratum);
    real_t infl_gain = wound_infl_rate * 0.25 * (1.0 - s_clamped);
    if (infl_gain > 1e-10) {
      infl_grid->ChangeConcentrationBy(snap.idx, infl_gain);
    }
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // INFLAMMATION_SOURCE_HOOK_H_
