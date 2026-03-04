#ifndef TEMPERATURE_SOURCE_HOOK_H_
#define TEMPERATURE_SOURCE_HOOK_H_

#include "core/grid_registry.h"
#include "core/voxel_snapshot.h"
#include "core/signal_board.h"

namespace bdm {
namespace skibidy {

// Source hook for temperature module.
// Perfusion warms wound tissue toward body temperature (dermal);
// exposed wound surface cools via evaporation (epidermal).
// Wound temperature affects enzyme kinetics (Q10), cell migration,
// and bacterial growth rate.
// Rumney et al. 2015 (doi:10.1111/iwj.12319)
struct TemperatureSourceHook {
  DiffusionGrid* temp_grid = nullptr;
  DiffusionGrid* vasc_grid = nullptr;
  const SimParam* sp_ = nullptr;
  bool active = false;

  inline void Init(const GridRegistry& reg, SignalBoard& sig) {
    sp_ = reg.Params();
    active = sp_->temperature.enabled;
    if (!active) return;
    temp_grid = reg.Get(fields::kTemperatureId);
    if (!temp_grid) { active = false; return; }
    vasc_grid = reg.Get(fields::kVascularId);
  }

  // Dermal: perfusion warms wound tissue toward body temp
  inline void ApplyDermal(const VoxelSnapshot& snap, SignalBoard& sig) {
    if (!snap.in_wound) return;
    real_t temp_cur = temp_grid->GetConcentration(snap.idx);
    if (temp_cur < 1.0 - 1e-10) {
      real_t warm = sp_->temperature.perfusion_warming * snap.vasc;
      real_t gain = std::min(1.0 - temp_cur, warm);
      if (gain > 1e-10) {
        temp_grid->ChangeConcentrationBy(snap.idx, gain);
      }
    }
  }

  // Epidermal wound: surface cooling + perfusion warming from below
  inline void ApplyEpiWound(const VoxelSnapshot& snap, SignalBoard& sig) {
    real_t temp_cur = temp_grid->GetConcentration(snap.idx);
    real_t temp_target = sp_->temperature.wound_surface;
    // Cooling proportional to exposed surface (1 - barrier)
    real_t expose = 1.0 - snap.barrier;
    if (expose > 1e-10 && temp_cur > temp_target + 1e-10) {
      real_t cool = sp_->temperature.surface_cooling * expose;
      real_t loss = std::min(temp_cur - temp_target, cool);
      temp_grid->ChangeConcentrationBy(snap.idx, -loss);
    }
    // Perfusion warming from below
    Real3 below_t = {snap.x, snap.y, -1.0};
    real_t perf_below = vasc_grid->GetValue(below_t);
    if (temp_cur < 1.0 - 1e-10 && perf_below > 0.1) {
      real_t warm = sp_->temperature.perfusion_warming * perf_below;
      real_t gain = std::min(1.0 - temp_cur, warm);
      if (gain > 1e-10) {
        temp_grid->ChangeConcentrationBy(snap.idx, gain);
      }
    }
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // TEMPERATURE_SOURCE_HOOK_H_
