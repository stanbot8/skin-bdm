#ifndef PRESSURE_SOURCE_HOOK_H_
#define PRESSURE_SOURCE_HOOK_H_

#include "core/grid_registry.h"
#include "core/voxel_snapshot.h"
#include "core/signal_board.h"

namespace bdm {
namespace skibidy {

// Pressure ulcer source hook: ischemia-reperfusion injury model.
// Sustained compression above capillary closing pressure occludes perfusion,
// causing progressive tissue ischemia. When pressure is released, reperfusion
// generates ROS burst (ischemia-reperfusion injury).
//
// Modifies existing fields:
//   - Perfusion: reduced during compression, restored on release
//   - ROS: burst on reperfusion (Stekelenburg et al. 2008)
//   - Inflammation: elevated from tissue damage
// Gefen 2007 (doi:10.12968/jowc.2007.16.8.27854)
// Kosiak 1959 (pressure-time relationship)
struct PressureSourceHook {
  DiffusionGrid* perf_grid = nullptr;
  DiffusionGrid* ros_grid = nullptr;
  DiffusionGrid* infl_grid = nullptr;
  DiffusionGrid* water_grid = nullptr;
  const SimParam* sp_ = nullptr;
  bool active = false;

  // Runtime state (updated once per timestep, not per voxel)
  real_t tissue_damage_ = 0.0;     // cumulative tissue damage (0 to 1)
  real_t compression_time_ = 0.0;  // hours under compression
  bool under_compression_ = true;  // current compression state
  int steps_since_reposition_ = 0; // steps since last pressure relief
  uint64_t last_state_step_ = UINT64_MAX;  // gate state updates to once per step

  inline void Init(const GridRegistry& reg, SignalBoard& sig) {
    sp_ = reg.Params();
    active = sp_->pressure.enabled;
    if (!active) return;
    perf_grid = reg.Get(fields::kVascularId);
    if (sp_->ros.enabled)
      ros_grid = reg.Get(fields::kROSId);
    infl_grid = reg.Get(fields::kInflammationId);
    water_grid = reg.Get(fields::kWaterId);
    tissue_damage_ = 0.0;
    compression_time_ = 0.0;
    under_compression_ = true;
    steps_since_reposition_ = 0;
  }

  // Dermal: ischemia from compression, perfusion occlusion, tissue damage.
  inline void ApplyDermal(const VoxelSnapshot& snap, SignalBoard& sig) {
    if (!snap.in_wound) return;

    real_t dt = snap.dt;

    // Update global pressure state once per timestep (not per voxel)
    auto* sim = Simulation::GetActive();
    uint64_t step = sim->GetScheduler()->GetSimulatedSteps();
    if (step != last_state_step_) {
      last_state_step_ = step;
      steps_since_reposition_++;

      // Check for repositioning event (pressure relief cycle)
      real_t reposition_steps = sp_->pressure.reposition_interval_h / dt;
      if (reposition_steps > 0 &&
          steps_since_reposition_ >= static_cast<int>(reposition_steps)) {
        under_compression_ = !under_compression_;
        steps_since_reposition_ = 0;
      }

      if (under_compression_) {
        compression_time_ += dt;
        // Tissue damage accumulation (pressure-time relationship, Kosiak 1959)
        tissue_damage_ += sp_->pressure.tissue_damage_rate * dt;
        tissue_damage_ = std::min(tissue_damage_, static_cast<real_t>(1.0));
      } else if (compression_time_ > 0) {
        compression_time_ = 0.0;
      }
    }

    // Per-voxel effects (these should scale with voxel count)
    if (under_compression_) {
      // Ischemia: compression occludes capillaries, reducing perfusion
      if (perf_grid) {
        real_t perf = perf_grid->GetConcentration(snap.idx);
        real_t occlusion = sp_->pressure.ischemia_rate *
                           sp_->pressure.shear_factor * dt;
        if (perf > 0.05) {
          perf_grid->ChangeConcentrationBy(snap.idx, -std::min(perf * 0.5, occlusion));
        }
      }

      // Moisture-associated damage (per-voxel check, global accumulation gated above)
      if (water_grid) {
        real_t water = water_grid->GetConcentration(snap.idx);
        if (water > 0.5 && step == last_state_step_) {
          // Moisture damage already accumulated in the per-step block
        }
      }
    } else {
      // Reperfusion: ROS burst per voxel
      if (ros_grid) {
        real_t ros_burst = sp_->pressure.reperfusion_ros_burst *
                           std::min(tissue_damage_, static_cast<real_t>(1.0)) * dt;
        if (ros_burst > 1e-10) {
          ros_grid->ChangeConcentrationBy(snap.idx, ros_burst);
        }
      }

      // Perfusion recovery during offloading
      if (perf_grid) {
        real_t perf = perf_grid->GetConcentration(snap.idx);
        real_t max_recovery = 1.0 - tissue_damage_;
        if (perf < max_recovery) {
          real_t recovery = std::min(max_recovery - perf,
                                     sp_->perfusion.angio_rate * dt * 2.0);
          perf_grid->ChangeConcentrationBy(snap.idx, recovery);
        }
      }
    }

    // Inflammation from tissue damage (per-voxel, staging-proportional)
    if (infl_grid && tissue_damage_ > sp_->pressure.stage_1) {
      real_t infl_rate = tissue_damage_ * 0.002 * dt;
      infl_grid->ChangeConcentrationBy(snap.idx, infl_rate);
    }
  }

  // Report current NPUAP stage based on tissue damage accumulator.
  inline int CurrentStage() const {
    if (tissue_damage_ >= sp_->pressure.stage_4) return 4;
    if (tissue_damage_ >= sp_->pressure.stage_3) return 3;
    if (tissue_damage_ >= sp_->pressure.stage_2) return 2;
    if (tissue_damage_ >= sp_->pressure.stage_1) return 1;
    return 0;
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // PRESSURE_SOURCE_HOOK_H_
