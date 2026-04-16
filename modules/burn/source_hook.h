#ifndef BURN_SOURCE_HOOK_H_
#define BURN_SOURCE_HOOK_H_

#include "core/grid_registry.h"
#include "core/voxel_snapshot.h"
#include "core/signal_board.h"

namespace bdm {
namespace skibidy {

// Burn module source hook: Jackson's three-zone burn wound model (1953).
// Zone of coagulation: irreversible necrosis (dead tissue)
// Zone of stasis: potentially salvageable tissue with impaired perfusion
//   that deteriorates without treatment (converts to coagulation over 24-48h)
// Zone of hyperemia: viable tissue with enhanced inflammatory response
//
// The burn module modifies existing fields:
//   - Perfusion: reduced in stasis zone, enhanced in hyperemia zone
//   - Inflammation: boosted in hyperemia zone
//   - Water: increased evaporative loss (destroyed barrier, TEWL)
//   - Scar: contracture risk accumulation
// Jackson 1953, Herndon 2012 (ISBN 978-1-4377-2786-9)
struct BurnSourceHook {
  DiffusionGrid* perf_grid = nullptr;
  DiffusionGrid* infl_grid = nullptr;
  DiffusionGrid* water_grid = nullptr;
  DiffusionGrid* scar_grid = nullptr;
  const SimParam* sp_ = nullptr;
  bool active = false;

  // Runtime state
  real_t stasis_viability_ = 0.5;  // current stasis zone viability (0 = dead, 1 = healthy)
  uint64_t last_decay_step_ = UINT64_MAX;  // gate: decay stasis once per step

  inline void Init(const GridRegistry& reg, SignalBoard& sig) {
    sp_ = reg.Params();
    active = sp_->burn.enabled;
    if (!active) return;
    perf_grid = reg.Get(fields::kVascularId);
    infl_grid = reg.Get(fields::kInflammationId);
    water_grid = reg.Get(fields::kWaterId);
    if (sp_->scar_enabled)
      scar_grid = reg.Get(fields::kScarId);
    stasis_viability_ = sp_->burn.stasis_viability;
  }

  // Classify voxel into burn zone based on z-depth and depth_fraction.
  // Returns: 0 = coagulation (deepest damage), 1 = stasis, 2 = hyperemia, 3 = unburned
  inline int BurnZone(real_t z, real_t wound_depth) const {
    if (wound_depth <= 0) return 3;
    real_t depth_frac = sp_->burn.depth_fraction;
    real_t norm_z = z / wound_depth;
    if (norm_z < depth_frac * 0.5) return 0;       // coagulation core
    if (norm_z < depth_frac) return 1;              // stasis
    if (norm_z < depth_frac + 0.2) return 2;        // hyperemia
    return 3;                                        // unburned
  }

  // Dermal: zone-specific perfusion and inflammation modifiers.
  inline void ApplyDermal(const VoxelSnapshot& snap, SignalBoard& sig) {
    if (!snap.in_wound) return;

    int zone = BurnZone(snap.z, sp_->volume_z_cornified);

    real_t dt = snap.dt;

    if (zone == 0) {
      // Zone of coagulation: perfusion near zero (dead tissue)
      if (perf_grid) {
        real_t perf = perf_grid->GetConcentration(snap.idx);
        real_t kill = perf * sp_->burn.coagulation_necrosis * 0.01;
        if (kill > 1e-10) {
          perf_grid->ChangeConcentrationBy(snap.idx, -kill);
        }
      }
    } else if (zone == 1) {
      // Zone of stasis: deteriorating perfusion without treatment
      // Decay stasis viability once per timestep (not per voxel)
      auto* sim = Simulation::GetActive();
      uint64_t step = sim->GetScheduler()->GetSimulatedSteps();
      if (step != last_decay_step_) {
        last_decay_step_ = step;
        stasis_viability_ -= sp_->burn.stasis_deterioration * dt;
        stasis_viability_ = std::max(stasis_viability_, static_cast<real_t>(0.0));
      }
      if (perf_grid) {
        real_t perf = perf_grid->GetConcentration(snap.idx);
        // Perfusion scales with remaining viability
        real_t target = perf * stasis_viability_;
        real_t delta = (target - perf) * 0.01;
        if (std::abs(delta) > 1e-10) {
          perf_grid->ChangeConcentrationBy(snap.idx, delta);
        }
      }
    } else if (zone == 2) {
      // Zone of hyperemia: enhanced inflammation
      if (infl_grid) {
        real_t boost = sp_->burn.hyperemia_boost * dt * 0.1;
        infl_grid->ChangeConcentrationBy(snap.idx, boost);
      }
    }

    // Contracture risk: scar accumulation in deep burns
    if (scar_grid && zone <= 1 && sp_->burn.depth_fraction >= 0.7) {
      scar_grid->ChangeConcentrationBy(snap.idx,
                                        sp_->burn.contracture_rate * dt);
    }
  }

  // Epidermal wound: barrier-dependent effects.
  inline void ApplyEpiWound(const VoxelSnapshot& snap, SignalBoard& sig) {
    real_t dt = snap.dt;

    // Transepidermal water loss: destroyed barrier loses water at TEWL multiplier
    if (water_grid) {
      real_t water = water_grid->GetConcentration(snap.idx);
      real_t loss = sp_->burn.fluid_loss_rate * sp_->burn.tewl_multiplier *
                    (1.0 - snap.barrier) * dt;
      if (loss > 0 && water > 0) {
        water_grid->ChangeConcentrationBy(snap.idx, -std::min(water, loss));
      }
    }

    // Eschar formation: dead tissue dries into eschar (reduces inflammation access)
    // Modeled as gradual inflammation suppression in coagulation zone
    int zone = BurnZone(snap.z, sp_->volume_z_cornified);
    if (zone == 0 && infl_grid) {
      real_t infl = infl_grid->GetConcentration(snap.idx);
      real_t eschar_suppress = sp_->burn.eschar_rate * dt;
      if (infl > eschar_suppress) {
        infl_grid->ChangeConcentrationBy(snap.idx, -eschar_suppress);
      }
    }
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // BURN_SOURCE_HOOK_H_
