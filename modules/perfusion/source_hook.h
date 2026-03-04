#ifndef PERFUSION_SOURCE_HOOK_H_
#define PERFUSION_SOURCE_HOOK_H_

#include "tissue/calcium.h"
#include "core/grid_registry.h"
#include "core/voxel_snapshot.h"
#include "core/signal_board.h"

namespace bdm {
namespace skibidy {

// Source hook for perfusion-dependent homeostasis.
// O2 pin (Bohr effect modulated), water pin (HA retention),
// water evaporation/recovery (epidermal), calcium gradient restoration.
// These are fundamental tissue physics driven by vascular perfusion.
struct PerfusionSourceHook {
  DiffusionGrid* o2_grid = nullptr;
  DiffusionGrid* water_grid = nullptr;
  DiffusionGrid* ca_grid = nullptr;
  DiffusionGrid* ph_grid = nullptr;
  DiffusionGrid* ha_grid = nullptr;
  const SimParam* sp_ = nullptr;
  bool active = false;
  bool do_ha_water = false;

  inline void Init(const GridRegistry& reg, SignalBoard& sig) {
    sp_ = reg.Params();
    active = true;  // perfusion always active
    o2_grid = reg.Get(fields::kOxygenId);
    water_grid = reg.Get(fields::kWaterId);
    ca_grid = reg.Get(fields::kCalciumId);
    ph_grid = reg.Get(fields::kPHId);
    do_ha_water = sp_->hyaluronan.enabled;
    if (do_ha_water) {
      ha_grid = reg.Get(fields::kHyaluronanId);
    }
  }

  // Dermal: O2 pin + water pin (always runs for all dermal voxels)
  inline void ApplyDermal(const VoxelSnapshot& snap, SignalBoard& sig) {
    // O2 pin: Bohr effect modulated by pH alkalinity
    real_t o2_target = sp_->oxygen_basal_conc * snap.vasc *
                       (1.0 - sp_->ph.bohr_factor * snap.ph);
    real_t o2_delta = o2_target - snap.o2;
    if (std::abs(o2_delta) > 1e-10) {
      o2_grid->ChangeConcentrationBy(snap.idx, o2_delta);
    }

    // Water pin: perfusion-proportional, HA retention
    real_t water_target = sp_->water_basal_conc * snap.vasc;
    if (do_ha_water && ha_grid) {
      real_t local_ha = ha_grid->GetConcentration(snap.coarse_si);
      water_target *= (1.0 + local_ha * sp_->hyaluronan.water_retention_factor);
    }
    real_t water_current = water_grid->GetConcentration(snap.idx);
    if (snap.in_wound) {
      if (water_current < water_target) {
        real_t gain = std::min(water_target - water_current,
                               sp_->water_recovery_rate);
        water_grid->ChangeConcentrationBy(snap.idx, gain);
      }
    } else {
      real_t delta = water_target - water_current;
      if (std::abs(delta) > 1e-10) {
        water_grid->ChangeConcentrationBy(snap.idx, delta);
      }
    }
  }

  // Epidermal wound: water evaporation + recovery + calcium restoration
  inline void ApplyEpiWound(const VoxelSnapshot& snap, SignalBoard& sig) {
    // Water evaporation from exposed surface
    real_t evap = sp_->water_surface_loss_rate * (1.0 - snap.barrier);
    real_t water_current = water_grid->GetConcentration(snap.idx);
    if (evap > 0 && water_current > 0) {
      real_t loss = std::min(water_current, evap);
      water_grid->ChangeConcentrationBy(snap.idx, -loss);
    }
    // Water serum recovery
    real_t water_tgt = GridContext::ExpDecay(snap.z, sp_->water_basal_conc,
                                             sp_->water_decay_length);
    water_current = water_grid->GetConcentration(snap.idx);
    if (water_current < water_tgt) {
      real_t gain = std::min(water_tgt - water_current,
                             sp_->water_recovery_rate);
      water_grid->ChangeConcentrationBy(snap.idx, gain);
    }

    // Calcium restoration (proportional to re-stratification)
    if (snap.stratum >= 0.5) {
      real_t ca_target = CalciumPDE::Profile(sp_, snap.z);
      real_t ca_current = ca_grid->GetConcentration(snap.idx);
      real_t ca_delta = (ca_target - ca_current) *
                         sp_->calcium_recovery_rate * snap.barrier;
      if (std::abs(ca_delta) > 1e-10) {
        ca_grid->ChangeConcentrationBy(snap.idx, ca_delta);
      }
    }
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // PERFUSION_SOURCE_HOOK_H_
