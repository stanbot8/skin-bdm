#ifndef GRID_REGISTRY_H_
#define GRID_REGISTRY_H_

#include <array>
#include "core/field_names.h"
#include "core/material.h"
#include "infra/sim_param.h"
#include "infra/util.h"

namespace bdm {
namespace skibidy {

// ---------------------------------------------------------------------------
// GridRegistry -- single-point grid pointer cache (Unreal Subsystem pattern).
// Filled once per step. All hooks read from this instead of calling
// rm->GetDiffusionGrid() independently. Eliminates 100+ redundant lookups.
// ---------------------------------------------------------------------------
struct GridRegistry {
  static constexpr int kMaxFields = 48;

  std::array<DiffusionGrid*, kMaxFields> grids_{};
  const SimParam* sp_ = nullptr;
  uint64_t step_ = 0;
  MaterialRegistry mat_reg_;
  bool mat_initialized_ = false;

  void Fill(Simulation* sim) {
    auto* rm = sim->GetResourceManager();
    sp_ = sim->GetParam()->Get<SimParam>();
    step_ = GetGlobalStep(sim);
    grids_.fill(nullptr);
    rm->ForEachContinuum([&](Continuum* cm) {
      auto* dg = dynamic_cast<DiffusionGrid*>(cm);
      if (dg) {
        int id = static_cast<int>(dg->GetContinuumId());
        if (id >= 0 && id < kMaxFields) grids_[id] = dg;
      }
    });
    // Lazy-init material registry (once, not per step)
    if (!mat_initialized_) {
      mat_reg_.RegisterSkinDefaults();
      if (sp_->photon.enabled) {
        mat_reg_.RegisterBrainDefaults();
      }
      if (sp_->ra.enabled) {
        mat_reg_.RegisterJointDefaults();
      }
      mat_initialized_ = true;
    }
  }

  // Field access (returns nullptr if field not registered).
  DiffusionGrid* Get(int field_id) const { return grids_[field_id]; }

  const SimParam* Params() const { return sp_; }
  uint64_t Step() const { return step_; }
  const MaterialRegistry& Materials() const { return mat_reg_; }

  // Resolve split vs unified inflammation grid.
  DiffusionGrid* InflammationGrid() const {
    if (sp_->inflammation.split_inflammation_enabled)
      return grids_[fields::kProInflammatoryId];
    return grids_[fields::kInflammationId];
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // GRID_REGISTRY_H_
