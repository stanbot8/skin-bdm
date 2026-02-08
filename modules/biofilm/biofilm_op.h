#ifndef BIOFILM_H_
#define BIOFILM_H_

#include <cmath>

#include "core/field_names.h"
#include "core/pde.h"
#include "infra/sim_param.h"
#include "core/operation/operation.h"

namespace bdm {
namespace skibidy {

// ---------------------------------------------------------------------------
// BiofilmGrowthOp -- bacterial biofilm dynamics in wound bed.
// Per step: logistic growth, PAMP-driven inflammation production, and
// one-time seeding at biofilm_seed_delay post-wound.
// Clearance is handled by immune cell behaviors (ClearBiofilm helper).
// ---------------------------------------------------------------------------
struct BiofilmGrowthOp : public StandaloneOperationImpl {
  BDM_OP_HEADER(BiofilmGrowthOp);

  void operator()() override {
    auto* sim = Simulation::GetActive();
    auto* sp = sim->GetParam()->Get<SimParam>();

    if (!sp->biofilm_enabled || !sp->wound_enabled) return;

    auto* scheduler = sim->GetScheduler();
    uint64_t step = scheduler->GetSimulatedSteps();
    uint64_t wound_step = static_cast<uint64_t>(sp->wound_trigger_step);
    if (step <= wound_step) return;
    uint64_t wound_age = step - wound_step;

    auto* rm = sim->GetResourceManager();
    auto* biofilm_grid = rm->GetDiffusionGrid(fields::kBiofilm);
    auto* stratum_grid = rm->GetDiffusionGrid(fields::kStratum);
    GridContext ctx(biofilm_grid, sp);

    // Determine inflammation target grid (split/single pattern)
    DiffusionGrid* infl_grid = nullptr;
    if (sp->split_inflammation_enabled) {
      infl_grid = rm->GetDiffusionGrid(fields::kProInflammatory);
    } else {
      infl_grid = rm->GetDiffusionGrid(fields::kInflammation);
    }

    real_t z_max = sp->volume_z_cornified + ctx.box_len;
    bool seeding_step =
        (wound_age == static_cast<uint64_t>(sp->biofilm_seed_delay) &&
         !seeded_);

    for (size_t idx = 0; idx < ctx.n; idx++) {
      real_t z = ctx.Z(idx);
      if (z < 0 || z > z_max) continue;
      real_t x = ctx.X(idx), y = ctx.Y(idx);
      if (!ctx.InWound(x, y)) continue;

      // Only grows where wound bed is open (stratum < 1.0)
      real_t stratum_val = stratum_grid->GetConcentration(idx);
      if (stratum_val >= 1.0) continue;

      real_t current = biofilm_grid->GetConcentration(idx);

      // Seeding: one-time inoculum at seed_delay
      if (seeding_step) {
        biofilm_grid->ChangeConcentrationBy(idx, sp->biofilm_seed_amount);
        current += sp->biofilm_seed_amount;
      }

      // Logistic growth
      if (current > 1e-10) {
        real_t K = sp->biofilm_carrying_capacity;
        real_t delta = sp->biofilm_growth_rate * current * (1.0 - current / K);
        if (delta > 1e-10) {
          biofilm_grid->ChangeConcentrationBy(idx, delta);
        }

        // PAMP-driven inflammation
        real_t infl_delta = sp->biofilm_inflammation_rate * current;
        if (infl_delta > 1e-10) {
          infl_grid->ChangeConcentrationBy(idx, infl_delta);
        }
      }
    }

    if (seeding_step) {
      seeded_ = true;
      std::cout << "[BiofilmGrowthOp] Seeded biofilm inoculum at step "
                << step << std::endl;
    }
  }

 private:
  bool seeded_ = false;
};

}  // namespace skibidy
}  // namespace bdm

#endif  // BIOFILM_H_
