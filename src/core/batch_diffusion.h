#ifndef BATCH_DIFFUSION_H_
#define BATCH_DIFFUSION_H_

// BatchDiffusionOp: replaces BDM's ContinuumOp with a leaner scheduler
// that pre-categorizes grids by sub-cycle rate and uses integer step
// counters instead of floating-point time accumulation.
//
// Eliminates per-step overhead of iterating all 29 grids via
// ForEachContinuum, IntegrateTimeAsynchronously bookkeeping, and
// dynamic_cast for gradient checks.

#include <vector>

#include "core/diffusion/diffusion_grid.h"
#include "core/pde.h"
#include "core/operation/operation.h"
#include "core/operation/operation_registry.h"
#include "core/resource_manager.h"
#include "core/simulation.h"
#include "infra/sim_param.h"

namespace bdm {
namespace skibidy {

struct BatchDiffusionOp : public StandaloneOperationImpl {
  BDM_OP_HEADER(BatchDiffusionOp);

  void operator()() override {
    // Lazy init on first call (grids are registered and SetTimeStep
    // values are already configured by RegisterFields).
    if (!initialized_) {
      CategorizeGrids();
      initialized_ = true;
    }

    step_counter_++;

    auto* sim = Simulation::GetActive();
    auto* sp = sim->GetParam()->Get<SimParam>();
    PerfTimer timer(sp->debug_perf);
    real_t base_dt = sim->GetParam()->simulation_time_step;

    // Fast grids: step every tick
    for (auto* g : fast_) {
      g->Step(base_dt);
    }

    // Medium grids: step every subcycle_medium ticks
    if (sp->subcycle_medium > 1) {
      if (step_counter_ % static_cast<uint64_t>(sp->subcycle_medium) == 0) {
        real_t dt = sp->subcycle_medium * base_dt;
        for (auto* g : medium_) {
          g->Step(dt);
        }
      }
    } else {
      for (auto* g : medium_) {
        g->Step(base_dt);
      }
    }

    // Slow grids: step every subcycle_slow ticks
    if (sp->subcycle_slow > 1) {
      if (step_counter_ % static_cast<uint64_t>(sp->subcycle_slow) == 0) {
        real_t dt = sp->subcycle_slow * base_dt;
        for (auto* g : slow_) {
          g->Step(dt);
        }
      }
    } else {
      for (auto* g : slow_) {
        g->Step(base_dt);
      }
    }

    timer.Print("batch_diffusion");
  }

 private:
  bool initialized_ = false;
  uint64_t step_counter_ = 0;
  std::vector<DiffusionGrid*> fast_;
  std::vector<DiffusionGrid*> medium_;
  std::vector<DiffusionGrid*> slow_;

  // Categorize all registered DiffusionGrids by their configured time step.
  // Must be called while BDM's ContinuumOp still exists in the scheduler
  // (needed for GetTimeStep() on grids without explicit SetTimeStep).
  void CategorizeGrids() {
    auto* sim = Simulation::GetActive();
    auto* rm = sim->GetResourceManager();
    auto* sp = sim->GetParam()->Get<SimParam>();
    real_t base_dt = sim->GetParam()->simulation_time_step;

    real_t slow_dt = sp->subcycle_slow * base_dt;
    real_t med_dt = sp->subcycle_medium * base_dt;

    rm->ForEachDiffusionGrid([&](DiffusionGrid* g) {
      // Fixed substance: D=0, decay=0. No physics to solve.
      if (g->IsFixedSubstance()) return;

      // Grids with very large time step are disabled (MarkPrescribed
      // or explicit SetTimeStep(1e30) for fields evolved by custom ops).
      real_t ts = g->GetTimeStep();
      if (ts >= 1e20) return;

      // Categorize by time step value
      if (sp->subcycle_slow > 1 && std::abs(ts - slow_dt) < 1e-6) {
        slow_.push_back(g);
      } else if (sp->subcycle_medium > 1 && std::abs(ts - med_dt) < 1e-6) {
        medium_.push_back(g);
      } else {
        fast_.push_back(g);
      }
    });
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // BATCH_DIFFUSION_H_
