#ifndef WOUND_H_
#define WOUND_H_

#include <cmath>

#include "tissue/basal_division.h"
#include "tissue/differentiation.h"
#include "tissue/migration.h"
#include "tissue/shedding.h"
#include "tissue/keratinocyte.h"
#include "core/composite_field.h"
#include "infra/sim_param.h"
#include "scar/scar_op.h"
#include "core/operation/operation.h"

namespace bdm {
namespace skibidy {

// ---------------------------------------------------------------------------
// WoundEvent -- standalone operation that fires once at a configurable
// timestep, disrupts continuum fields in the wound cylinder, and seeds basal
// keratinocytes at the wound margin for re-epithelialization.
// ---------------------------------------------------------------------------
struct WoundEvent : public StandaloneOperationImpl {
  BDM_OP_HEADER(WoundEvent);

  explicit WoundEvent(CompositeField* fields) : fields_(fields) {}

  void operator()() override {
    auto* sim = Simulation::GetActive();
    auto* scheduler = sim->GetScheduler();
    auto* sp = sim->GetParam()->Get<SimParam>();

    if (!sp->wound_enabled) return;

    uint64_t step = scheduler->GetSimulatedSteps();
    if (fired_ || step != static_cast<uint64_t>(sp->wound_trigger_step)) {
      return;
    }
    fired_ = true;

    real_t cx = sp->wound_center_x;
    real_t cy = sp->wound_center_y;
    real_t r = sp->wound_radius;

    // Delegate field disruption to each PDE channel
    fields_->ApplyWoundAll(sim, cx, cy, r);
    SpawnMarginCells(sim, sp, cx, cy, r);

    if (sp->debug_wound) {
      std::cout << "[wound] punch biopsy step=" << step
                << " center=(" << cx << "," << cy << ") r=" << r << std::endl;
    }
  }

 private:
  CompositeField* fields_;
  bool fired_ = false;

  // Spawn basal keratinocytes evenly around the wound perimeter.
  void SpawnMarginCells(Simulation* sim, const SimParam* sp,
                        real_t cx, real_t cy, real_t r) {
    auto* ctxt = sim->GetExecutionContext();

    real_t d = sp->division_diameter;
    real_t z = d / 2.0;  // on basement membrane

    int n_cells = static_cast<int>(std::round(2.0 * M_PI * r / d));
    if (n_cells < 6) n_cells = 6;

    real_t d_theta = 2.0 * M_PI / n_cells;

    for (int i = 0; i < n_cells; i++) {
      real_t theta = i * d_theta;
      real_t x = cx + r * std::cos(theta);
      real_t y = cy + r * std::sin(theta);

      // Clamp to tissue extent
      x = std::max(sp->tissue_min + d, std::min(x, sp->tissue_max - d));
      y = std::max(sp->tissue_min + d, std::min(y, sp->tissue_max - d));

      auto* cell = new Keratinocyte({x, y, z});
      cell->SetDiameter(d);
      cell->SetStratum(kBasal);
      cell->SetTADivisionsMax(sp->max_ta_divisions);
      cell->InitFromFields(sim, sp);

      cell->AddBehavior(new BasalDivision());
      cell->AddBehavior(new Differentiation());
      cell->AddBehavior(new Migration());
      cell->AddBehavior(new Shedding());
      ctxt->AddAgent(cell);
    }

    if (sp->debug_wound) {
      std::cout << "[wound] spawned " << n_cells << " margin cells" << std::endl;
    }
  }
};

// ---------------------------------------------------------------------------
// WoundResolution -- safety timeout that removes any agents still present
// after the wound has healed. Most cells dissolve individually via per-cell
// UWYN handoff in Differentiation; this catches stragglers (stem cells,
// cells that never stabilized) to guarantee a clean continuum state.
// Resolves when either:
//   - simulation reaches num_steps - 200 (end-of-sim cleanup), OR
//   - >90% of wound-cylinder Stratum voxels are non-zero (coverage)
// ---------------------------------------------------------------------------
struct WoundResolution : public StandaloneOperationImpl {
  BDM_OP_HEADER(WoundResolution);

  void operator()() override {
    auto* sim = Simulation::GetActive();
    auto* scheduler = sim->GetScheduler();
    auto* sp = sim->GetParam()->Get<SimParam>();

    if (!sp->wound_enabled || resolved_) return;

    uint64_t step = scheduler->GetSimulatedSteps();
    if (step <= static_cast<uint64_t>(sp->wound_trigger_step)) return;

    uint64_t wound_age = step -
        static_cast<uint64_t>(sp->wound_trigger_step);

    bool should_resolve = false;
    bool safety_timeout = false;

    // Criterion 1: end-of-sim safety cleanup (fires 200 steps before end)
    uint64_t total = static_cast<uint64_t>(sp->num_steps);
    if (total > 200 && step >= total - 200) {
      should_resolve = true;
      safety_timeout = true;
      if (sp->debug_wound) {
        std::cout << "[wound] safety timeout step=" << step << std::endl;
      }
    }

    // Criterion 2: coverage-based (check every 100 steps).
    // 80% threshold matches the geometric coverage plateau (~83%) --
    // discrete voxel rounding prevents full coverage, so 80% indicates
    // the wound bed is effectively re-epithelialized.
    if (!should_resolve && wound_age > 100 && wound_age % 100 == 0) {
      real_t coverage = ComputeCoverage(sim, sp);
      if (coverage >= 0.80) {
        should_resolve = true;
        if (sp->debug_wound) {
          std::cout << "[wound] coverage " << (coverage * 100)
                    << "% step=" << step << std::endl;
        }
      }
    }

    if (should_resolve) {
      WriteScarStrata(sim);
      if (safety_timeout) {
        // End-of-sim: remove ALL agents (fibroblasts, immune, etc.)
        RemoveAllAgents(sim);
      } else {
        // Coverage-based: only dissolve keratinocytes.
        // Fibroblasts and immune cells continue their natural lifecycle
        // (collagen remodeling continues after re-epithelialization).
        RemoveKeratinocytes(sim);
      }
      resolved_ = true;
      int removed = last_removed_;
      if (sp->debug_wound) {
        std::cout << "[wound] dissolved " << removed
                  << (safety_timeout ? " agents (safety)" : " keratinocytes")
                  << std::endl;
      }
    }
  }

 private:
  bool resolved_ = false;
  int last_removed_ = 0;

  real_t ComputeCoverage(Simulation* sim, const SimParam* sp) {
    auto* rm = sim->GetResourceManager();
    auto* grid = rm->GetDiffusionGrid(fields::kStratum);

    real_t cx = sp->wound_center_x;
    real_t cy = sp->wound_center_y;
    real_t r2 = sp->wound_radius * sp->wound_radius;

    size_t res = grid->GetResolution();
    size_t n = grid->GetNumBoxes();
    auto dims = grid->GetDimensions();
    real_t lo = static_cast<real_t>(dims[0]);
    real_t box_len = grid->GetBoxLength();
    real_t z_max = sp->volume_z_cornified + box_len;

    int total = 0;
    int filled = 0;
    for (size_t idx = 0; idx < n; idx++) {
      uint32_t ix = static_cast<uint32_t>(idx % res);
      uint32_t iy = static_cast<uint32_t>((idx / res) % res);
      uint32_t iz = static_cast<uint32_t>(idx / (res * res));

      real_t z = lo + iz * box_len + box_len / 2.0;
      if (z < 0 || z > z_max) continue;  // dermis or above tissue, skip

      real_t x = lo + ix * box_len + box_len / 2.0;
      real_t y = lo + iy * box_len + box_len / 2.0;

      real_t dx = x - cx;
      real_t dy = y - cy;
      if (dx * dx + dy * dy <= r2) {
        total++;
        if (std::abs(grid->GetConcentration(idx)) > 0.5) {
          filled++;
        }
      }
    }
    if (total == 0) return 1.0;
    return static_cast<real_t>(filled) / static_cast<real_t>(total);
  }

  void WriteScarStrata(Simulation* sim) {
    auto* rm = sim->GetResourceManager();
    rm->ForEachAgent([](Agent* agent) {
      auto* cell = dynamic_cast<Keratinocyte*>(agent);
      if (cell) WriteScarValue(cell);
    });
  }

  void RemoveKeratinocytes(Simulation* sim) {
    auto* rm = sim->GetResourceManager();
    int count = 0;
    rm->ForEachAgent([&](Agent* agent) {
      if (dynamic_cast<Keratinocyte*>(agent)) {
        agent->RemoveFromSimulation();
        count++;
      }
    });
    last_removed_ = count;
  }

  void RemoveAllAgents(Simulation* sim) {
    auto* rm = sim->GetResourceManager();
    int count = 0;
    rm->ForEachAgent([&](Agent* agent) {
      agent->RemoveFromSimulation();
      count++;
    });
    last_removed_ = count;
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // WOUND_H_
