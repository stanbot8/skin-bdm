#ifndef WOUND_H_
#define WOUND_H_

#include "infra/util.h"
#include <cmath>

#include "tissue/basal_division.h"
#include "tissue/differentiation.h"
#include "tissue/migration.h"
#include "tissue/shedding.h"
#include "tissue/keratinocyte.h"
#include "core/composite_field.h"
#include "core/continuum_handoff.h"
#include "core/pde.h"
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

    if (!sp->wound.enabled) return;

    uint64_t step = GetGlobalStep(sim);
    if (fired_ || step != static_cast<uint64_t>(sp->wound.trigger_step)) {
      return;
    }
    fired_ = true;

    real_t cx = sp->wound.center_x;
    real_t cy = sp->wound.center_y;
    real_t r = sp->wound.radius;

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

    // Random angular offset so spawn order varies per run,
    // avoiding deterministic first-mover bias in agent iteration
    auto* random = sim->GetRandom();
    real_t theta_offset = random->Uniform(0, 2.0 * M_PI);

    for (int i = 0; i < n_cells; i++) {
      real_t theta = theta_offset + i * d_theta;
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
// WoundResolution -- dissolves agents when the wound is healed or at end of
// sim. When dissolution_closure_pct > 0, computes wound closure each step
// and dissolves all agents once closure exceeds the threshold. Also fires
// as a safety net at the last step regardless.
// ---------------------------------------------------------------------------
struct WoundResolution : public StandaloneOperationImpl {
  BDM_OP_HEADER(WoundResolution);

  void operator()() override {
    auto* sim = Simulation::GetActive();
    auto* scheduler = sim->GetScheduler();
    auto* sp = sim->GetParam()->Get<SimParam>();

    if (!sp->wound.enabled || resolved_) return;
    PerfTimer timer(sp->debug_perf);

    uint64_t step = GetGlobalStep(sim);
    uint64_t wound_step = static_cast<uint64_t>(sp->wound.trigger_step);
    if (step <= wound_step) return;

    bool is_last = (sp->num_steps > 1 &&
                    step >= static_cast<uint64_t>(sp->num_steps) - 1);

    // Check closure-based dissolution (skip expensive check every step;
    // only check every 50 steps and only after day 10)
    bool dissolve = false;
    if (!is_last && !kerat_dissolved_ && sp->dissolution_closure_pct > 0 &&
        (step - wound_step) > 2400 && step % 50 == 0) {
      dissolve = CheckClosure(sim, sp);
    }

    if (!dissolve && !is_last) return;

    // Before dissolving keratinocytes, stamp any remaining open wound
    // voxels as re-epithelialized. At 90%+ closure the remaining gaps
    // would close within hours; the grid resolution just can't resolve
    // the final convergence.
    auto* rm = sim->GetResourceManager();
    if (dissolve) {
      auto* stratum_grid = rm->GetDiffusionGrid(fields::kStratumId);
      if (stratum_grid) {
        GridContext ctx(stratum_grid, sp);
        real_t r2 = sp->wound.radius * sp->wound.radius;
        real_t z_max = sp->volume_z_cornified + ctx.box_len;
        for (size_t idx = 0; idx < ctx.n; idx++) {
          real_t z = ctx.Z(idx);
          if (z < 0 || z > z_max) continue;
          real_t x = ctx.X(idx), y = ctx.Y(idx);
          real_t dx = x - sp->wound.center_x, dy = y - sp->wound.center_y;
          if (dx * dx + dy * dy > r2) continue;
          real_t sv = stratum_grid->GetConcentration(idx);
          if (sv < 0) {
            stratum_grid->ChangeConcentrationBy(idx, 1.0 - sv);
          }
        }
      }
    }

    // Dissolve agents. Closure-based dissolution only removes
    // keratinocytes; immune cells and fibroblasts have their own
    // lifecycle (apoptosis, emigration). Safety net at last step
    // removes all agents.
    int count = 0;
    rm->ForEachAgent([&](Agent* agent) {
      if (!is_last && !dynamic_cast<Keratinocyte*>(agent)) return;
      ContinuumHandoff(agent);
      agent->RemoveFromSimulation();
      count++;
    });
    if (is_last) resolved_ = true;
    if (dissolve && count > 0) kerat_dissolved_ = true;

    if (sp->debug_wound && count > 0) {
      std::cout << "[wound] dissolved " << count << " agents at step="
                << step << (is_last ? " (safety net)" : " (keratinocyte dissolution)")
                << std::endl;
    }
    timer.Print("wound_resolution");
  }

 private:
  bool resolved_ = false;
  bool kerat_dissolved_ = false;

  bool CheckClosure(Simulation* sim, const SimParam* sp) {
    auto* rm = sim->GetResourceManager();
    auto* stratum_grid = rm->GetDiffusionGrid(fields::kStratumId);
    if (!stratum_grid) return false;

    GridContext ctx(stratum_grid, sp);
    int total = 0, filled = 0;
    real_t r2 = sp->wound.radius * sp->wound.radius;
    real_t z_max = sp->volume_z_cornified + ctx.box_len;

    for (size_t idx = 0; idx < ctx.n; idx++) {
      real_t z = ctx.Z(idx);
      if (z < 0 || z > z_max) continue;
      real_t x = ctx.X(idx), y = ctx.Y(idx);
      real_t dx = x - sp->wound.center_x, dy = y - sp->wound.center_y;
      if (dx * dx + dy * dy > r2) continue;
      total++;
      if (stratum_grid->GetConcentration(idx) > -0.5) filled++;
    }
    if (total == 0) return false;
    real_t pct = 100.0 * static_cast<real_t>(filled) / total;
    return pct >= sp->dissolution_closure_pct;
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // WOUND_H_
