#ifndef TUMOR_BEHAVIOR_H_
#define TUMOR_BEHAVIOR_H_

#include <cmath>

#include "tumor/tumor_cell.h"
#include "core/continuum_handoff.h"
#include "core/field_names.h"
#include "infra/sim_param.h"
#include "infra/util.h"

namespace bdm {
namespace skibidy {

// ---------------------------------------------------------------------------
// Behavior: G1->S->G2->M cell cycle for tumor cells.
// Same state machine as BasalDivision but with tumor-specific overrides:
//   - Phase durations scaled by tumor_cycle_factor (> 1 = slower for BCC)
//   - Contact inhibition threshold = tumor_max_neighbors (higher = less)
//   - Unlimited divisions (no stem/TA distinction)
//   - No stratum gating (divides at any position)
//   - Random division direction (no wound bias)
//   - O2-dependent: hypoxia still suppresses proliferation
// ---------------------------------------------------------------------------
struct TumorBehavior : public Behavior {
  BDM_BEHAVIOR_HEADER(TumorBehavior, Behavior, 1);

  TumorBehavior() { AlwaysCopyToNew(); }
  virtual ~TumorBehavior() {}

  // Cached tumor count: updated once per step, shared across all instances
  static inline uint64_t cached_step_ = UINT64_MAX;
  static inline int cached_count_ = 0;

  void Run(Agent* agent) override {
    auto* cell = dynamic_cast<TumorCell*>(agent);
    if (!cell) return;

    auto* sim = Simulation::GetActive();
    auto* sp = sim->GetParam()->Get<SimParam>();
    real_t dt = sim->GetParam()->simulation_time_step;

    // --- Apoptosis check (any phase) ---
    if (sp->tumor_apoptosis_rate > 0) {
      auto* rng = sim->GetRandom();
      if (rng->Uniform(0, 1) < sp->tumor_apoptosis_rate) {
        ContinuumHandoff(cell);
        cell->RemoveFromSimulation();
        return;
      }
    }

    cell->SetCellAge(cell->GetCellAge() + dt);
    cell->SetPhaseElapsed(cell->GetPhaseElapsed() + dt);
    auto* random = sim->GetRandom();
    real_t ran = random->Uniform(0, 1);

    // Scaled phase durations
    real_t cf = sp->tumor_cycle_factor;

    switch (cell->GetCyclePhase()) {
      case kG1: {
        // BCC cells have an extended G1 phase (quiescent checkpoint),
        // resulting in Ki-67 ~27.4% vs the ~59% from uniform scaling.
        // Khoo et al. 2019 (doi:10.2340/00015555-3325)
        // Constant hazard rate: p = dt/T per step gives exponentially
        // distributed exit times with mean T (reliably adjustable).
        real_t T_g1 = sp->g1_duration * cf * sp->tumor_g1_factor;
        real_t p = dt / T_g1;

        // O2-dependent proliferation: hypoxia suppresses G1->S
        auto* rm = sim->GetResourceManager();
        auto* o2_grid = rm->GetDiffusionGrid(fields::kOxygenId);
        Real3 qpos = ClampToBounds(cell->GetPosition(), sim->GetParam());
        real_t o2 = o2_grid->GetValue(qpos);
        real_t o2_factor = std::min(1.0, o2 / sp->tumor_o2_threshold);
        if (o2_factor < 0) o2_factor = 0;
        p *= o2_factor;

        if (p > ran) {
          cell->SetCyclePhase(kS);
          cell->SetPhaseElapsed(0);
        }
        break;
      }
      case kS: {
        if (cell->GetDiameter() < sp->division_diameter) {
          cell->ChangeVolume(sp->tumor_growth_rate);
        }
        real_t p = dt / (sp->s_duration * cf);
        if (p > ran) {
          cell->SetCyclePhase(kG2);
          cell->SetPhaseElapsed(0);
        }
        break;
      }
      case kG2: {
        real_t p = dt / (sp->g2_duration * cf);
        if (p > ran) {
          cell->SetCyclePhase(kM);
          cell->SetPhaseElapsed(0);
        }
        break;
      }
      case kM: {
        real_t p = dt / (sp->m_duration * cf);
        if (p > ran) {
          bool divided = false;

          // Carrying capacity (cached once per step, not per cell)
          bool at_capacity = false;
          if (sp->tumor_max_cells > 0) {
            auto step = GetGlobalStep(sim);
            if (step != cached_step_) {
              cached_step_ = step;
              cached_count_ = 0;
              sim->GetResourceManager()->ForEachAgent([&](Agent* a) {
                if (dynamic_cast<TumorCell*>(a)) cached_count_++;
              });
            }
            at_capacity = (cached_count_ >= sp->tumor_max_cells);
          }

          if (!at_capacity && cell->GetDiameter() >= sp->division_diameter) {
            // Contact inhibition with tumor-specific threshold
            // Search within 1.5 diameters to catch close but non-touching
            // neighbors in a random (non-ideal) packing.
            int neighbors = 0;
            auto* ctxt = sim->GetExecutionContext();
            auto count = L2F([&](Agent*, real_t) { neighbors++; });
            real_t d = cell->GetDiameter();
            ctxt->ForEachNeighbor(count, *cell, 2.25 * d * d);

            // Soft contact inhibition (YAP/TAZ mechanotransduction):
            // division probability decreases smoothly with neighbor count.
            bool can_divide = false;
            if (sp->tumor_ci_steepness > 0) {
              real_t ratio = static_cast<real_t>(neighbors) /
                             static_cast<real_t>(sp->tumor_max_neighbors);
              if (ratio < 1.0) {
                real_t p = 1.0 - std::pow(ratio, sp->tumor_ci_steepness);
                can_divide = (random->Uniform(0, 1) < p);
              }
            } else {
              can_divide = (neighbors < sp->tumor_max_neighbors);
            }

            if (can_divide) {
              // Random division direction (no wound bias)
              real_t theta = random->Uniform(0, 2.0 * M_PI);
              real_t dx = std::cos(theta);
              real_t dy = std::sin(theta);

              // Vertical extrusion: base component ensures 3D nodular
              // growth even at low density; increases with crowding.
              real_t crowding = static_cast<real_t>(neighbors) /
                                static_cast<real_t>(sp->tumor_max_neighbors);
              real_t dz;
              if (crowding > 0.3) {
                real_t extrusion = (crowding - 0.3) / 0.7;
                dz = 0.6 + 1.4 * extrusion;
              } else {
                dz = random->Uniform(0.3, 0.6);
              }

              cell->Divide({dx, dy, dz});
              divided = true;
            }
          }
          cell->SetCyclePhase(divided ? kG1 : kG0);
          cell->SetPhaseElapsed(0);
        }
        break;
      }
      default: {  // kG0 — quiescent; recheck contact inhibition
        cell->IncrementG0Steps();
        int g0 = cell->GetG0Steps();

        // Recheck neighbors every 10 steps (not every step)
        if (g0 % 10 == 0) {
          int neighbors = 0;
          auto* ctxt = sim->GetExecutionContext();
          auto count_fn = L2F([&](Agent*, real_t) { neighbors++; });
          real_t d = cell->GetDiameter();
          ctxt->ForEachNeighbor(count_fn, *cell, 2.25 * d * d);

          // Soft CI re-entry: probability based on neighbor count,
          // evaluated against 2/3 effective max (hysteresis).
          bool re_enter = false;
          if (sp->tumor_ci_steepness > 0) {
            real_t eff_max = sp->tumor_max_neighbors * 2.0 / 3.0;
            real_t ratio = static_cast<real_t>(neighbors) / eff_max;
            if (ratio < 1.0) {
              real_t p = 1.0 - std::pow(ratio, sp->tumor_ci_steepness);
              re_enter = (random->Uniform(0, 1) < p);
            }
          } else {
            re_enter = (neighbors < sp->tumor_max_neighbors * 2 / 3);
          }

          if (re_enter) {
            cell->SetG0Steps(0);
            cell->SetCyclePhase(kG1);
            cell->SetPhaseElapsed(0);
            break;
          }
        }

        // UWYN handoff: quiescent long enough -> convert to field
        if (sp->tumor_handoff_delay > 0 && g0 > sp->tumor_handoff_delay) {
          ContinuumHandoff(cell);
          cell->RemoveFromSimulation();
          return;
        }
        break;
      }
    }
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // TUMOR_BEHAVIOR_H_
