#ifndef BASAL_DIVISION_H_
#define BASAL_DIVISION_H_

#include <cmath>

#include "tissue/keratinocyte.h"
#include "core/field_names.h"
#include "infra/sim_param.h"
#include "infra/util.h"

namespace bdm {
namespace skibidy {

// ---------------------------------------------------------------------------
// Behavior: explicit G1->S->G2->M cell cycle for basal keratinocytes
// with contact inhibition and asymmetric stem cell division.
// Non-basal cells and exhausted TA cells enter quiescence (G0).
// ---------------------------------------------------------------------------
struct BasalDivision : public Behavior {
  BDM_BEHAVIOR_HEADER(BasalDivision, Behavior, 1);

  BasalDivision() { AlwaysCopyToNew(); }
  virtual ~BasalDivision() {}

  void Run(Agent* agent) override {
    auto* cell = dynamic_cast<Keratinocyte*>(agent);
    if (!cell) return;

    auto* sim = Simulation::GetActive();
    auto* sp = sim->GetParam()->Get<SimParam>();
    real_t dt = sim->GetParam()->simulation_time_step;

    // Track cell age
    cell->SetCellAge(cell->GetCellAge() + dt);

    // Non-basal cells or exhausted TA cells are quiescent
    if (cell->GetStratum() != kBasal ||
        (!cell->IsStem() && cell->GetDivisionsLeft() <= 0)) {
      cell->SetCyclePhase(kG0);
      return;
    }

    cell->SetPhaseElapsed(cell->GetPhaseElapsed() + dt);
    auto* random = sim->GetRandom();
    real_t ran = random->Uniform(0, 1);

    switch (cell->GetCyclePhase()) {
      case kG1: {
        // Probabilistic transition: cumulative probability over time
        real_t p = cell->GetPhaseElapsed() / sp->g1_duration;

        // KGF-dependent proliferation boost (Michaelis-Menten kinetics).
        // KGF is high at basement membrane (dermis source) and decays
        // upward, so only basal cells near z=0 get the full signal.
        auto* rm = sim->GetResourceManager();
        auto* kgf_grid = rm->GetDiffusionGrid(fields::kKGF);
        Real3 qpos = ClampToBounds(cell->GetPosition(), sim->GetParam());
        real_t kgf = kgf_grid->GetValue(qpos);
        real_t kgf_boost = 1.0 + sp->kgf_max_boost *
                           kgf / (sp->kgf_half_maximal + kgf);
        p *= kgf_boost;

        // O2-dependent proliferation: hypoxia suppresses G1->S transition.
        // Below oxygen_prolif_threshold, probability scales linearly to zero.
        auto* o2_grid = rm->GetDiffusionGrid(fields::kOxygen);
        real_t o2 = o2_grid->GetValue(qpos);
        real_t o2_factor = std::min(1.0, o2 / sp->oxygen_prolif_threshold);
        if (o2_factor < 0) o2_factor = 0;
        p *= o2_factor;

        // Water-dependent proliferation: dehydrated cells cannot divide.
        // Below water_prolif_threshold, probability scales linearly to zero.
        auto* water_grid = rm->GetDiffusionGrid(fields::kWater);
        real_t water = water_grid->GetValue(qpos);
        real_t water_factor = std::min(1.0, water / sp->water_prolif_threshold);
        if (water_factor < 0) water_factor = 0;
        p *= water_factor;

        // Inflammation gating: Hill-function dose-response.
        // Normal: quadratic Hill K^2/(K^2+I^2) for sharp threshold.
        // Diabetic: linear Hill K/(K+I) for gradual, persistent
        // suppression reflecting chronic sensitivity (Rasik & Shukla 2000).
        real_t infl = GetImmunePressure(sim, qpos);
        real_t eff_infl = infl;
        if (sp->diabetic_mode) {
          eff_infl *= sp->diabetic_inflammation_sensitivity;
        }
        real_t infl_factor = 1.0;
        if (sp->diabetic_mode) {
          real_t K = sp->diabetic_prolif_infl_K;
          if (K + eff_infl > 1e-12) {
            infl_factor = K / (K + eff_infl);
          }
        } else {
          real_t K = sp->inflammation_prolif_threshold;
          real_t K2 = K * K;
          real_t denom = K2 + eff_infl * eff_infl;
          if (denom > 1e-12) {
            infl_factor = K2 / denom;
          }
        }
        p *= infl_factor;

        // Diabetic mode: hyperglycemia impairs keratinocyte proliferation
        if (sp->diabetic_mode) {
          p *= sp->diabetic_prolif_factor;
        }

        if (p > ran) {
          cell->SetCyclePhase(kS);
          cell->SetPhaseElapsed(0);
        }
        break;
      }
      case kS: {
        // S phase: DNA synthesis + cell growth (capped at division size)
        if (cell->GetDiameter() < sp->division_diameter) {
          cell->ChangeVolume(sp->growth_rate);
        }
        real_t p = cell->GetPhaseElapsed() / sp->s_duration;
        if (p > ran) {
          cell->SetCyclePhase(kG2);
          cell->SetPhaseElapsed(0);
        }
        break;
      }
      case kG2: {
        real_t p = cell->GetPhaseElapsed() / sp->g2_duration;
        if (p > ran) {
          cell->SetCyclePhase(kM);
          cell->SetPhaseElapsed(0);
        }
        break;
      }
      case kM: {
        real_t p = cell->GetPhaseElapsed() / sp->m_duration;
        if (p > ran) {
          // Mitosis complete — attempt division
          bool divided = false;
          if (cell->GetDiameter() >= sp->division_diameter) {
            // Contact inhibition check
            int neighbors = 0;
            auto* ctxt = sim->GetExecutionContext();
            auto count = L2F([&](Agent*, real_t) { neighbors++; });
            real_t d = cell->GetDiameter();
            ctxt->ForEachNeighbor(count, *cell, d * d);

            if (neighbors < sp->max_neighbors) {
              // Inward-biased division toward wound center.
              // Blend random theta with center-directed theta.
              real_t theta_random = random->Uniform(0, 2.0 * M_PI);
              real_t theta;
              if (sp->wound_enabled && sp->wound_inward_bias > 0) {
                Real3 cpos = cell->GetPosition();
                real_t to_cx = sp->wound_center_x - cpos[0];
                real_t to_cy = sp->wound_center_y - cpos[1];
                real_t dist_c = std::sqrt(to_cx * to_cx + to_cy * to_cy);
                if (dist_c > 1e-6) {
                  real_t theta_center = std::atan2(to_cy, to_cx);
                  real_t diff = theta_random - theta_center;
                  while (diff > M_PI) diff -= 2.0 * M_PI;
                  while (diff < -M_PI) diff += 2.0 * M_PI;
                  theta = theta_center +
                          diff * (1.0 - sp->wound_inward_bias);
                } else {
                  theta = theta_random;
                }
              } else {
                theta = theta_random;
              }
              real_t dx = std::cos(theta);
              real_t dy = std::sin(theta);

              // Crowding-driven vertical extrusion: when basal layer
              // is crowded, push daughters upward for stratification.
              // Threshold 0.3 so extrusion triggers in narrow wounds
              // where cells see only 4-6 neighbors out of max 14.
              real_t crowding = static_cast<real_t>(neighbors) /
                                static_cast<real_t>(sp->max_neighbors);
              real_t dz;
              if (crowding > 0.3) {
                real_t extrusion = (crowding - 0.3) / 0.7;
                dz = 0.5 + 1.5 * extrusion;
              } else {
                dz = random->Uniform(0, 0.2);
              }

              cell->Divide({dx, dy, dz});
              divided = true;
              // Daughter inherits via Initialize():
              //   stem mother -> TA daughter (divisions_left = ta_max)
              //   TA mother -> TA daughter (divisions_left - 1)
              // Decrement mother's count (stem cells stay at -1)
              if (!cell->IsStem()) {
                cell->SetDivisionsLeft(cell->GetDivisionsLeft() - 1);
              }
            }
          }
          if (divided) {
            cell->SetCyclePhase(kG1);
          } else {
            // Contact inhibition or insufficient size — enter quiescence
            cell->SetCyclePhase(kG0);
          }
          cell->SetPhaseElapsed(0);
        }
        break;
      }
      default:  // kG0 — shouldn't reach here, but re-enter cycle
        cell->SetCyclePhase(kG1);
        cell->SetPhaseElapsed(0);
        break;
    }
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // BASAL_DIVISION_H_
