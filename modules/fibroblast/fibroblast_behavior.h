#ifndef FIBROBLAST_BEHAVIOR_H_
#define FIBROBLAST_BEHAVIOR_H_

#include <cmath>
#include <iostream>

#include "fibroblast/fibroblast.h"
#include "immune/immune_helpers.h"
#include "core/field_names.h"
#include "infra/sim_param.h"
#include "infra/util.h"

namespace bdm {
namespace skibidy {

// ---------------------------------------------------------------------------
// Behavior: fibroblast lifecycle, TGF-beta-driven state transitions,
// collagen deposition, and migration.
//
// States: quiescent -> activated -> myofibroblast -> removed
// Transitions driven by local TGF-beta concentration.
// Only myofibroblasts produce TGF-beta (positive feedback) and deposit
// collagen (= mechanistic scar output).
// ---------------------------------------------------------------------------
struct FibroblastBehavior : public Behavior {
  BDM_BEHAVIOR_HEADER(FibroblastBehavior, Behavior, 1);

  FibroblastBehavior() {}
  virtual ~FibroblastBehavior() {}

  void Run(Agent* agent) override {
    auto* cell = dynamic_cast<Fibroblast*>(agent);
    if (!cell) return;

    auto* sim = Simulation::GetActive();
    auto* sp = sim->GetParam()->Get<SimParam>();
    auto* rm = sim->GetResourceManager();

    cell->IncrementAge();
    cell->IncrementStateAge();
    int age = cell->GetAge();

    // --- Hard lifespan limit ---
    int eff_lifespan = sp->fibroblast_lifespan;
    if (sp->diabetic_mode) {
      eff_lifespan = static_cast<int>(
          sp->fibroblast_lifespan * sp->diabetic_fibroblast_lifespan_factor);
    }
    if (age > eff_lifespan) {
      cell->RemoveFromSimulation();
      return;
    }

    // --- Read local TGF-beta ---
    Real3 pos = cell->GetPosition();
    Real3 qpos = ClampToBounds(pos, sim->GetParam());

    auto* tgfb_grid = rm->GetDiffusionGrid(fields::kTGFBeta);
    real_t local_tgfb = tgfb_grid->GetValue(qpos);

    // --- State transitions ---
    auto state = cell->GetFibroblastState();
    int state_age = cell->GetStateAge();

    if (state == kFibroQuiescent) {
      int eff_activation_delay = sp->fibroblast_activation_delay;
      if (sp->diabetic_mode) {
        eff_activation_delay = static_cast<int>(
            sp->fibroblast_activation_delay *
            sp->diabetic_fibroblast_activation_factor);
      }
      if (local_tgfb > sp->fibroblast_activation_threshold ||
          state_age > eff_activation_delay) {
        cell->SetFibroblastState(kFibroActivated);
        if (sp->debug_fibroblast) {
          std::cout << "[fibroblast] quiescent->activated tgfb=" << local_tgfb
                    << " age=" << state_age << std::endl;
        }
      }
    } else if (state == kFibroActivated) {
      if (local_tgfb > sp->fibroblast_myofibroblast_threshold &&
          state_age > sp->fibroblast_myofibroblast_delay) {
        cell->SetFibroblastState(kMyofibroblast);
        if (sp->debug_fibroblast) {
          std::cout << "[fibroblast] activated->myofibroblast tgfb=" << local_tgfb
                    << " age=" << state_age << std::endl;
        }
      }
    } else if (state == kMyofibroblast) {
      // Immediate apoptosis if TGF-beta drops to near-zero
      if (local_tgfb < sp->fibroblast_apoptosis_threshold &&
          age > sp->fibroblast_min_lifespan) {
        cell->RemoveFromSimulation();
        return;
      }
      // Stochastic apoptosis program (Desmouliere 1995): after onset delay,
      // myofibroblasts undergo apoptosis at a per-step probability
      if (state_age > sp->fibroblast_apoptosis_onset) {
        auto* random = sim->GetRandom();
        if (random->Uniform(0, 1) < sp->fibroblast_apoptosis_rate) {
          cell->RemoveFromSimulation();
          return;
        }
      }
    }

    // --- TGF-beta production (myofibroblasts only) ---
    // Production tapers with state age (mechanical tension resolves as wound
    // contracts; Tomasek et al. 2002, doi:10.1038/nrm809).
    if (cell->GetFibroblastState() == kMyofibroblast) {
      real_t tgfb_rate = sp->fibroblast_tgfb_rate;
      // Exponential taper: halves every ~580 steps (~2.4 days)
      real_t taper = std::exp(-sp->fibroblast_tgfb_taper * cell->GetStateAge());
      tgfb_rate *= taper;
      ScaledGrid sg(tgfb_grid, sp);
      sg.AgentDeposit(sg.Index(qpos), tgfb_rate);
      if (sp->debug_fibroblast) {
        std::cout << "[fibroblast] tgfb pos=(" << pos[0] << "," << pos[1]
                  << "," << pos[2] << ") rate=" << tgfb_rate
                  << " scaled=" << tgfb_rate * sg.agent_factor << std::endl;
      }
    }

    // --- Collagen deposition (myofibroblasts only) ---
    // Myofibroblasts constitutively synthesize type I/III collagen.
    // TGF-beta drives the differentiation (gated by threshold), not the
    // synthesis rate -- once a cell IS a myofibroblast, it deposits at a
    // constant rate (Murphy et al. 2012).
    // Collagen is non-diffusing (structural ECM) -- no ScaledGrid.
    if (cell->GetFibroblastState() == kMyofibroblast) {
      auto* col_grid = rm->GetDiffusionGrid(fields::kCollagen);
      size_t col_idx = col_grid->GetBoxIndex(qpos);
      real_t deposit = sp->collagen_deposition_rate;
      if (sp->diabetic_mode) {
        deposit *= sp->diabetic_collagen_factor;
      }
      col_grid->ChangeConcentrationBy(col_idx, deposit);
    }

    // --- MMP production (activated and myofibroblast states) ---
    // Fibroblasts produce MMP-1 (collagenase) and MMP-3 (stromelysin)
    // to remodel ECM during wound repair (Nagase et al. 1999).
    if (sp->mmp_enabled &&
        cell->GetFibroblastState() != kFibroQuiescent) {
      ScaledGrid sg(rm->GetDiffusionGrid(fields::kMMP), sp);
      real_t mmp_rate = sp->mmp_fibroblast_rate;
      if (sp->diabetic_mode) {
        mmp_rate *= sp->diabetic_mmp_factor;
      }
      sg.AgentDeposit(sg.Index(qpos), mmp_rate);
    }

    // --- Fibronectin deposition (activated and myofibroblast states) ---
    // Fibroblasts deposit cellular fibronectin as provisional matrix
    // scaffold for keratinocyte migration (Clark 1990).
    // Deposition slows as local FN approaches carrying capacity
    // and as collagen accumulates (fibroblasts transition from FN
    // to collagen-dominant ECM production; Clark 1996).
    // Fibronectin is non-diffusing ECM -- no ScaledGrid.
    // Only activated fibroblasts deposit FN; myofibroblasts shift to
    // collagen-dominant production (Clark 1996).
    if (sp->fibronectin_enabled &&
        cell->GetFibroblastState() == kFibroActivated) {
      auto* fn_grid = rm->GetDiffusionGrid(fields::kFibronectin);
      size_t fn_idx = fn_grid->GetBoxIndex(qpos);
      real_t fn_cur = fn_grid->GetConcentration(fn_idx);
      real_t fn_cap = sp->fibronectin_carrying_capacity;
      real_t deposit = sp->fibronectin_deposition_rate;
      if (fn_cap > 0 && fn_cur > 0) {
        deposit *= std::max(static_cast<real_t>(0), 1.0 - fn_cur / fn_cap);
      }
      // Collagen-dependent suppression: as collagen cross-links mature,
      // fibroblasts shift from FN secretion to collagen synthesis
      real_t col_fn_switch = sp->collagen_fn_transition;
      if (col_fn_switch > 0) {
        auto* col_grid = rm->GetDiffusionGrid(fields::kCollagen);
        real_t col_local = col_grid->GetConcentration(col_grid->GetBoxIndex(qpos));
        deposit *= std::max(static_cast<real_t>(0), 1.0 - col_local / col_fn_switch);
      }
      if (deposit > 1e-10) {
        fn_grid->ChangeConcentrationBy(fn_idx, deposit);
      }
    }

    // --- Tropoelastin production (activated and myofibroblast states) ---
    // Fibroblasts secrete tropoelastin that self-assembles into elastic
    // fibers via lysyl oxidase cross-linking (Sage 1982, Kielty 2002).
    // Elastin is non-diffusing structural ECM -- no ScaledGrid.
    if (sp->elastin_enabled &&
        cell->GetFibroblastState() != kFibroQuiescent) {
      auto* el_grid = rm->GetDiffusionGrid(fields::kElastin);
      el_grid->ChangeConcentrationBy(el_grid->GetBoxIndex(qpos),
                                     sp->elastin_production_rate);
    }

    // --- Hyaluronan synthesis (activated and myofibroblast states) ---
    // Fibroblasts produce HA via HAS2 synthase. HA provides hydrated
    // ground substance for cell migration (Toole 2004, Stern 2006).
    if (sp->hyaluronan_enabled &&
        cell->GetFibroblastState() != kFibroQuiescent) {
      ScaledGrid sg(rm->GetDiffusionGrid(fields::kHyaluronan), sp);
      sg.AgentDeposit(sg.Index(qpos), sp->hyaluronan_production_rate);
    }

    // --- Migration (activated and myofibroblast states) ---
    if (cell->GetFibroblastState() == kFibroQuiescent) return;
    if (!sp->wound_enabled) return;

    // Chemotaxis: follow TGF-beta gradient
    if (sp->chemotaxis_enabled) {
      Real3 gradient = {0, 0, 0};
      tgfb_grid->GetGradient(qpos, &gradient, false);
      real_t grad_mag = std::sqrt(
          gradient[0] * gradient[0] +
          gradient[1] * gradient[1] +
          gradient[2] * gradient[2]);

      if (grad_mag > 1e-8) {
        real_t speed = sp->fibroblast_migration_speed *
                       sp->chemotaxis_speed_scale * grad_mag;
        real_t inv_mag = 1.0 / grad_mag;
        cell->SetTractorForce({
            gradient[0] * inv_mag * speed,
            gradient[1] * inv_mag * speed,
            0});
        return;
      }
    }

    // Geometric fallback: migrate toward wound center
    real_t dx = pos[0] - sp->wound_center_x;
    real_t dy = pos[1] - sp->wound_center_y;
    real_t dist = std::sqrt(dx * dx + dy * dy);
    if (dist < 1e-6) return;

    real_t inv_dist = 1.0 / dist;
    real_t dir_x = -dx * inv_dist;
    real_t dir_y = -dy * inv_dist;

    real_t r = sp->wound_radius;
    real_t speed = sp->fibroblast_migration_speed * (dist / r);

    cell->SetTractorForce({dir_x * speed, dir_y * speed, 0});
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // FIBROBLAST_BEHAVIOR_H_
