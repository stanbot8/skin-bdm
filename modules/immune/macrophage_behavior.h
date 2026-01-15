#ifndef MACROPHAGE_BEHAVIOR_H_
#define MACROPHAGE_BEHAVIOR_H_

#include <cmath>

#include "immune/immune_cell.h"
#include "immune/immune_helpers.h"
#include "core/field_names.h"
#include "infra/sim_param.h"
#include "infra/util.h"

namespace bdm {
namespace skibidy {

// ---------------------------------------------------------------------------
// Macrophage behavior: M1 (pro-inflammatory) -> M2 (pro-resolution)
// transition via timer or efferocytosis, cytokine production, TGF-beta
// secretion, and migration toward wound center.
// ---------------------------------------------------------------------------
struct MacrophageBehavior : public Behavior {
  BDM_BEHAVIOR_HEADER(MacrophageBehavior, Behavior, 1);

  MacrophageBehavior() {}
  virtual ~MacrophageBehavior() {}

  void Run(Agent* agent) override {
    auto* cell = dynamic_cast<ImmuneCell*>(agent);
    if (!cell) return;

    auto* sim = Simulation::GetActive();
    auto* sp = sim->GetParam()->Get<SimParam>();

    cell->IncrementAge();
    cell->IncrementStateAge();

    // Hard ceiling lifespan
    if (cell->GetAge() > sp->macrophage_lifespan) {
      cell->RemoveFromSimulation();
      return;
    }

    // Stochastic apoptosis after minimum survival period
    if (cell->GetAge() > sp->macrophage_min_survival) {
      if (sim->GetRandom()->Uniform(0, 1) < sp->macrophage_apoptosis_rate) {
        cell->RemoveFromSimulation();
        return;
      }
    }

    // M2 emigration: once the wound is re-epithelialized locally,
    // M2 macrophages exit via lymphatics (Rodero et al. 2010).
    if (cell->GetState() == kM2Resolving) {
      auto* sg = sim->GetResourceManager()->GetDiffusionGrid(
          fields::kStratum);
      if (sg) {
        Real3 pos = ClampToBounds(cell->GetPosition(), sim->GetParam());
        real_t sv = sg->GetValue(pos);
        real_t emig_rate = sp->macrophage_emigration_rate;
        if (sp->diabetic_mode) {
          emig_rate *= sp->diabetic_macrophage_emigration_factor;
        }
        if (sv >= 1.0 &&
            sim->GetRandom()->Uniform(0, 1) < emig_rate) {
          cell->RemoveFromSimulation();
          return;
        }
      }
    }

    // M1 -> M2 state transition
    if (cell->GetState() == kM1Active) {
      // Efferocytosis: engulf dying neutrophils for early M2 transition
      if (sp->efferocytosis_enabled) {
        TryEfferocytosis(cell, sim, sp);
      }

      if (cell->GetState() == kM1Active) {
        int eff_m1_duration = sp->macrophage_m1_duration;
        if (sp->diabetic_mode) {
          eff_m1_duration = static_cast<int>(
              sp->macrophage_m1_duration * sp->diabetic_m1_duration_factor);
        }

        // Cytokine-driven: transition when local inflammation drops
        Real3 qpos = ClampToBounds(cell->GetPosition(), sim->GetParam());
        real_t infl = GetNetInflammation(sim, qpos);
        bool cytokine_trigger =
            (cell->GetStateAge() > sp->m1_transition_min_age &&
             infl < sp->m1_transition_threshold);

        // Timer ceiling fallback
        bool timer_trigger = (cell->GetStateAge() > eff_m1_duration);

        // Biofilm gate: persistent biofilm blocks M1->M2
        bool biofilm_blocked = false;
        if (sp->biofilm_enabled) {
          auto* bg = sim->GetResourceManager()->GetDiffusionGrid(
              fields::kBiofilm);
          if (bg->GetValue(qpos) > sp->biofilm_m1_block_threshold) {
            biofilm_blocked = true;
          }
        }

        if ((cytokine_trigger || timer_trigger) && !biofilm_blocked) {
          cell->SetState(kM2Resolving);
        }
      }
    }

    // Cytokine output (M1 production tapers with state age -- initial
    // NF-kB burst is strongest; chronic M1 cells are less productive).
    Real3 qpos = ClampToBounds(cell->GetPosition(), sim->GetParam());
    if (cell->GetState() == kM1Active) {
      real_t m1_taper = std::exp(-sp->m1_cytokine_taper * cell->GetStateAge());
      immune::ProduceProInflammatory(qpos, sim, sp, m1_taper);
      immune::ProduceMMP(qpos, sim, sp);
    } else {
      immune::ProduceAntiInflammatory(qpos, sim, sp);
      real_t m2_taper = std::exp(-sp->m2_tgfb_taper * cell->GetStateAge());
      immune::ProduceTGFBeta(qpos, sim, sp, m2_taper);
      immune::ProduceVEGF(qpos, sim, sp);
    }

    // Biofilm clearance (both M1 and M2 contribute)
    immune::ClearBiofilm(qpos, sim, sp, kMacrophage);

    // Migration
    immune::Migrate(cell, sim, sp);
  }

 private:
  void TryEfferocytosis(ImmuneCell* cell, Simulation* sim,
                         const SimParam* sp) {
    auto* ctxt = sim->GetExecutionContext();
    int age_threshold = static_cast<int>(
        sp->neutrophil_lifespan * sp->efferocytosis_age_fraction);
    bool engulfed = false;
    real_t r2 = sp->efferocytosis_radius * sp->efferocytosis_radius;

    auto detect = L2F([&](Agent* neighbor, real_t) {
      if (engulfed) return;
      auto* neutrophil = dynamic_cast<ImmuneCell*>(neighbor);
      if (!neutrophil) return;
      if (neutrophil->GetImmuneCellType() != kNeutrophil) return;
      if (neutrophil->GetAge() < age_threshold) return;

      // Diabetic mode reduces efferocytosis probability
      if (sp->diabetic_mode) {
        auto* random = sim->GetRandom();
        if (random->Uniform(0, 1) > sp->diabetic_efferocytosis_factor) {
          return;
        }
      }

      // Don't call RemoveFromSimulation() from a neighbor callback --
      // it races with the neutrophil's own behavior in the same parallel
      // step, causing double-removal crashes in TypeIndex::Remove.
      // Instead, force the neutrophil past its lifespan so it removes
      // itself through its own behavior.
      neutrophil->SetAge(999999);
      engulfed = true;
    });
    ctxt->ForEachNeighbor(detect, *cell, r2);

    if (engulfed) {
      cell->SetState(kM2Resolving);
    }
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // MACROPHAGE_BEHAVIOR_H_
