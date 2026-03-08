#ifndef MACROPHAGE_BEHAVIOR_H_
#define MACROPHAGE_BEHAVIOR_H_

#include <cmath>

#include "immune/immune_cell.h"
#include "immune/immune_helpers.h"
#include "core/continuum_handoff.h"
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
      ContinuumHandoff(cell);
      cell->RemoveFromSimulation();
      return;
    }

    // Stochastic apoptosis after minimum survival period
    if (cell->GetAge() > sp->macrophage_min_survival) {
      real_t apop_rate = sp->macrophage_apoptosis_rate;
      if (sp->diabetic_mode) {
        apop_rate *= sp->diabetic_macrophage_apoptosis_factor;
      }
      if (sim->GetRandom()->Uniform(0, 1) < apop_rate) {
        ContinuumHandoff(cell);
        cell->RemoveFromSimulation();
        return;
      }
    }

    // M2 emigration: once the wound is re-epithelialized locally,
    // M2 macrophages exit via lymphatics (Rodero et al. 2010).
    if (cell->GetState() == kM2Resolving) {
      auto* sg = sim->GetResourceManager()->GetDiffusionGrid(
          fields::kStratumId);
      if (sg) {
        Real3 pos = ClampToBounds(cell->GetPosition(), sim->GetParam());
        real_t sv = sg->GetValue(pos);
        real_t emig_rate = sp->macrophage_emigration_rate;
        if (sp->diabetic_mode) {
          emig_rate *= sp->diabetic_macrophage_emigration_factor;
        }
        if (sv >= 1.0 &&
            sim->GetRandom()->Uniform(0, 1) < emig_rate) {
          ContinuumHandoff(cell);
          cell->RemoveFromSimulation();
          return;
        }
      }
    }

    // M1 -> M2 state transition
    if (cell->GetState() == kM1Active) {
      // Efferocytosis: engulf dying neutrophils
      if (sp->efferocytosis_enabled) {
        TryEfferocytosis(cell, sim, sp);
      }

      if (cell->GetState() == kM1Active) {
        // Biofilm gate: persistent biofilm blocks M1->M2
        Real3 qpos = ClampToBounds(cell->GetPosition(), sim->GetParam());
        bool biofilm_blocked = false;
        if (sp->biofilm_enabled) {
          auto* bg = sim->GetResourceManager()->GetDiffusionGrid(
              fields::kBiofilmId);
          if (bg->GetValue(qpos) > sp->biofilm_m1_block_threshold) {
            biofilm_blocked = true;
          }
        }

        // Effective M1 duration (shared by both modes as ceiling)
        int eff_m1_duration = sp->macrophage_m1_duration;
        if (sp->diabetic_mode) {
          eff_m1_duration = static_cast<int>(
              sp->macrophage_m1_duration * sp->diabetic_m1_duration_factor);
        }
        // AGE-RAGE prolongation
        if (sp->diabetic_mode && sp->glucose_enabled &&
            sp->age_m1_prolongation > 0) {
          auto* age_grid = sim->GetResourceManager()->GetDiffusionGrid(
              fields::kAGEId);
          if (age_grid) {
            real_t age_val = age_grid->GetValue(qpos);
            eff_m1_duration = static_cast<int>(
                eff_m1_duration * (1.0 + sp->age_m1_prolongation * age_val));
          }
        }

        bool should_transition = false;
        if (sp->mech_m1_m2_transition) {
          // Mechanistic: efferocytosis-count driven with timer ceiling.
          // Engulfing apoptotic neutrophils is the primary M2 signal
          // (PS receptor -> NR4A1 -> anti-inflammatory program). Once
          // quota is met, transition fires stochastically. Timer ceiling
          // models alternative M2 signals (IL-4, IL-10, glucocorticoids).
          if (cell->GetEngulfCount() >= sp->mech_efferocytosis_quota) {
            should_transition =
                (sim->GetRandom()->Uniform(0, 1) < sp->mech_m2_transition_rate);
          }
          // Timer ceiling fallback (alternative M2 polarization signals)
          if (cell->GetStateAge() > eff_m1_duration) {
            should_transition = true;
          }
        } else {
          // Parametric: cytokine threshold + timer ceiling
          real_t infl = GetNetInflammation(sim, qpos);
          bool cytokine_trigger =
              (cell->GetStateAge() > sp->m1_transition_min_age &&
               infl < sp->m1_transition_threshold);
          bool timer_trigger = (cell->GetStateAge() > eff_m1_duration);
          should_transition = (cytokine_trigger || timer_trigger);
        }

        if (should_transition && !biofilm_blocked) {
          cell->SetState(kM2Resolving);
        }
      }
    }

    // Cytokine output (M1 production tapers with state age -- initial
    // NF-kB burst is strongest; chronic M1 cells are less productive).
    Real3 qpos = ClampToBounds(cell->GetPosition(), sim->GetParam());
    auto* rm = sim->GetResourceManager();
    if (cell->GetState() == kM1Active) {
      real_t m1_taper = std::exp(-sp->m1_cytokine_taper * cell->GetStateAge());
      immune::ProduceProInflammatory(qpos, sim, sp, m1_taper);
      immune::ProduceMMP(qpos, sim, sp);
      // M1 macrophages produce NO via iNOS (Witte & Barbul 2002)
      if (sp->nitric_oxide_enabled) {
        auto* no_grid = rm->GetDiffusionGrid(fields::kNitricOxideId);
        if (no_grid) {
          real_t no_rate = sp->no_m1_production * m1_taper;
          if (sp->diabetic_mode) no_rate *= sp->diabetic_no_factor;
          ScaledGrid no_sg(no_grid, sp);
          no_sg.AgentDeposit(no_sg.grid->GetBoxIndex(qpos), no_rate);
        }
      }
    } else {
      immune::ProduceAntiInflammatory(qpos, sim, sp);
      real_t m2_taper = std::exp(-sp->m2_tgfb_taper * cell->GetStateAge());
      immune::ProduceTGFBeta(qpos, sim, sp, m2_taper);
      immune::ProduceVEGF(qpos, sim, sp);
      // M2 macrophages produce TIMP-1 as part of pro-resolution program.
      // Mantovani et al. 2004 (doi:10.1016/j.it.2004.09.015)
      immune::ProduceTIMP(qpos, sim, sp, sp->timp_macrophage_rate * m2_taper);
    }

    // Receptor-mediated TGF-beta endocytosis (both M1 and M2 have TbRII/TbRI)
    if (sp->fibroblast_enabled && sp->tgfb_receptor_consumption > 0) {
      real_t local_tgfb = rm->GetDiffusionGrid(fields::kTGFBetaId)->GetValue(qpos);
      immune::ConsumeTGFBeta(qpos, sim, sp, local_tgfb);
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
      cell->IncrementEngulfCount();
      if (!sp->mech_m1_m2_transition) {
        // Parametric mode: immediate M2 transition on first engulfment
        cell->SetState(kM2Resolving);
      }
      // Mechanistic mode: count accumulates, transition handled above
    }
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // MACROPHAGE_BEHAVIOR_H_
