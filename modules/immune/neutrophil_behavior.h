#ifndef NEUTROPHIL_BEHAVIOR_H_
#define NEUTROPHIL_BEHAVIOR_H_

#include <cmath>

#include "immune/immune_cell.h"
#include "immune/immune_helpers.h"
#include "core/continuum_handoff.h"
#include "infra/sim_param.h"
#include "infra/util.h"

namespace bdm {
namespace skibidy {

// ---------------------------------------------------------------------------
// Neutrophil behavior: pro-inflammatory cytokine production, fixed lifespan,
// migration toward wound center (geometric or chemotactic).
// ---------------------------------------------------------------------------
struct NeutrophilBehavior : public Behavior {
  BDM_BEHAVIOR_HEADER(NeutrophilBehavior, Behavior, 1);

  NeutrophilBehavior() {}
  virtual ~NeutrophilBehavior() {}

  void Run(Agent* agent) override {
    auto* cell = dynamic_cast<ImmuneCell*>(agent);
    if (!cell) return;

    auto* sim = Simulation::GetActive();
    auto* sp = sim->GetParam()->Get<SimParam>();

    cell->IncrementAge();
    cell->IncrementStateAge();

    // Hard ceiling lifespan
    int lifespan = sp->immune.neutrophil_lifespan;
    if (sp->diabetic.mode) {
      lifespan = static_cast<int>(
          lifespan * sp->diabetic.neutrophil_lifespan_factor);
    }
    if (cell->GetAge() > lifespan) {
      ContinuumHandoff(cell);
      cell->RemoveFromSimulation();
      return;
    }

    // Stochastic apoptosis after minimum survival period
    if (cell->GetAge() > sp->immune.neutrophil_min_survival) {
      real_t apop_rate = sp->immune.neutrophil_apoptosis_rate;
      if (sp->diabetic.mode) {
        apop_rate *= sp->diabetic.neutrophil_apoptosis_factor;
      }
      if (sim->GetRandom()->Uniform(0, 1) < apop_rate) {
        ContinuumHandoff(cell);
        cell->RemoveFromSimulation();
        return;
      }
    }

    // Pro-inflammatory cytokine output (tapers with age -- young
    // neutrophils at the wound site are the most active producers;
    // NF-kB signaling declines as they age toward apoptosis).
    Real3 qpos = ClampToBounds(cell->GetPosition(), sim->GetParam());
    real_t cytokine_taper = std::exp(-sp->immune.neutrophil_cytokine_taper * cell->GetAge());
    immune::ProduceProInflammatory(qpos, sim, sp, cytokine_taper);

    // MMP-8 production (dominant MMP source in acute wounds;
    // Nwomeh 1998: 100-200x more abundant than MMP-1).
    // Tapered with age like cytokine output.
    immune::ProduceMMP(qpos, sim, sp,
                       sp->mmp.neutrophil_rate * cytokine_taper);

    // MMP-9 degranulation: pre-formed gelatinase B released from secondary
    // granules directly into the active MMP pool on arrival (no pro-MMP
    // activation delay). Kolaczkowska & Kubes 2013.
    if (sp->mmp.enabled && cell->GetAge() < sp->mmp.degranulation_window) {
      auto* mmp_grid = sim->GetResourceManager()->GetDiffusionGrid(
          fields::kMMPId);
      if (mmp_grid) {
        real_t burst = sp->mmp.neutrophil_degranulation;
        if (sp->diabetic.mode) burst *= sp->diabetic.mmp_factor;
        ScaledGrid sg(mmp_grid, sp);
        sg.AgentDeposit(sg.Index(qpos), burst);
      }
    }

    // NO production via iNOS (Witte & Barbul 2002)
    if (sp->nitric_oxide.enabled) {
      auto* rm = sim->GetResourceManager();
      auto* no_grid = rm->GetDiffusionGrid(fields::kNitricOxideId);
      if (no_grid) {
        real_t no_rate = sp->no_neutrophil_production * cytokine_taper;
        if (sp->diabetic.mode) no_rate *= sp->diabetic.no_factor;
        ScaledGrid no_sg(no_grid, sp);
        no_sg.AgentDeposit(no_sg.grid->GetBoxIndex(qpos), no_rate);
      }
    }

    // ROS production via NADPH oxidase respiratory burst (Babior 2000)
    if (sp->ros.enabled) {
      auto* rm2 = sim->GetResourceManager();
      auto* ros_grid = rm2->GetDiffusionGrid(fields::kROSId);
      if (ros_grid) {
        real_t ros_rate = sp->ros.neutrophil_burst * cytokine_taper;
        ScaledGrid ros_sg(ros_grid, sp);
        ros_sg.AgentDeposit(ros_sg.grid->GetBoxIndex(qpos), ros_rate);
      }
    }

    // Biofilm clearance
    immune::ClearBiofilm(qpos, sim, sp, kNeutrophil);

    // Migration
    immune::Migrate(cell, sim, sp);
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // NEUTROPHIL_BEHAVIOR_H_
