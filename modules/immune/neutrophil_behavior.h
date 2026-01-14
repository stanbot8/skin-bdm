#ifndef NEUTROPHIL_BEHAVIOR_H_
#define NEUTROPHIL_BEHAVIOR_H_

#include <cmath>

#include "immune/immune_cell.h"
#include "immune/immune_helpers.h"
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
    int lifespan = sp->neutrophil_lifespan;
    if (sp->diabetic_mode) {
      lifespan = static_cast<int>(
          lifespan * sp->diabetic_neutrophil_lifespan_factor);
    }
    if (cell->GetAge() > lifespan) {
      cell->RemoveFromSimulation();
      return;
    }

    // Stochastic apoptosis after minimum survival period
    if (cell->GetAge() > sp->neutrophil_min_survival) {
      real_t apop_rate = sp->neutrophil_apoptosis_rate;
      if (sp->diabetic_mode) {
        apop_rate *= sp->diabetic_neutrophil_apoptosis_factor;
      }
      if (sim->GetRandom()->Uniform(0, 1) < apop_rate) {
        cell->RemoveFromSimulation();
        return;
      }
    }

    // Pro-inflammatory cytokine output (tapers with age -- young
    // neutrophils at the wound site are the most active producers;
    // NF-kB signaling declines as they age toward apoptosis).
    Real3 qpos = ClampToBounds(cell->GetPosition(), sim->GetParam());
    real_t cytokine_taper = std::exp(-sp->neutrophil_cytokine_taper * cell->GetAge());
    immune::ProduceProInflammatory(qpos, sim, sp, cytokine_taper);

    // MMP-8 production (dominant MMP source in acute wounds;
    // Nwomeh 1998: 100-200x more abundant than MMP-1).
    // Tapered with age like cytokine output.
    immune::ProduceMMP(qpos, sim, sp,
                       sp->mmp_neutrophil_rate * cytokine_taper);

    // Biofilm clearance
    immune::ClearBiofilm(qpos, sim, sp, kNeutrophil);

    // Migration
    immune::Migrate(cell, sim, sp);
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // NEUTROPHIL_BEHAVIOR_H_
