#ifndef SHEDDING_H_
#define SHEDDING_H_

#include "tissue/keratinocyte.h"
#include "infra/sim_param.h"

namespace bdm {
namespace skibidy {

// ---------------------------------------------------------------------------
// Behavior: desquamation of cornified cells after residence time, and
// apoptosis of exhausted TA cells that remain in the basal layer.
// ---------------------------------------------------------------------------
struct Shedding : public Behavior {
  BDM_BEHAVIOR_HEADER(Shedding, Behavior, 1);

  Shedding() { AlwaysCopyToNew(); }
  virtual ~Shedding() {}

  void Run(Agent* agent) override {
    auto* cell = dynamic_cast<Keratinocyte*>(agent);
    if (!cell) return;

    auto* sim = Simulation::GetActive();
    auto* sp = sim->GetParam()->Get<SimParam>();

    // Tissue boundary: remove cells pushed past lateral (x/y) edges
    auto pos = cell->GetPosition();
    if (pos[0] < sp->tissue_min || pos[0] > sp->tissue_max ||
        pos[1] < sp->tissue_min || pos[1] > sp->tissue_max) {
      cell->RemoveFromSimulation();
      return;
    }

    if (cell->GetStratum() == kCornified) {
      steps_cornified_++;
      if (steps_cornified_ > sp->shedding_delay) {
        cell->RemoveFromSimulation();
        return;
      }
    } else {
      steps_cornified_ = 0;
    }

    // Apoptosis: exhausted TA cells lingering in basal layer are cleared
    if (cell->GetStratum() == kBasal && !cell->IsStem() &&
        cell->GetDivisionsLeft() <= 0 &&
        cell->GetStratumSteps() > sp->apoptosis_delay) {
      cell->RemoveFromSimulation();
    }
  }

 private:
  int steps_cornified_ = 0;
};

}  // namespace skibidy
}  // namespace bdm

#endif  // SHEDDING_H_
