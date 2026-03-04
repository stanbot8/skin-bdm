#ifndef LACTATE_PDE_H_
#define LACTATE_PDE_H_

#include "core/pde.h"
#include "infra/sim_param.h"

namespace bdm {
namespace skibidy {

// Tissue lactate field: produced by anaerobic glycolysis in hypoxic wound
// tissue. Signals angiogenesis (stabilizes HIF-1a) and boosts collagen
// synthesis by fibroblasts. Cleared by vascular perfusion.
struct LactatePDE : public PDE {
  explicit LactatePDE(const SimParam* sp)
      : diffusion_(sp->lactate.diffusion),
        decay_(sp->lactate.decay) {}

  const char* GetName() const override { return fields::kLactate; }
  int GetId() const override { return fields::kLactateId; }

  void Init(Simulation* sim) override {
    DefineGrid(sim, diffusion_, decay_);
    // Starts at 0; builds from hypoxia after wounding
  }

 private:
  real_t diffusion_;
  real_t decay_;
};

}  // namespace skibidy
}  // namespace bdm

#endif  // LACTATE_PDE_H_
