#ifndef TCELL_DENSITY_PDE_H_
#define TCELL_DENSITY_PDE_H_

#include "core/pde.h"
#include "infra/sim_param.h"

namespace bdm {
namespace skibidy {

// T cell density field: adaptive immunity proxy representing Th1/Th17
// infiltration into RA synovial tissue. Recruited by inflammation,
// proliferates via antigen-driven clonal expansion (flare-gated),
// and modulates autoimmune cytokine production in source_hook.
// Suppressed by methotrexate (purine synthesis inhibition) and
// JAK inhibitors (IL-2/IL-15 signaling blockade).
// McInnes & Schett 2011 (doi:10.1056/NEJMra1004965)
struct TCellDensityPDE : public PDE {
  explicit TCellDensityPDE(const SimParam* sp)
      : diffusion_(sp->ra.tcell_diffusion),
        decay_(sp->ra.tcell_decay) {}

  const char* GetName() const override { return fields::kTCellDensity; }
  int GetId() const override { return fields::kTCellDensityId; }

  void Init(Simulation* sim) override {
    DefineGrid(sim, diffusion_, decay_);
    // Starts at 0; T cells recruited post-wound by inflammation
  }

 private:
  real_t diffusion_;
  real_t decay_;
};

}  // namespace skibidy
}  // namespace bdm

#endif  // TCELL_DENSITY_PDE_H_
