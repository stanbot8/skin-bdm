#ifndef STIFFNESS_PDE_H_
#define STIFFNESS_PDE_H_

#include "core/pde.h"
#include "infra/sim_param.h"

namespace bdm {
namespace skibidy {

// Tissue stiffness field: derived from ECM composition (collagen + elastin).
// Non-diffusing structural property. Computed each step in fused_source from
// local collagen and elastin concentrations. Drives YAP/TAZ mechanosensing
// in fibroblasts (Tomasek et al. 2002) and modulates scar formation.
struct StiffnessPDE : public PDE {
  explicit StiffnessPDE(const SimParam*) {}

  const char* GetName() const override { return fields::kStiffness; }
  int GetId() const override { return fields::kStiffnessId; }

  void Init(Simulation* sim) override {
    // D=0, decay=0: stiffness is a derived structural property, not a PDE.
    // Updated each step from collagen + elastin in fused_source.
    DefineGrid(sim, 0.0, 0.0);
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // STIFFNESS_PDE_H_
