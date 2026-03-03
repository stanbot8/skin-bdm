#ifndef EDEMA_PDE_H_
#define EDEMA_PDE_H_

#include "core/pde.h"
#include "infra/sim_param.h"

namespace bdm {
namespace skibidy {

// Interstitial edema field: accumulates from inflammatory vascular leak
// (Starling forces) and drains via lymphatic vessels. Impairs O2 delivery
// (diffusion barrier) and cell migration (hydrostatic resistance).
// Replaces the old parametric lymphatic sigmoid with a mechanistic
// source/sink PDE that couples inflammation, perfusion, and lymphatics.
struct EdemaPDE : public PDE {
  explicit EdemaPDE(const SimParam*) {}

  const char* GetName() const override { return fields::kEdema; }
  int GetId() const override { return fields::kEdemaId; }

  void Init(Simulation* sim) override {
    // D=0: edema is local interstitial fluid, does not diffuse laterally.
    // Decay=0: drained mechanistically by lymphatic field in fused_source.
    DefineGrid(sim, 0.0, 0.0);
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // EDEMA_PDE_H_
