#ifndef HYALURONAN_PDE_H_
#define HYALURONAN_PDE_H_

#include "core/pde.h"
#include "infra/sim_param.h"

namespace bdm {
namespace skibidy {

// Hyaluronan (HA) PDE channel.
// Primary ground substance of the dermis. Papillary dermis is HA-rich
// (water retention, cell migration scaffold); reticular dermis less so.
// None in epidermis. Wound destroys local HA.
// Production/degradation source terms are stubs for now (wired when
// fibroblast behavior is extended to produce HA).
struct HyaluronanPDE : public PDE {
  const char* GetName() const override { return fields::kHyaluronan; }
  int GetId() const override { return fields::kHyaluronanId; }

  void Init(Simulation* sim) override {
    auto* sp = sim->GetParam()->Get<SimParam>();
    DefineStructuralGrid(sim, sp->hyaluronan_diffusion, sp->hyaluronan_decay);
    // Sub-layer-aware initial profile (papillary > reticular)
    real_t papillary_z = sp->dermal_z_papillary;
    real_t reticular_z = sp->dermal_z_reticular;
    real_t papillary_d = sp->hyaluronan_basal_density;
    real_t reticular_d = sp->hyaluronan_reticular_density;
    ModelInitializer::InitializeSubstance(GetId(),
        [papillary_z, reticular_z, papillary_d, reticular_d](
            real_t x, real_t y, real_t z) -> real_t {
          if (z >= 0) return 0;                        // no HA in epidermis
          if (z >= papillary_z) return papillary_d;    // papillary: HA-rich
          if (z >= reticular_z) return reticular_d;    // reticular: moderate
          return reticular_d * 0.3;                    // hypodermis: minimal
        });
  }

  void ApplyWound(Simulation* sim, real_t, real_t, real_t) override {
    auto* grid = Grid(sim);
    auto* sp = sim->GetParam()->Get<SimParam>();
    GridContext ctx(grid, sp);
    GridContext::ZeroWoundDermis(grid, ctx);
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // HYALURONAN_PDE_H_
