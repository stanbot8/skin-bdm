#ifndef DERMIS_PDE_H_
#define DERMIS_PDE_H_

#include "core/pde.h"
#include "infra/sim_param.h"

namespace bdm {
namespace skibidy {

// Dermis PDE channel.
// Represents dermal tissue integrity (ECM structural state) with
// layer-specific healthy baselines: papillary (1.0), reticular (0.8),
// hypodermis (0.5). Wound destroys tissue; recovery is coupled to
// collagen deposition; MMP degrades it.
struct DermisPDE : public PDE {
  const char* GetName() const override { return fields::kDermis; }
  int GetId() const override { return fields::kDermisId; }

  void Init(Simulation* sim) override {
    auto* sp = sim->GetParam()->Get<SimParam>();
    DefineStructuralGrid(sim, sp->dermis_diffusion, sp->dermis_decay);
    // Sub-layer-aware initial profile (like ElastinPDE)
    real_t papillary_z = sp->dermal_z_papillary;
    real_t reticular_z = sp->dermal_z_reticular;
    real_t papillary_d = sp->dermis_papillary_density;
    real_t reticular_d = sp->dermis_reticular_density;
    real_t hypodermis_d = sp->dermis_hypodermis_density;
    ModelInitializer::InitializeSubstance(GetId(),
        [papillary_z, reticular_z, papillary_d, reticular_d, hypodermis_d](
            real_t x, real_t y, real_t z) -> real_t {
          if (z >= 0) return 0;                          // no dermis in epidermis
          if (z >= papillary_z) return papillary_d;      // papillary
          if (z >= reticular_z) return reticular_d;      // reticular
          return hypodermis_d;                            // hypodermis
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

#endif  // DERMIS_PDE_H_
