#ifndef ELASTIN_PDE_H_
#define ELASTIN_PDE_H_

#include "core/pde.h"
#include "infra/sim_param.h"

namespace bdm {
namespace skibidy {

// Elastin PDE channel.
// Elastic fiber network in the dermis. Reticular dermis has denser fibers
// than papillary dermis; none in epidermis. Wound destroys local elastin.
// Production/degradation source terms are stubs for now (wired when
// fibroblast behavior is extended to produce tropoelastin).
struct ElastinPDE : public PDE {
  const char* GetName() const override { return fields::kElastin; }
  int GetId() const override { return fields::kElastinId; }

  void Init(Simulation* sim) override {
    auto* sp = sim->GetParam()->Get<SimParam>();
    DefineStructuralGrid(sim, sp->elastin_diffusion, sp->elastin_decay);
    // Sub-layer-aware initial profile
    real_t papillary_z = sp->dermal_z_papillary;
    real_t reticular_z = sp->dermal_z_reticular;
    real_t papillary_d = sp->elastin_papillary_density;
    real_t reticular_d = sp->elastin_basal_density;
    ModelInitializer::InitializeSubstance(GetId(),
        [papillary_z, reticular_z, papillary_d, reticular_d](
            real_t x, real_t y, real_t z) -> real_t {
          if (z >= 0) return 0;                        // no elastin in epidermis
          if (z >= papillary_z) return papillary_d;    // papillary: thin fibers
          if (z >= reticular_z) return reticular_d;    // reticular: dense fibers
          return reticular_d * 0.5;                    // hypodermis: sparse
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

#endif  // ELASTIN_PDE_H_
