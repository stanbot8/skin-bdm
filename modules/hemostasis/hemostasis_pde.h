#ifndef HEMOSTASIS_PDE_H_
#define HEMOSTASIS_PDE_H_

#include "core/pde.h"
#include "infra/sim_param.h"

namespace bdm {
namespace skibidy {

// Fibrin clot field: provisional wound matrix from platelet aggregation.
// D=0 (structural ECM, immobile), decay via fused_source manual decay.
// Seeded on wound creation; degraded by MMP (fibrinolysis).
struct FibrinPDE : public SimpleStructuralPDE {
  explicit FibrinPDE(const SimParam* sp)
      : SimpleStructuralPDE(fields::kFibrin, fields::kFibrinId, 0, 0) {}

  // Wound creates fibrin clot: seed in wound cylinder.
  void ApplyWound(Simulation* sim, real_t, real_t, real_t) override {
    auto* grid = Grid(sim);
    auto* sp = sim->GetParam()->Get<SimParam>();
    GridContext ctx(grid, sp);
    real_t seed = sp->hemostasis_wound_seed;
    for (size_t idx = 0; idx < ctx.n; idx++) {
      if (ctx.InWound(ctx.X(idx), ctx.Y(idx))) {
        real_t z = ctx.Z(idx);
        if (z >= 0 && z <= sp->volume_z_cornified) {
          grid->ChangeConcentrationBy(idx, seed);
        }
      }
    }
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // HEMOSTASIS_PDE_H_
