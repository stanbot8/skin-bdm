#ifndef FIBRONECTIN_PDE_H_
#define FIBRONECTIN_PDE_H_

#include "core/pde.h"
#include "infra/sim_param.h"

namespace bdm {
namespace skibidy {

struct FibronectinPDE : public SimpleStructuralPDE {
  explicit FibronectinPDE(const SimParam* sp)
      : SimpleStructuralPDE(fields::kFibronectin, fields::kFibronectinId,
                            0, sp->fibronectin_decay) {}

  // Plasma fibronectin co-deposits with fibrin in the provisional clot
  // (Clark 1996, doi:10.1007/978-1-4615-1795-5).
  void ApplyWound(Simulation* sim, real_t, real_t, real_t) override {
    auto* grid = Grid(sim);
    auto* sp = sim->GetParam()->Get<SimParam>();
    real_t seed = sp->fibronectin_wound_seed;
    if (seed <= 0) return;
    GridContext ctx(grid, sp);
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

#endif  // FIBRONECTIN_PDE_H_
