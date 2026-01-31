#ifndef PH_PDE_H_
#define PH_PDE_H_

#include "core/pde.h"
#include "infra/sim_param.h"

namespace bdm {
namespace skibidy {

// Wound alkalinity field: 0 = normal acidic skin (pH ~5.5),
// 1.0 = fresh wound (alkaline, pH ~7.4).
// D=0 (no lateral diffusion), decay=0 (recovery via fused source).
struct PHPDE : public SimpleStructuralPDE {
  explicit PHPDE(const SimParam*)
      : SimpleStructuralPDE(fields::kPH, fields::kPHId, 0, 0) {}

  // Wound disrupts acid mantle: set alkalinity to 1.0 in wound cylinder.
  void ApplyWound(Simulation* sim, real_t cx, real_t cy, real_t r) override {
    auto* grid = Grid(sim);
    auto* sp = sim->GetParam()->Get<SimParam>();
    GridContext ctx(grid, sp);
    for (size_t idx = 0; idx < ctx.n; idx++) {
      if (ctx.InWound(ctx.X(idx), ctx.Y(idx))) {
        real_t current = grid->GetConcentration(idx);
        if (current < 1.0) {
          grid->ChangeConcentrationBy(idx, 1.0 - current);
        }
      }
    }
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // PH_PDE_H_
