#ifndef SCAB_PDE_H_
#define SCAB_PDE_H_

#include "core/pde.h"
#include "infra/sim_param.h"

namespace bdm {
namespace skibidy {

// Scab field: dried fibrin/platelet/exudate crust over wound surface.
// D=0 (structural, immobile), no FTCS decay (manual decay via source hook).
// Seeded on wound creation in epidermal wound voxels only.
// Eming et al. 2007 (doi:10.1038/sj.jid.5700701)
struct ScabPDE : public SimpleStructuralPDE {
  explicit ScabPDE(const SimParam* sp)
      : SimpleStructuralPDE(fields::kScab, fields::kScabId, 0, 0) {}

  // Wound creates scab: seed in epidermal wound cylinder only (surface crust).
  void ApplyWound(Simulation* sim, real_t, real_t, real_t) override {
    auto* grid = Grid(sim);
    auto* sp = sim->GetParam()->Get<SimParam>();
    GridContext ctx(grid, sp);
    real_t seed = sp->scab.wound_seed;
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

#endif  // SCAB_PDE_H_
