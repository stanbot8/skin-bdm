#ifndef TGFBETA_PDE_H_
#define TGFBETA_PDE_H_

#include "core/pde.h"
#include "infra/sim_param.h"

namespace bdm {
namespace skibidy {

struct TGFBetaPDE : public SimplePDE {
  explicit TGFBetaPDE(const SimParam* sp)
      : SimplePDE(fields::kTGFBeta, fields::kTGFBetaId,
                  sp->tgfb_diffusion, sp->tgfb_decay) {}

  // Platelet alpha-granule degranulation: TGF-beta1 released from clot
  // within minutes of hemostasis (Shah et al. 1995, doi:10.1016/S0140-6736(95)90124-8).
  void ApplyWound(Simulation* sim, real_t, real_t, real_t) override {
    auto* grid = Grid(sim);
    auto* sp = sim->GetParam()->Get<SimParam>();
    real_t seed = sp->tgfb_wound_seed;
    if (seed <= 0) return;
    GridContext ctx(grid, sp);
    for (size_t idx = 0; idx < ctx.n; idx++) {
      if (ctx.InWound(ctx.X(idx), ctx.Y(idx))) {
        real_t z = ctx.Z(idx);
        if (z >= sp->dermal_z_reticular && z <= sp->volume_z_cornified) {
          grid->ChangeConcentrationBy(idx, seed);
        }
      }
    }
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // TGFBETA_PDE_H_
