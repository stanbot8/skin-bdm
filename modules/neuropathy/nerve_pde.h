#ifndef NERVE_PDE_H_
#define NERVE_PDE_H_

#include "core/pde.h"

namespace bdm {
namespace skibidy {

// Nerve density field: represents cutaneous sensory innervation density.
// Initialized at basal level in healthy tissue, zeroed in wound area
// (Wallerian degeneration), and slowly regenerates from wound margins
// guided by Schwann cells and neurotrophic factors.
// Boulton et al. 2005 (doi:10.2337/diacare.28.4.956)
struct NervePDE : public SimpleStructuralPDE {
  NervePDE(const SimParam* sp)
      : SimpleStructuralPDE(fields::kNerve, fields::kNerveId,
                            sp->neuropathy_diffusion, sp->neuropathy_decay),
        sp_(sp) {}

  void Init(Simulation* sim) override {
    SimpleStructuralPDE::Init(sim);
    // Set basal nerve density in dermal tissue
    auto* grid = Grid(sim);
    GridContext ctx(grid, sp_);
    real_t basal = sp_->neuropathy_basal_density;
    if (sp_->diabetic_mode) {
      basal *= sp_->diabetic_nerve_factor;
    }
    for (size_t idx = 0; idx < ctx.n; idx++) {
      real_t z = ctx.Z(idx);
      // Nerves in dermis only (papillary and reticular)
      if (z < 0 && z >= sp_->dermal_z_reticular) {
        grid->ChangeConcentrationBy(idx, basal);
      }
    }
  }

  // Wound event: denervation (Wallerian degeneration) in wound area
  void ApplyWound(Simulation* sim, real_t, real_t, real_t) override {
    auto* grid = Grid(sim);
    GridContext ctx(grid, sp_);
    for (size_t idx = 0; idx < ctx.n; idx++) {
      if (ctx.InWound(ctx.X(idx), ctx.Y(idx))) {
        real_t current = grid->GetConcentration(idx);
        if (current > 1e-10) {
          grid->ChangeConcentrationBy(idx, -current);
        }
      }
    }
  }

 private:
  const SimParam* sp_;
};

}  // namespace skibidy
}  // namespace bdm

#endif  // NERVE_PDE_H_
