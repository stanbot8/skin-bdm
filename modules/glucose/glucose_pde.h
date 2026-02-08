#ifndef GLUCOSE_PDE_H_
#define GLUCOSE_PDE_H_

#include "core/pde.h"
#include "infra/sim_param.h"

namespace bdm {
namespace skibidy {

// Interstitial glucose field: normalized to healthy fasting level (1.0).
// Supplied by vascular perfusion, consumed by cells and bacteria.
// In diabetic mode, basal level elevated (hyperglycemia) driving AGE formation.
struct GlucosePDE : public PDE {
  explicit GlucosePDE(const SimParam* sp)
      : diffusion_(sp->glucose_diffusion),
        decay_(sp->glucose_decay) {}

  const char* GetName() const override { return fields::kGlucose; }
  int GetId() const override { return fields::kGlucoseId; }

  void Init(Simulation* sim) override {
    DefineGrid(sim, diffusion_, decay_);
    // Initialize everywhere to basal glucose level
    auto* sp = sim->GetParam()->Get<SimParam>();
    auto* grid = Grid(sim);
    size_t n = grid->GetNumBoxes();
    for (size_t i = 0; i < n; i++) {
      grid->ChangeConcentrationBy(i, sp->glucose_basal_conc);
    }
  }

  // Wound disrupts vasculature: glucose supply drops in wound bed
  void ApplyWound(Simulation* sim, real_t cx, real_t cy, real_t r) override {
    auto* grid = Grid(sim);
    auto* sp = sim->GetParam()->Get<SimParam>();
    GridContext ctx(grid, sp);
    for (size_t idx = 0; idx < ctx.n; idx++) {
      if (!ctx.InWound(ctx.X(idx), ctx.Y(idx))) continue;
      // Dermal wound voxels lose glucose supply proportional to vascular damage
      real_t current = grid->GetConcentration(idx);
      real_t target = sp->glucose_basal_conc * 0.5;  // 50% of basal
      if (current > target) {
        grid->ChangeConcentrationBy(idx, target - current);
      }
    }
  }

 private:
  real_t diffusion_;
  real_t decay_;
};

}  // namespace skibidy
}  // namespace bdm

#endif  // GLUCOSE_PDE_H_
