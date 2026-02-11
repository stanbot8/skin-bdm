#ifndef TEMPERATURE_PDE_H_
#define TEMPERATURE_PDE_H_

#include "core/pde.h"
#include "infra/sim_param.h"

namespace bdm {
namespace skibidy {

// Wound temperature field: normalized to body core (37C = 1.0).
// Healthy skin ~37C (perfused). Open wound surface cools to ~33C (exposed).
// D>0: thermal conduction through tissue. No decay (energy conserved,
// sources/sinks handled in fused_source).
struct TemperaturePDE : public PDE {
  explicit TemperaturePDE(const SimParam* sp)
      : diffusion_(sp->temperature_diffusion),
        decay_(sp->temperature_decay) {}

  const char* GetName() const override { return fields::kTemperature; }
  int GetId() const override { return fields::kTemperatureId; }

  void Init(Simulation* sim) override {
    DefineGrid(sim, diffusion_, decay_);
    // Initialize everywhere to body temperature (normalized 1.0)
    auto* grid = Grid(sim);
    size_t n = grid->GetNumBoxes();
    for (size_t i = 0; i < n; i++) {
      grid->ChangeConcentrationBy(i, 1.0);
    }
  }

  // Wound disrupts thermal regulation: exposed surface cools.
  void ApplyWound(Simulation* sim, real_t cx, real_t cy, real_t r) override {
    auto* grid = Grid(sim);
    auto* sp = sim->GetParam()->Get<SimParam>();
    GridContext ctx(grid, sp);
    real_t wound_temp = sp->temperature_wound_surface;
    for (size_t idx = 0; idx < ctx.n; idx++) {
      if (!ctx.InWound(ctx.X(idx), ctx.Y(idx))) continue;
      // Epidermal wound voxels cool to surface temperature
      if (ctx.Z(idx) >= 0) {
        real_t current = grid->GetConcentration(idx);
        real_t delta = wound_temp - current;
        if (std::abs(delta) > 1e-10) {
          grid->ChangeConcentrationBy(idx, delta);
        }
      }
      // Dermal wound voxels partially cooled (depth-dependent)
      else {
        real_t depth_factor = std::min(1.0, -ctx.Z(idx) / 10.0);
        real_t target = wound_temp + (1.0 - wound_temp) * depth_factor;
        real_t current = grid->GetConcentration(idx);
        real_t delta = target - current;
        if (std::abs(delta) > 1e-10) {
          grid->ChangeConcentrationBy(idx, delta);
        }
      }
    }
  }

 private:
  real_t diffusion_;
  real_t decay_;
};

}  // namespace skibidy
}  // namespace bdm

#endif  // TEMPERATURE_PDE_H_
