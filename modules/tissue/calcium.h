#ifndef FIELDS_CALCIUM_H_
#define FIELDS_CALCIUM_H_

#include "core/pde.h"
#include "core/composite_field.h"
#include "infra/sim_param.h"

namespace bdm {
namespace skibidy {

// Calcium PDE channel.
// Extracellular Ca2+ drives differentiation (Grabe & Neuber 2005).
// Maintained by tight junctions and ion pumps - not free diffusion.
// Source term restores the calcium gradient inside the wound cylinder
// after wound-driven disruption, gated on Stratum recovery.
struct CalciumPDE : public PDE {
  const char* GetName() const override { return fields::kCalcium; }
  int GetId() const override { return fields::kCalciumId; }

  static real_t Profile(const SimParam* sp, real_t z) {
    return GridContext::Sigmoid(z, sp->calcium_basal, sp->calcium_peak,
                               sp->calcium_midpoint_z, sp->calcium_steepness);
  }

  void Init(Simulation* sim) override {
    auto* sp = sim->GetParam()->Get<SimParam>();
    DefineGrid(sim, sp->calcium_diffusion, sp->calcium_decay);
    ModelInitializer::InitializeSubstance(GetId(),
        [sp](real_t x, real_t y, real_t z) { return Profile(sp, z); });
    MarkPrescribed(sim);  // ion-pump-maintained, not free diffusion
  }

  void ApplyWound(Simulation* sim, real_t, real_t, real_t) override {
    ZeroInWound(sim);
  }

  // Calcium recovery: tight junctions reform as cells re-populate the wound,
  // gradually re-establishing the extracellular calcium gradient.
  // Gated on Stratum recovery (cells must be present to rebuild tight junctions).
  void ApplySource(Simulation* sim, const CompositeField& fields) override {
    auto* sp = sim->GetParam()->Get<SimParam>();
    if (!sp->wound_enabled) return;

    auto* ca_grid = Grid(sim);
    auto* stratum_grid = fields.Grid(fields::kStratum, sim);
    GridContext ctx(ca_grid, sp);
    real_t rate = sp->calcium_recovery_rate;

    for (size_t idx = 0; idx < ctx.n; idx++) {
      real_t x = ctx.X(idx), y = ctx.Y(idx);
      if (!ctx.InWound(x, y)) continue;

      real_t stratum_val = stratum_grid->GetConcentration(idx);
      if (stratum_val < 0.5) continue;

      real_t barrier = GridContext::StratumGate(stratum_val);
      real_t target = Profile(sp, ctx.Z(idx));
      real_t current = ca_grid->GetConcentration(idx);
      real_t delta = (target - current) * rate * barrier;
      if (std::abs(delta) > 1e-10) {
        ca_grid->ChangeConcentrationBy(idx, delta);
      }
    }
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // FIELDS_CALCIUM_H_
