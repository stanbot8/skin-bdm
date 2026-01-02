#ifndef FIELDS_WATER_H_
#define FIELDS_WATER_H_

#include "core/pde.h"
#include "core/composite_field.h"
#include "infra/sim_param.h"

namespace bdm {
namespace skibidy {

// Water PDE channel.
// Models tissue moisture content. Healthy skin maintains a hydration gradient
// from dermis (wet) to corneum (dry, barrier function). After wounding, the
// barrier is breached and tissue water drops due to trans-epidermal water loss
// (TEWL). Dermal vasculature supplies serum/exudate that re-hydrates the wound
// bed over hours, gating keratinocyte migration and proliferation.
struct WaterPDE : public PDE {
  const char* GetName() const override { return fields::kWater; }
  int GetId() const override { return fields::kWaterId; }

  void Init(Simulation* sim) override {
    auto* sp = sim->GetParam()->Get<SimParam>();
    DefineGrid(sim, sp->water_diffusion, sp->water_decay);
    ModelInitializer::InitializeSubstance(GetId(),
        [sp](real_t x, real_t y, real_t z) {
          if (z < 0) return sp->water_basal_conc;
          return GridContext::ExpDecay(z, sp->water_basal_conc,
                                      sp->water_decay_length);
        });
  }

  void ApplyWound(Simulation* sim, real_t cx, real_t cy, real_t r) override {
    ZeroInWound(sim);
  }

  // Source term: dermal supply + surface evaporation.
  // Healthy dermis: hard-pin to vascular perfusion (instant equilibrium).
  // Wound dermis: let diffusion fill from intact margins (exudate), with
  // gentle recovery toward vascular target as angiogenesis restores vessels.
  // Epidermal wound: evaporation + serum recovery from below.
  void ApplySource(Simulation* sim, const CompositeField& fields) override {
    auto* sp = sim->GetParam()->Get<SimParam>();
    auto* water_grid = Grid(sim);
    auto* vasc_grid = fields.Grid(fields::kVascular, sim);
    auto* stratum_grid = fields.Grid(fields::kStratum, sim);
    GridContext ctx(water_grid, sp);

    // Dermal supply: healthy tissue hard-pins; wound lets diffusion fill
    for (size_t idx = 0; idx < ctx.n; idx++) {
      if (ctx.Z(idx) >= 0) continue;
      real_t target = sp->water_basal_conc * vasc_grid->GetConcentration(idx);
      real_t current = water_grid->GetConcentration(idx);

      if (sp->wound_enabled && ctx.InWound(ctx.X(idx), ctx.Y(idx))) {
        // Wound dermis: gentle recovery toward vascular target.
        // Diffusion from margins provides the rapid initial fill (exudate).
        if (current < target) {
          real_t gain = std::min(target - current, sp->water_recovery_rate);
          water_grid->ChangeConcentrationBy(idx, gain);
        }
      } else {
        // Healthy dermis: hard pin to vascular perfusion
        real_t delta = target - current;
        if (std::abs(delta) > 1e-10) {
          water_grid->ChangeConcentrationBy(idx, delta);
        }
      }
    }

    // Epidermal: evaporation + serum recovery (wound voxels only)
    if (!sp->wound_enabled) return;
    for (size_t idx = 0; idx < ctx.n; idx++) {
      real_t z = ctx.Z(idx);
      if (z < 0) continue;

      real_t x = ctx.X(idx), y = ctx.Y(idx);
      if (!ctx.InWound(x, y)) continue;

      // Barrier integrity from Stratum field: high stratum = intact barrier
      real_t barrier = GridContext::StratumGate(
          stratum_grid->GetConcentration(idx));

      // Evaporation: exposed surface loses water, barrier reduces loss
      real_t evap = sp->water_surface_loss_rate * (1.0 - barrier);
      real_t current = water_grid->GetConcentration(idx);
      if (evap > 0 && current > 0) {
        real_t loss = std::min(current, evap);
        water_grid->ChangeConcentrationBy(idx, -loss);
      }

      // Serum recovery from below: push water up toward target
      real_t target = GridContext::ExpDecay(z, sp->water_basal_conc,
                                           sp->water_decay_length);
      current = water_grid->GetConcentration(idx);
      if (current < target) {
        real_t gain = std::min(target - current, sp->water_recovery_rate);
        water_grid->ChangeConcentrationBy(idx, gain);
      }
    }
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // FIELDS_WATER_H_
