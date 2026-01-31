#ifndef FIELDS_VASCULAR_H_
#define FIELDS_VASCULAR_H_

#include <cmath>

#include "core/pde.h"
#include "core/composite_field.h"
#include "infra/sim_param.h"

namespace bdm {
namespace skibidy {

// Vascular PDE channel.
// Represents vessel density/integrity in the dermis with layer-specific
// perfusion profiles: papillary (dense capillary loops, highest), reticular
// (arterioles/venules, moderate), hypodermis (sparse vessels, lowest).
// O2 and Water source terms read this field instead of hardcoding damage.
struct VascularPDE : public PDE {
  const char* GetName() const override { return fields::kVascular; }
  int GetId() const override { return fields::kVascularId; }

  void Init(Simulation* sim) override {
    auto* sp = sim->GetParam()->Get<SimParam>();
    DefineGrid(sim, sp->perfusion_diffusion, sp->perfusion_decay);

    // Layer-specific perfusion profile with smooth sigmoid transitions.
    // Sigmoid blending avoids artificial diffusion gradients at boundaries.
    real_t basal = sp->perfusion_basal;
    real_t papillary_val = basal * sp->perfusion_papillary_fraction;
    real_t reticular_val = basal * sp->perfusion_reticular_fraction;
    real_t hypodermis_val = basal * sp->perfusion_hypodermis_fraction;
    real_t z_papillary = sp->dermal_z_papillary;
    real_t z_reticular = sp->dermal_z_reticular;

    ModelInitializer::InitializeSubstance(GetId(),
        [=](real_t x, real_t y, real_t z) -> real_t {
          if (z > 0) return static_cast<real_t>(0);
          // Smooth sigmoid blend between layers (width 0.5 ~ 1 voxel)
          real_t w = 0.5;
          real_t f_papillary = 1.0 / (1.0 + std::exp(-(z - z_papillary) / w));
          real_t f_reticular = 1.0 / (1.0 + std::exp(-(z - z_reticular) / w));
          return f_papillary * papillary_val +
                 (1.0 - f_papillary) * (f_reticular * reticular_val +
                                         (1.0 - f_reticular) * hypodermis_val);
        });
  }

  // Wound disruption: zero vessels inside wound cylinder (dermis only).
  void ApplyWound(Simulation* sim, real_t, real_t, real_t) override {
    auto* grid = Grid(sim);
    auto* sp = sim->GetParam()->Get<SimParam>();
    GridContext ctx(grid, sp);
    GridContext::ZeroWoundDermis(grid, ctx);
  }

  // Angiogenesis source term: recover perfusion toward layer-specific baseline.
  // Gated by time delay (capillary sprouting onset) and tissue presence
  // above (Stratum field proxy for tissue demand / growth factor signals).
  void ApplySource(Simulation* sim, const CompositeField& fields) override {
    auto* sp = sim->GetParam()->Get<SimParam>();
    if (!sp->wound_enabled) return;

    auto* scheduler = sim->GetScheduler();
    uint64_t step = scheduler->GetSimulatedSteps();
    uint64_t wound_step = static_cast<uint64_t>(sp->wound_trigger_step);
    if (step <= wound_step) return;

    uint64_t wound_age = step - wound_step;
    if (wound_age < static_cast<uint64_t>(sp->perfusion_angio_delay)) return;

    auto* vasc_grid = Grid(sim);
    auto* stratum_grid = fields.Grid(fields::kStratum, sim);
    GridContext ctx(vasc_grid, sp);

    real_t basal = sp->perfusion_basal;
    real_t rate = sp->perfusion_angio_rate;
    real_t z_papillary = sp->dermal_z_papillary;
    real_t z_reticular = sp->dermal_z_reticular;
    real_t papillary_fraction = sp->perfusion_papillary_fraction;
    real_t reticular_fraction = sp->perfusion_reticular_fraction;
    real_t hypodermis_fraction = sp->perfusion_hypodermis_fraction;
    real_t angio_papillary = sp->angio_papillary_factor;
    real_t angio_reticular = sp->angio_reticular_factor;
    real_t angio_hypodermis = sp->angio_hypodermis_factor;

    // VEGF-modulated angiogenesis: cache grid pointer outside loop
    DiffusionGrid* vegf_grid = nullptr;
    if (sp->angiogenesis_enabled) {
      vegf_grid = fields.Grid(fields::kVEGF, sim);
    }

    for (size_t idx = 0; idx < ctx.n; idx++) {
      real_t z = ctx.Z(idx);
      if (z >= 0) continue;
      real_t x = ctx.X(idx), y = ctx.Y(idx);
      if (!ctx.InWound(x, y)) continue;

      // Per-layer target and rate
      real_t layer_fraction, layer_angio;
      if (z >= z_papillary) {
        layer_fraction = papillary_fraction;
        layer_angio = angio_papillary;
      } else if (z >= z_reticular) {
        layer_fraction = reticular_fraction;
        layer_angio = angio_reticular;
      } else {
        layer_fraction = hypodermis_fraction;
        layer_angio = angio_hypodermis;
      }
      real_t target = basal * layer_fraction;

      real_t current = vasc_grid->GetConcentration(idx);
      if (current >= target - 1e-10) continue;

      // Tissue demand gate: Stratum field at corresponding epidermal voxel
      // above this dermal position. More tissue above = stronger angio signal.
      // Floor of 0.2 represents wound-bed VEGF from platelets and inflammatory
      // cells, which drive basal angiogenesis even without epithelium above.
      Real3 above = {x, y, 1.0};
      real_t stratum_val = stratum_grid->GetValue(above);
      real_t demand = std::max(static_cast<real_t>(0.2),
                               GridContext::StratumGate(stratum_val));

      // Recovery: layer-scaled rate, VEGF-modulated when enabled
      real_t eff_rate = rate * layer_angio;
      if (vegf_grid) {
        real_t local_vegf = vegf_grid->GetConcentration(idx);
        eff_rate = sp->angio_vegf_rate * local_vegf * layer_angio;
      }

      real_t delta = eff_rate * demand;
      if (current + delta > target) {
        delta = target - current;
      }
      if (delta > 1e-10) {
        vasc_grid->ChangeConcentrationBy(idx, delta);

        // VEGF consumption: angiogenesis consumes VEGF proportional to recovery
        if (vegf_grid) {
          real_t vegf_current = vegf_grid->GetConcentration(idx);
          real_t consume = std::min(vegf_current,
                                    sp->vegf_consumption_rate * delta);
          if (consume > 1e-10) {
            vegf_grid->ChangeConcentrationBy(idx, -consume);
          }
        }
      }
    }
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // FIELDS_VASCULAR_H_
