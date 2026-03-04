#ifndef BASAL_DENSITY_PDE_H_
#define BASAL_DENSITY_PDE_H_

#include <cmath>

#include "core/pde.h"
#include "core/composite_field.h"
#include "core/field_names.h"
#include "core/operation/operation.h"
#include "infra/sim_param.h"
#include "infra/util.h"

namespace bdm {
namespace skibidy {

// ---------------------------------------------------------------------------
// BasalDensityPDE -- continuum basal keratinocyte population field.
//
// Represents the normalized density of homeostatic basal cells OUTSIDE the
// active wound zone. This is the central idea of the extended UWYN paradigm:
// continuum is not just cells at rest -- it is any spatially homogeneous
// population whose individual variation doesn't affect macro-scale outcomes.
//
// Paradigm rules:
//   1. Homeostatic basale (far from wound): represented here as density ρ ∈ [0,1].
//   2. Wound-front basale (within demotion_radius): explicit agents.
//   3. Agents and continuum co-exist -- an agent IS a local deviation above
//      the density floor, not a replacement for it.
//   4. Demotion: agents that drift outside demotion_radius fold back into ρ.
//
// D=0: basal cells do not diffuse laterally (they are anchored to the basement
// membrane by integrin alpha6-beta4). Independent tissue columns.
// decay=0: the population self-maintains at carrying capacity under homeostasis.
// ---------------------------------------------------------------------------
struct BasalDensityPDE : public PDE {
  const char* GetName() const override { return fields::kBasalDensity; }
  int GetId() const override { return fields::kBasalDensityId; }

  void Init(Simulation* sim) override {
    auto* sp = sim->GetParam()->Get<SimParam>();
    // D=0 (no lateral spread), decay=0 (homeostatic self-maintenance)
    DefineGrid(sim, 0.0, 0.0);

    // Initial condition: homeostatic density = 1.0 everywhere in the tissue
    // basal layer (x/y inside tissue domain, z in basale zone [0, spinous]).
    // This represents the healthy, intact basale before any wound.
    real_t spinous_z = sp->volume_z_spinous;
    real_t tmin = sp->tissue_min;
    real_t tmax = sp->tissue_max;
    ModelInitializer::InitializeSubstance(GetId(),
        [tmin, tmax, spinous_z](real_t x, real_t y, real_t z) -> real_t {
          if (x < tmin || x > tmax || y < tmin || y > tmax) return 0.0;
          if (z < 0.0 || z > spinous_z) return 0.0;
          return 1.0;  // homeostatic density
        });
  }

  // Wound creation: zero out basal density inside the punch biopsy cylinder.
  // Agents at the wound margin take over from this point.
  void ApplyWound(Simulation* sim, real_t cx, real_t cy, real_t r) override {
    ZeroInWound(sim);
  }
};

// ---------------------------------------------------------------------------
// BasalDensityOp -- logistic PDE evolution for homeostatic basale.
//
// Evolution equation per voxel (outside wound zone, z in basale):
//
//   dρ/dt = r_div(O₂, KGF, water) · ρ · (1 - ρ/ρ_max)   ← logistic growth
//          - r_diff · ρ                                    ← differentiation loss
//
// Where r_div is the intrinsic division rate modulated by the same
// resource signals as BasalDivision (so continuum and agents see the same
// biology). Under homeostasis (O₂=1, KGF=2, water=1), ρ stays at ρ_max.
// In a wound (ρ→0 at wound event), ρ grows back as O₂ recovers.
// The density field is NOT used by agents -- it is a parallel representation
// that costs zero agent dispatches.
// ---------------------------------------------------------------------------
struct BasalDensityOp : public StandaloneOperationImpl {
  BDM_OP_HEADER(BasalDensityOp);

  void operator()() override {
    auto* sim = Simulation::GetActive();
    auto* sp = sim->GetParam()->Get<SimParam>();
    if (!sp->basal_density_enabled || !sp->wound.enabled) return;
    PerfTimer timer(sp->debug_perf);

    // Sub-cycling: density changes slowly (logistic timescale >> dt).
    uint64_t step = GetGlobalStep(sim);
    if (sp->basal_density_subcycle > 1 &&
        step % static_cast<uint64_t>(sp->basal_density_subcycle) != 0) return;

    auto* rm = sim->GetResourceManager();
    auto* dens_grid = rm->GetDiffusionGrid(fields::kBasalDensityId);
    auto* o2_grid   = rm->GetDiffusionGrid(fields::kOxygenId);
    auto* kgf_grid  = rm->GetDiffusionGrid(fields::kKGFId);
    auto* water_grid = rm->GetDiffusionGrid(fields::kWaterId);
    if (!dens_grid || !o2_grid || !kgf_grid || !water_grid) return;

    GridContext ctx(dens_grid, sp);
    real_t dt_eff = sim->GetParam()->simulation_time_step *
                    static_cast<real_t>(sp->basal_density_subcycle);
    real_t rho_max    = sp->basal_density_max;
    real_t r_div_base = sp->basal_density_div_rate;
    real_t r_diff     = sp->basal_density_diff_rate;
    real_t spinous_z  = sp->volume_z_spinous;
    real_t tmin       = sp->tissue_min;
    real_t tmax       = sp->tissue_max;
    real_t kgf_km     = sp->kgf_half_maximal;
    real_t kgf_boost  = sp->kgf_max_boost;
    real_t o2_thresh  = sp->oxygen_prolif_threshold;
    real_t w_thresh   = sp->water_prolif_threshold;

    for (size_t idx = 0; idx < ctx.n; idx++) {
      real_t z = ctx.Z(idx);
      // Only evolve voxels in the tissue basale zone
      if (z < 0.0 || z > spinous_z) continue;
      real_t x = ctx.X(idx), y = ctx.Y(idx);
      if (x < tmin || x > tmax || y < tmin || y > tmax) continue;
      // Skip wound zone -- agents handle this
      if (ctx.InWound(x, y)) continue;

      real_t rho = dens_grid->GetConcentration(idx);
      if (rho <= 0.0 && rho_max <= 0.0) continue;

      // Resource-modulated division rate (same biology as BasalDivision G1)
      Real3 pos = {x, y, z};
      real_t kgf   = kgf_grid->GetValue(pos);
      real_t o2    = o2_grid->GetValue(pos);
      real_t water = water_grid->GetValue(pos);

      real_t kgf_factor = 1.0 + kgf_boost * kgf / (kgf_km + kgf + 1e-12);
      real_t o2_factor  = std::min(1.0, o2 / (o2_thresh + 1e-12));
      real_t w_factor   = std::min(1.0, water / (w_thresh + 1e-12));
      if (o2_factor < 0) o2_factor = 0;
      if (w_factor < 0) w_factor = 0;

      real_t r_div = r_div_base * kgf_factor * o2_factor * w_factor;

      // Logistic growth with differentiation outflux
      real_t d_rho = dt_eff * (r_div * rho * (1.0 - rho / rho_max) -
                                r_diff * rho);

      if (std::abs(d_rho) > 1e-12) {
        real_t new_rho = std::max(0.0, rho + d_rho);
        dens_grid->ChangeConcentrationBy(idx, new_rho - rho);
      }
    }
    timer.Print("basal_density");
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // BASAL_DENSITY_PDE_H_
