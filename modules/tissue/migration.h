#ifndef MIGRATION_H_
#define MIGRATION_H_

#include <cmath>

#include "tissue/keratinocyte.h"
#include "core/field_names.h"
#include "infra/sim_param.h"
#include "infra/util.h"

namespace bdm {
namespace skibidy {

// ---------------------------------------------------------------------------
// Behavior: active cell migration (crawling) toward wound center.
// Basal keratinocytes extend lamellipodia on provisional wound matrix and
// crawl inward. Uses BioDynaMo's tractor force mechanism so migration
// composes naturally with mechanical repulsion (contact inhibition).
//
// Speed scales with distance from wound center: leader cells at the edge
// move fastest, cells near center decelerate as they converge.
// ---------------------------------------------------------------------------
struct Migration : public Behavior {
  BDM_BEHAVIOR_HEADER(Migration, Behavior, 1);

  Migration() { AlwaysCopyToNew(); }
  virtual ~Migration() {}

  void Run(Agent* agent) override {
    auto* cell = dynamic_cast<Keratinocyte*>(agent);
    if (!cell) return;

    // Once a cell differentiates, clear any residual migration force
    // so it can stratify vertically without being pulled inward.
    if (cell->GetStratum() != kBasal) {
      cell->SetTractorForce({0, 0, 0});
      return;
    }

    auto* sim = Simulation::GetActive();
    auto* sp = sim->GetParam()->Get<SimParam>();
    if (!sp->wound_enabled || !sp->migration_enabled) return;

    // Only migrate inside wound cylinder
    Real3 pos = cell->GetPosition();
    real_t dx = pos[0] - sp->wound_center_x;
    real_t dy = pos[1] - sp->wound_center_y;
    real_t dist2 = dx * dx + dy * dy;
    real_t r = sp->wound_radius;
    if (dist2 > r * r) return;

    real_t dist = std::sqrt(dist2);
    if (dist < 1e-6) return;  // at center, no direction

    // Direction: toward wound center
    real_t inv_dist = 1.0 / dist;
    real_t dir_x = -dx * inv_dist;
    real_t dir_y = -dy * inv_dist;

    // Speed: base speed, scaled by distance from center.
    // Floor at 30% prevents convergence stall near wound center --
    // leading-edge cells drive the sheet; followers maintain momentum.
    real_t speed = sp->migration_speed *
                   std::max(static_cast<real_t>(0.3), dist / r);

    // Moisture gating: cells cannot crawl on dry substrate
    auto* rm = sim->GetResourceManager();
    auto* water_grid = rm->GetDiffusionGrid(fields::kWater);
    Real3 qpos = ClampToBounds(pos, sim->GetParam());
    real_t water = water_grid->GetValue(qpos);
    real_t water_factor = std::min(1.0, water / sp->water_migration_threshold);
    if (water_factor < 0) water_factor = 0;
    speed *= water_factor;

    // ECM scaffold boost: precomputed composite (collagen + fibronectin +
    // elastin + fibrin) enhances keratinocyte crawling.
    if (g_ecm_quality) {
      real_t ecm = g_ecm_quality->GetValue(qpos);
      speed *= (1.0 + sp->ecm_migration_boost * ecm / (ecm + 1.0));
    }

    // Inflammation gating: Hill-function dose-response.
    // Read immune pressure (excludes wound DAMPs) to avoid feedback loop.
    // Normal: quadratic Hill K^2/(K^2+I^2) for sharp threshold.
    // Diabetic: linear Hill K/(K+I) for gradual, persistent suppression
    // reflecting chronic inflammatory hypersensitivity (Brownlee 2005).
    real_t infl = GetImmunePressure(sim, qpos);
    real_t eff_infl = infl;
    if (sp->diabetic_mode) {
      eff_infl *= sp->diabetic_inflammation_sensitivity;
    }
    real_t infl_factor = 1.0;
    if (sp->diabetic_mode) {
      real_t K = sp->diabetic_migration_infl_K;
      if (K + eff_infl > 1e-12) {
        infl_factor = K / (K + eff_infl);
      }
    } else {
      real_t K = sp->inflammation_migration_threshold;
      real_t K2 = K * K;
      real_t denom = K2 + eff_infl * eff_infl;
      if (denom > 1e-12) {
        infl_factor = K2 / denom;
      }
    }
    speed *= infl_factor;

    // pH gating: alkaline wound bed suppresses keratinocyte migration.
    // Schneider et al. 2007: acidic pH (5.5-6.5) promotes crawling.
    {
      auto* ph_grid = rm->GetDiffusionGrid(fields::kPH);
      if (ph_grid) {
        real_t alkalinity = ph_grid->GetValue(qpos);
        if (alkalinity > 0) {
          speed *= 1.0 - sp->ph_migration_suppression *
                          std::min(static_cast<real_t>(1.0), alkalinity);
        }
      }
    }

    // Diabetic mode: impaired keratinocyte crawling
    if (sp->diabetic_mode) {
      speed *= sp->diabetic_migration_factor;
    }

    cell->SetTractorForce({dir_x * speed, dir_y * speed, 0});
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // MIGRATION_H_
