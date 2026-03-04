#ifndef MIGRATION_H_
#define MIGRATION_H_

#include <cmath>

#include "tissue/keratinocyte.h"
#include "core/field_names.h"
#include "core/voxel_env.h"
#include "infra/sim_param.h"
#include "infra/util.h"

namespace bdm {
namespace skibidy {

// ---------------------------------------------------------------------------
// Behavior: mechanistic keratinocyte migration toward wound center.
//
// Direction integrates three orthogonal biophysical cues (weighted sum,
// then normalized):
//   1. O2 gradient (hypoxia chemotaxis, weight=1.0, always active):
//      HIF-1alpha activation drives keratinocytes toward the hypoxic wound
//      bed. Primary and most robust signal -- wound vasculature disruption
//      creates a steep, wound-geometry-encoding O2 gradient regardless of
//      wound size (Simon 2004, Tandara & Mustoe 2004).
//   2. VEGF chemotaxis (weight=vegf_migration_weight, when field active):
//      Keratinocyte VEGFR-1 (Flt-1) mediates directed migration toward VEGF
//      secreted by hypoxic macrophages and fibroblasts in the wound bed
//      (Bauer et al. 2005, Rossiter et al. 2004). Reinforces O2 direction.
//   3. Fibronectin haptotaxis (weight=fibronectin_haptotaxis_weight, when
//      field active): integrin alpha-5-beta-1 engages fibronectin in the
//      provisional matrix; cells crawl up the FN gradient deposited by the
//      platelet clot and activated fibroblasts (Grinnell 1984, Clark 1996).
//
// Speed: driven by O2 gradient magnitude (steepest at wound edge -> fastest
//   leader cells; flat near wound center -> natural deceleration).
//   migration_speed calibrated to keratinocyte crawling at ~30 um/h
//   (Larsen et al. 2006): 30 um/h x 0.1h/step / 5 um/unit = 0.6 units/step
//   net; migration_speed = 3.0 provides the tractor force amplitude that
//   yields this net displacement against contact mechanics.
//
// Contact inhibition of locomotion (CIL):
//   When the leading-edge voxel (one cell diameter ahead) already has
//   differentiated tissue (stratum >= 1 = spinosum or higher), Eph-ephrin
//   signaling suppresses forward tractor force. This provides the natural
//   closure stop signal -- wound coverage plateaus when cells can no longer
//   advance into already-healed tissue (Coburn et al. 2013).
//
// Additional gating from wound-bed microenvironment:
//   moisture, ECM scaffold boost, inflammation (normal: quadratic Hill;
//   diabetic: linear Hill for chronic hypersensitivity), pH, temperature Q10.
//   See individual comments for literature references.
// ---------------------------------------------------------------------------
struct Migration : public Behavior {
  BDM_BEHAVIOR_HEADER(Migration, Behavior, 1);

  Migration() { AlwaysCopyToNew(); }
  virtual ~Migration() {}

  // Per-cell gradient cache: recomputed every migration_subcycle steps.
  // Caching normalized direction + O2 speed factor eliminates 3 GetGradient
  // calls on skip steps. Microenvironment scalars (water, inflammation, pH,
  // temperature) are cheap GetValue calls and are still evaluated every step.
  real_t cached_dir_x_ = 0;
  real_t cached_dir_y_ = 0;
  real_t cached_speed_factor_ = 0.3;
  uint64_t last_grad_step_ = UINT64_MAX;

  void Run(Agent* agent) override {
    auto* cell = dynamic_cast<Keratinocyte*>(agent);
    if (!cell) return;

    // Once a cell differentiates, clear residual migration force so it
    // stratifies vertically without being pulled inward.
    if (cell->GetStratum() != kBasal) {
      cell->SetTractorForce({0, 0, 0});
      return;
    }

    auto* sim = Simulation::GetActive();
    auto* sp = sim->GetParam()->Get<SimParam>();
    if (!sp->wound.enabled || !sp->migration_enabled) return;

    // Only migrate inside wound cylinder
    Real3 pos = cell->GetPosition();
    real_t dx = pos[0] - sp->wound.center_x;
    real_t dy = pos[1] - sp->wound.center_y;
    real_t dist2 = dx * dx + dy * dy;
    real_t r = sp->wound.radius;
    if (dist2 > r * r) return;

    real_t dist = std::sqrt(dist2);
    if (dist < 1e-6) return;  // at center, no direction

    Real3 qpos = ClampToBounds(pos, sim->GetParam());
    size_t vi = g_voxel_env.Index(qpos);
    const auto& env = g_voxel_env[vi];

    // ------------------------------------------------------------------
    // Gradient cues with temporal sub-cycling.
    // O2/VEGF/FN gradients change slowly (PDE time scale >> 0.1h/step),
    // so we recompute direction and speed factor only every
    // migration_subcycle steps and reuse the cached values otherwise.
    // Scalar reads come from the fused VoxelEnv (no redundant pos-to-idx).
    // ------------------------------------------------------------------
    uint64_t cur_step = GetGlobalStep(sim);
    bool need_grad = (last_grad_step_ == UINT64_MAX ||
        cur_step - last_grad_step_ >=
            static_cast<uint64_t>(sp->migration_subcycle));

    real_t dir_x, dir_y, speed_factor;
    if (need_grad) {
      auto* rm = sim->GetResourceManager();

      // Cue 1: O2 gradient -- hypoxia chemotaxis (primary, always active).
      auto* o2_grid = rm->GetDiffusionGrid(fields::kOxygenId);
      Real3 o2_grad = {0, 0, 0};
      if (o2_grid) o2_grid->GetGradient(qpos, &o2_grad, false);
      real_t gx = -o2_grad[0];
      real_t gy = -o2_grad[1];
      real_t gmag_o2 = std::sqrt(gx * gx + gy * gy);

      // Cue 2: VEGF chemotaxis.
      real_t vx = 0, vy = 0;
      {
        auto* vegf_grid = rm->GetDiffusionGrid(fields::kVEGFId);
        if (vegf_grid) {
          Real3 vegf_grad = {0, 0, 0};
          vegf_grid->GetGradient(qpos, &vegf_grad, false);
          vx = vegf_grad[0];
          vy = vegf_grad[1];
        }
      }

      // Cue 3: Fibronectin haptotaxis.
      real_t fx = 0, fy = 0;
      {
        auto* fn_grid = rm->GetDiffusionGrid(fields::kFibronectinId);
        if (fn_grid) {
          Real3 fn_grad = {0, 0, 0};
          fn_grid->GetGradient(qpos, &fn_grad, false);
          fx = fn_grad[0];
          fy = fn_grad[1];
        }
      }

      // Direction: weighted combination, then normalize.
      real_t w_vegf = sp->angiogenesis.vegf_migration_weight;
      real_t w_fn   = sp->fibronectin.haptotaxis_weight;
      real_t cx_sum = gx + w_vegf * vx + w_fn * fx;
      real_t cy_sum = gy + w_vegf * vy + w_fn * fy;
      real_t gmag_sum = std::sqrt(cx_sum * cx_sum + cy_sum * cy_sum);

      if (gmag_sum > 1e-8) {
        dir_x = cx_sum / gmag_sum;
        dir_y = cy_sum / gmag_sum;
        if (gmag_o2 > 1e-8) {
          speed_factor = std::max(static_cast<real_t>(0.3),
                                  std::min(static_cast<real_t>(1.0),
                                           gmag_o2 * sp->migration_gradient_scale));
        } else {
          speed_factor = static_cast<real_t>(0.3);
        }
      } else {
        // All gradients flat: geometric fallback toward wound center
        real_t inv_dist = 1.0 / dist;
        dir_x = -dx * inv_dist;
        dir_y = -dy * inv_dist;
        speed_factor = static_cast<real_t>(0.3);
      }

      cached_dir_x_ = dir_x;
      cached_dir_y_ = dir_y;
      cached_speed_factor_ = speed_factor;
      last_grad_step_ = cur_step;
    } else {
      dir_x = cached_dir_x_;
      dir_y = cached_dir_y_;
      speed_factor = cached_speed_factor_;
    }

    real_t speed = sp->migration_speed * speed_factor;

    // ------------------------------------------------------------------
    // Contact inhibition of locomotion (CIL): Eph/ephrin stop signal.
    // Sample the stratum field one cell diameter ahead in the migration
    // direction. If that voxel already has differentiated epithelium
    // (stratum >= 1.0 = spinosum or scar-encoded), suppress tractor force.
    // ------------------------------------------------------------------
    {
      Real3 fwd = {pos[0] + dir_x * sp->division_diameter,
                   pos[1] + dir_y * sp->division_diameter,
                   pos[2]};
      fwd = ClampToBounds(fwd, sim->GetParam());
      size_t fwd_vi = g_voxel_env.Index(fwd);
      if (g_voxel_env[fwd_vi].stratum >= 1.0) {
        speed *= sp->cil_suppression;
      }
    }

    // Moisture gating: cells cannot crawl on dry substrate
    real_t water_factor = std::min(1.0,
        env.water / sp->water_migration_threshold);
    if (water_factor < 0) water_factor = 0;
    speed *= water_factor;

    // ECM scaffold boost: precomputed composite (collagen + fibronectin +
    // elastin + fibrin) enhances keratinocyte crawling on ECM substrate.
    {
      real_t ecm = env.ecm_quality;
      if (ecm > 0) {
        speed *= (1.0 + sp->ecm_migration_boost * ecm / (ecm + 1.0));
      }
    }

    // Inflammation gating: Hill-function dose-response.
    // Normal: quadratic Hill K^2/(K^2+I^2) for sharp threshold.
    // Diabetic: linear Hill K/(K+I) for gradual, persistent suppression
    // reflecting chronic inflammatory hypersensitivity (Brownlee 2005).
    real_t eff_infl = env.immune_pressure;
    if (sp->diabetic.mode) {
      eff_infl *= sp->diabetic.inflammation_sensitivity;
    }
    real_t infl_factor = 1.0;
    if (sp->diabetic.mode) {
      real_t K = sp->diabetic.migration_infl_K;
      if (K + eff_infl > 1e-12) {
        infl_factor = K / (K + eff_infl);
      }
    } else {
      real_t K = sp->inflammation.migration_threshold;
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
      real_t alkalinity = env.ph;
      if (alkalinity > 0) {
        speed *= 1.0 - sp->ph.migration_suppression *
                        std::min(static_cast<real_t>(1.0), alkalinity);
      }
    }

    // Temperature Q10: wound cooling reduces migration speed
    // Q10 ~ 2.0 for keratinocyte crawling (Kanokwan & Bhattacharya 2012)
    if (sp->temperature.enabled &&
        std::abs(sp->temperature.q10_migration - 1.0) > 1e-6) {
      real_t temp_val = env.temperature;
      real_t temp_c = temp_val * 37.0;  // denormalize to Celsius
      speed *= std::pow(sp->temperature.q10_migration,
                        (temp_c - 37.0) / 10.0);
    }

    // Diabetic mode: impaired keratinocyte crawling
    if (sp->diabetic.mode) {
      speed *= sp->diabetic.migration_factor;
    }

    cell->SetTractorForce({dir_x * speed, dir_y * speed, 0});
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // MIGRATION_H_
