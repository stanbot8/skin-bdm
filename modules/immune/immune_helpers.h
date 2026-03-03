#ifndef IMMUNE_HELPERS_H_
#define IMMUNE_HELPERS_H_

#include <cmath>

#include "immune/immune_cell.h"
#include "core/field_names.h"
#include "core/pde.h"
#include "infra/sim_param.h"
#include "infra/util.h"

namespace bdm {
namespace skibidy {
namespace immune {

// Write pro-inflammatory cytokines at the given grid position.
// Used by neutrophils and M1 macrophages.
// scale: age-based taper (1.0 = full output, decays as cell ages).
// This makes the inflammation peak-and-decay shape emerge from cell
// aging dynamics rather than requiring precise resolution_rate tuning.
inline void ProduceProInflammatory(const Real3& qpos, Simulation* sim,
                                   const SimParam* sp,
                                   real_t scale = 1.0) {
  real_t amount = sp->immune.cytokine_rate * scale;
  if (amount < 1e-12) return;
  auto* rm = sim->GetResourceManager();
  if (sp->inflammation.split_inflammation_enabled) {
    ScaledGrid sg(rm->GetDiffusionGrid(fields::kProInflammatoryId), sp);
    sg.AgentDeposit(sg.Index(qpos), amount);
  } else {
    ScaledGrid sg(rm->GetDiffusionGrid(fields::kInflammationId), sp);
    if (sp->debug_immune) {
      real_t dx = qpos[0] - sp->wound.center_x;
      real_t dy = qpos[1] - sp->wound.center_y;
      std::cout << "[immune] infl qpos=(" << qpos[0] << "," << qpos[1]
                << "," << qpos[2] << ") dist=" << std::sqrt(dx*dx + dy*dy)
                << " r=" << sp->wound.radius
                << " scaled=" << amount * sg.agent_factor << std::endl;
    }
    sg.AgentDeposit(sg.Index(qpos), amount);
  }
  // Mirror deposit to immune pressure field (keratinocyte suppression signal)
  if (auto* ip_grid = rm->GetDiffusionGrid(fields::kImmunePressureId)) {
    ScaledGrid sg(ip_grid, sp);
    sg.AgentDeposit(sg.Index(qpos), amount);
  }
}

// Write anti-inflammatory mediators (M2 macrophages).
// Split mode: add to anti-inflammatory grid.
// Single mode: subtract from inflammation grid.
inline void ProduceAntiInflammatory(const Real3& qpos, Simulation* sim,
                                    const SimParam* sp) {
  auto* rm = sim->GetResourceManager();
  real_t eff_rate = sp->immune.resolution_rate;
  if (sp->diabetic.mode) {
    eff_rate *= sp->diabetic.resolution_factor;
  }

  if (sp->inflammation.split_inflammation_enabled) {
    ScaledGrid sg(rm->GetDiffusionGrid(fields::kAntiInflammatoryId), sp);
    sg.AgentDeposit(sg.Index(qpos), eff_rate);
  } else {
    ScaledGrid sg(rm->GetDiffusionGrid(fields::kInflammationId), sp);
    sg.AgentRemove(sg.Index(qpos), eff_rate);
  }
  // Mirror resolution to immune pressure field
  if (auto* ip_grid = rm->GetDiffusionGrid(fields::kImmunePressureId)) {
    ScaledGrid sg(ip_grid, sp);
    sg.AgentRemove(sg.Index(qpos), eff_rate);
  }
}

// Write TGF-beta at the given grid position (fibroblast module coupling).
// scale: age-based taper (1.0 = full output, decays as M2 ages).
// M2 TGF-b output declines as resolution progresses, matching
// literature TGF-b peak day 5-7 followed by decline.
inline void ProduceTGFBeta(const Real3& qpos, Simulation* sim,
                           const SimParam* sp, real_t scale = 1.0) {
  if (!sp->fibroblast.enabled) return;
  real_t amount = sp->m2_tgfb_rate * scale;
  if (amount < 1e-12) return;
  ScaledGrid sg(sim->GetResourceManager()->GetDiffusionGrid(fields::kTGFBetaId), sp);
  sg.AgentDeposit(sg.Index(qpos), amount);
}

// Receptor-mediated TGF-beta endocytosis. Cells with TGF-beta receptors
// (TbRII/TbRI) internalize ligand proportional to local concentration.
// Clearance scales with cell density, replacing the abstract decay constant.
// (Vilar et al. 2006, doi:10.1016/j.jtbi.2006.03.024;
//  Wakefield et al. 1990, doi:10.1172/JCI114647)
inline void ConsumeTGFBeta(const Real3& qpos, Simulation* sim,
                           const SimParam* sp, real_t local_tgfb) {
  if (!sp->fibroblast.enabled) return;
  if (sp->fibroblast.tgfb_receptor_consumption <= 0) return;
  real_t consumed = sp->fibroblast.tgfb_receptor_consumption * local_tgfb;
  if (consumed < 1e-12) return;
  ScaledGrid sg(sim->GetResourceManager()->GetDiffusionGrid(fields::kTGFBetaId), sp);
  sg.AgentRemove(sg.Index(qpos), consumed);
}

// Clear biofilm at the given grid position (immune cell phagocytosis).
// Neutrophils are more effective than macrophages (Jesaitis et al. 2003).
// Biofilm is non-diffusing -- no ScaledGrid.
inline void ClearBiofilm(const Real3& qpos, Simulation* sim,
                         const SimParam* sp, ImmuneCellType cell_type) {
  if (!sp->biofilm.enabled) return;
  auto* biofilm_grid =
      sim->GetResourceManager()->GetDiffusionGrid(fields::kBiofilmId);
  if (!biofilm_grid) return;
  size_t idx = biofilm_grid->GetBoxIndex(qpos);
  real_t current = biofilm_grid->GetConcentration(idx);
  if (current < 1e-10) return;

  real_t rate = (cell_type == kNeutrophil) ? sp->biofilm.neutrophil_clearance
                                           : sp->biofilm.macrophage_clearance;
  if (sp->diabetic.mode) {
    rate *= sp->diabetic.biofilm_clearance_factor;
  }

  real_t remove = std::min(current, rate);
  if (remove > 1e-10) {
    biofilm_grid->ChangeConcentrationBy(idx, -remove);
  }
}

// Write VEGF at the given grid position (M2 macrophage angiogenic signaling).
inline void ProduceVEGF(const Real3& qpos, Simulation* sim,
                        const SimParam* sp) {
  if (!sp->angiogenesis.enabled) return;
  ScaledGrid sg(sim->GetResourceManager()->GetDiffusionGrid(fields::kVEGFId), sp);
  real_t rate;
  if (sp->mech_vegf_production) {
    // Mechanistic: HIF-1alpha stabilization under hypoxia. O2 below
    // threshold stabilizes HIF-1a, which transactivates VEGF promoter.
    // Rate = max_rate * max(0, 1 - O2/threshold)
    auto* o2_grid = sim->GetResourceManager()->GetDiffusionGrid(fields::kOxygenId);
    real_t o2 = o2_grid ? std::max(static_cast<real_t>(0), o2_grid->GetValue(qpos))
                        : 1.0;
    real_t hypoxia_signal = std::max(static_cast<real_t>(0),
        1.0 - o2 / sp->mech_hif_o2_threshold);
    rate = sp->mech_hif_vegf_rate * hypoxia_signal;
  } else {
    // Parametric: flat rate per M2 macrophage
    rate = sp->m2_vegf_rate;
  }
  if (sp->diabetic.mode) {
    rate *= sp->diabetic.vegf_factor;
  }
  sg.AgentDeposit(sg.Index(qpos), rate);
}

// Write pro-MMP (zymogen) at the given grid position. All cellular MMP
// production deposits into the inactive pro-MMP pool; activation to active
// MMP occurs in fused_post.h via plasmin and autocatalytic cleavage.
// Visse & Nagase 2003 (doi:10.1161/01.res.0000070112.80711.3d)
inline void ProduceMMP(const Real3& qpos, Simulation* sim,
                       const SimParam* sp, real_t base_rate) {
  if (!sp->mmp.enabled) return;
  auto* prommp_grid = sim->GetResourceManager()->GetDiffusionGrid(fields::kProMMPId);
  if (!prommp_grid) return;
  ScaledGrid sg(prommp_grid, sp);
  real_t rate = base_rate;
  if (sp->diabetic.mode) {
    rate *= sp->diabetic.mmp_factor;
  }
  sg.AgentDeposit(sg.Index(qpos), rate);
}

// M1 macrophage MMP-9 production (backward-compatible overload).
inline void ProduceMMP(const Real3& qpos, Simulation* sim,
                       const SimParam* sp) {
  ProduceMMP(qpos, sim, sp, sp->mmp.m1_rate);
}

// Write TIMP at the given grid position. TIMP-1/2 production from fibroblasts,
// M2 macrophages, and keratinocytes. Diabetic tissue has reduced TIMP expression.
// Brew et al. 2000 (doi:10.1016/S0167-4838(99)00252-5)
inline void ProduceTIMP(const Real3& qpos, Simulation* sim,
                        const SimParam* sp, real_t base_rate) {
  if (!sp->mmp.enabled) return;
  auto* timp_grid = sim->GetResourceManager()->GetDiffusionGrid(fields::kTIMPId);
  if (!timp_grid) return;
  ScaledGrid sg(timp_grid, sp);
  real_t rate = base_rate;
  if (sp->diabetic.mode) {
    rate *= sp->diabetic.timp_production_factor;
  }
  sg.AgentDeposit(sg.Index(qpos), rate);
}

// Chemotaxis with geometric fallback. Sets tractor force on cell.
inline void Migrate(ImmuneCell* cell, Simulation* sim, const SimParam* sp) {
  if (!sp->wound.enabled) return;

  Real3 pos = cell->GetPosition();
  Real3 qpos = ClampToBounds(pos, sim->GetParam());

  // Chemotaxis: follow inflammation gradient
  if (sp->chemotaxis_enabled) {
    auto* rm = sim->GetResourceManager();
    DiffusionGrid* grad_grid = nullptr;
    if (sp->inflammation.split_inflammation_enabled) {
      grad_grid = rm->GetDiffusionGrid(fields::kProInflammatoryId);
    } else {
      grad_grid = rm->GetDiffusionGrid(fields::kInflammationId);
    }

    Real3 gradient = {0, 0, 0};
    grad_grid->GetGradient(qpos, &gradient, false);
    real_t grad_mag = std::sqrt(gradient[0] * gradient[0] +
                                gradient[1] * gradient[1] +
                                gradient[2] * gradient[2]);

    if (grad_mag > 1e-8) {
      real_t speed =
          sp->immune.migration_speed * sp->chemotaxis_speed_scale * grad_mag;
      real_t inv_mag = 1.0 / grad_mag;
      cell->SetTractorForce(
          {gradient[0] * inv_mag * speed, gradient[1] * inv_mag * speed, 0});
      return;
    }
    // Flat gradient: fall through to geometric migration
  }

  // Geometric migration toward wound center
  real_t dx = pos[0] - sp->wound.center_x;
  real_t dy = pos[1] - sp->wound.center_y;
  real_t dist = std::sqrt(dx * dx + dy * dy);
  if (dist < 1e-6) return;

  real_t inv_dist = 1.0 / dist;
  real_t dir_x = -dx * inv_dist;
  real_t dir_y = -dy * inv_dist;

  real_t speed = sp->immune.migration_speed * (dist / sp->wound.radius);
  cell->SetTractorForce({dir_x * speed, dir_y * speed, 0});
}

}  // namespace immune
}  // namespace skibidy
}  // namespace bdm

#endif  // IMMUNE_HELPERS_H_
