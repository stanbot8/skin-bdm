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
  real_t amount = sp->immune_cytokine_rate * scale;
  if (amount < 1e-12) return;
  auto* rm = sim->GetResourceManager();
  if (sp->split_inflammation_enabled) {
    ScaledGrid sg(rm->GetDiffusionGrid(fields::kProInflammatory), sp);
    sg.AgentDeposit(sg.Index(qpos), amount);
  } else {
    ScaledGrid sg(rm->GetDiffusionGrid(fields::kInflammation), sp);
    if (sp->debug_immune) {
      real_t dx = qpos[0] - sp->wound_center_x;
      real_t dy = qpos[1] - sp->wound_center_y;
      std::cout << "[immune] infl qpos=(" << qpos[0] << "," << qpos[1]
                << "," << qpos[2] << ") dist=" << std::sqrt(dx*dx + dy*dy)
                << " r=" << sp->wound_radius
                << " scaled=" << amount * sg.agent_factor << std::endl;
    }
    sg.AgentDeposit(sg.Index(qpos), amount);
  }
  // Mirror deposit to immune pressure field (keratinocyte suppression signal)
  if (auto* ip_grid = rm->GetDiffusionGrid(fields::kImmunePressure)) {
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
  real_t eff_rate = sp->immune_resolution_rate;
  if (sp->diabetic_mode) {
    eff_rate *= sp->diabetic_resolution_factor;
  }

  if (sp->split_inflammation_enabled) {
    ScaledGrid sg(rm->GetDiffusionGrid(fields::kAntiInflammatory), sp);
    sg.AgentDeposit(sg.Index(qpos), eff_rate);
  } else {
    ScaledGrid sg(rm->GetDiffusionGrid(fields::kInflammation), sp);
    sg.AgentRemove(sg.Index(qpos), eff_rate);
  }
  // Mirror resolution to immune pressure field
  if (auto* ip_grid = rm->GetDiffusionGrid(fields::kImmunePressure)) {
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
  if (!sp->fibroblast_enabled) return;
  real_t amount = sp->m2_tgfb_rate * scale;
  if (amount < 1e-12) return;
  ScaledGrid sg(sim->GetResourceManager()->GetDiffusionGrid(fields::kTGFBeta), sp);
  sg.AgentDeposit(sg.Index(qpos), amount);
}

// Clear biofilm at the given grid position (immune cell phagocytosis).
// Neutrophils are more effective than macrophages (Jesaitis et al. 2003).
// Biofilm is non-diffusing -- no ScaledGrid.
inline void ClearBiofilm(const Real3& qpos, Simulation* sim,
                         const SimParam* sp, ImmuneCellType cell_type) {
  if (!sp->biofilm_enabled) return;
  auto* biofilm_grid =
      sim->GetResourceManager()->GetDiffusionGrid(fields::kBiofilm);
  size_t idx = biofilm_grid->GetBoxIndex(qpos);
  real_t current = biofilm_grid->GetConcentration(idx);
  if (current < 1e-10) return;

  real_t rate = (cell_type == kNeutrophil) ? sp->biofilm_neutrophil_clearance
                                           : sp->biofilm_macrophage_clearance;
  if (sp->diabetic_mode) {
    rate *= sp->diabetic_biofilm_clearance_factor;
  }

  real_t remove = std::min(current, rate);
  if (remove > 1e-10) {
    biofilm_grid->ChangeConcentrationBy(idx, -remove);
  }
}

// Write VEGF at the given grid position (M2 macrophage angiogenic signaling).
inline void ProduceVEGF(const Real3& qpos, Simulation* sim,
                        const SimParam* sp) {
  if (!sp->angiogenesis_enabled) return;
  ScaledGrid sg(sim->GetResourceManager()->GetDiffusionGrid(fields::kVEGF), sp);
  real_t rate = sp->m2_vegf_rate;
  if (sp->diabetic_mode) {
    rate *= sp->diabetic_vegf_factor;
  }
  sg.AgentDeposit(sg.Index(qpos), rate);
}

// Write MMP at the given grid position (M1 macrophage MMP-9 production).
// Lobmann et al. 2002: MMP-9 14-fold elevated in diabetic wounds.
inline void ProduceMMP(const Real3& qpos, Simulation* sim,
                       const SimParam* sp, real_t base_rate) {
  if (!sp->mmp_enabled) return;
  ScaledGrid sg(sim->GetResourceManager()->GetDiffusionGrid(fields::kMMP), sp);
  real_t rate = base_rate;
  if (sp->diabetic_mode) {
    rate *= sp->diabetic_mmp_factor;
  }
  sg.AgentDeposit(sg.Index(qpos), rate);
}

// M1 macrophage MMP-9 production (backward-compatible overload).
inline void ProduceMMP(const Real3& qpos, Simulation* sim,
                       const SimParam* sp) {
  ProduceMMP(qpos, sim, sp, sp->mmp_m1_rate);
}

// Chemotaxis with geometric fallback. Sets tractor force on cell.
inline void Migrate(ImmuneCell* cell, Simulation* sim, const SimParam* sp) {
  if (!sp->wound_enabled) return;

  Real3 pos = cell->GetPosition();
  Real3 qpos = ClampToBounds(pos, sim->GetParam());

  // Chemotaxis: follow inflammation gradient
  if (sp->chemotaxis_enabled) {
    auto* rm = sim->GetResourceManager();
    DiffusionGrid* grad_grid = nullptr;
    if (sp->split_inflammation_enabled) {
      grad_grid = rm->GetDiffusionGrid(fields::kProInflammatory);
    } else {
      grad_grid = rm->GetDiffusionGrid(fields::kInflammation);
    }

    Real3 gradient = {0, 0, 0};
    grad_grid->GetGradient(qpos, &gradient, false);
    real_t grad_mag = std::sqrt(gradient[0] * gradient[0] +
                                gradient[1] * gradient[1] +
                                gradient[2] * gradient[2]);

    if (grad_mag > 1e-8) {
      real_t speed =
          sp->immune_migration_speed * sp->chemotaxis_speed_scale * grad_mag;
      real_t inv_mag = 1.0 / grad_mag;
      cell->SetTractorForce(
          {gradient[0] * inv_mag * speed, gradient[1] * inv_mag * speed, 0});
      return;
    }
    // Flat gradient: fall through to geometric migration
  }

  // Geometric migration toward wound center
  real_t dx = pos[0] - sp->wound_center_x;
  real_t dy = pos[1] - sp->wound_center_y;
  real_t dist = std::sqrt(dx * dx + dy * dy);
  if (dist < 1e-6) return;

  real_t inv_dist = 1.0 / dist;
  real_t dir_x = -dx * inv_dist;
  real_t dir_y = -dy * inv_dist;

  real_t speed = sp->immune_migration_speed * (dist / sp->wound_radius);
  cell->SetTractorForce({dir_x * speed, dir_y * speed, 0});
}

}  // namespace immune
}  // namespace skibidy
}  // namespace bdm

#endif  // IMMUNE_HELPERS_H_
