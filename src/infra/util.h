#ifndef INFRA_UTIL_H_
#define INFRA_UTIL_H_

#include <cmath>

#include "biodynamo.h"
#include "core/derived_field.h"
#include "core/field_names.h"
#include "infra/sim_param.h"

namespace bdm {
namespace skibidy {

// Step offset for checkpoint-fork runs. When loading a checkpoint saved at
// step N, g_step_offset is set to N so that all step-dependent logic
// (immune spawn, wound event, metrics export) uses the correct global time.
inline uint64_t g_step_offset = 0;

// Get the global simulation step (local step + checkpoint offset).
inline uint64_t GetGlobalStep(Simulation* sim) {
  return sim->GetScheduler()->GetSimulatedSteps() + g_step_offset;
}

// Global derived field pointers, set once during simulation setup (skibidy.h).
// Same pattern as BDM's Simulation::GetActive() global accessor.
inline DerivedField* g_ecm_quality = nullptr;
inline DerivedField* g_tissue_viability = nullptr;
inline DerivedField* g_wound_microenv = nullptr;

// Derived field helpers (return 0 when fields not enabled).
inline real_t GetECMQuality(const Real3& qpos) {
  return g_ecm_quality ? g_ecm_quality->GetValue(qpos) : 0;
}
inline real_t GetTissueViability(const Real3& qpos) {
  return g_tissue_viability ? g_tissue_viability->GetValue(qpos) : 0;
}

// Clamp a position inside the simulation bounds to avoid out-of-bounds
// lookups. Uses per-axis bounds from SimParam when available, otherwise
// falls back to Param scalar bounds.
inline Real3 ClampToBounds(const Real3& pos, const Param* p) {
  constexpr real_t margin = 1e-3;
  auto* sp = p->Get<SimParam>();
  return {std::max(sp->bounds_min[0] + margin,
                   std::min(pos[0], sp->bounds_max[0] - margin)),
          std::max(sp->bounds_min[1] + margin,
                   std::min(pos[1], sp->bounds_max[1] - margin)),
          std::max(sp->bounds_min[2] + margin,
                   std::min(pos[2], sp->bounds_max[2] - margin))};
}

// Read net inflammation at a position, handling split vs single field.
// Returns max(0, pro - weight*anti) when split, or raw value otherwise.
inline real_t GetNetInflammation(Simulation* sim, const Real3& qpos) {
  auto* sp = sim->GetParam()->Get<SimParam>();
  auto* rm = sim->GetResourceManager();
  if (sp->inflammation.split_inflammation_enabled) {
    auto* pro_grid = rm->GetDiffusionGrid(fields::kProInflammatoryId);
    auto* anti_grid = rm->GetDiffusionGrid(fields::kAntiInflammatoryId);
    real_t pro = pro_grid->GetValue(qpos);
    real_t anti = anti_grid->GetValue(qpos);
    return std::max(real_t{0},
                    pro - sp->anti_inflammation_weight * anti);
  } else {
    auto* infl_grid = rm->GetDiffusionGrid(fields::kInflammationId);
    return infl_grid->GetValue(qpos);
  }
}

// Immune pressure at a position -- same dynamics as inflammation but excludes
// wound-derived DAMPs. Used by keratinocytes for proliferation/migration
// suppression to avoid feedback with wound source.
inline real_t GetImmunePressure(Simulation* sim, const Real3& qpos) {
  auto* rm = sim->GetResourceManager();
  auto* grid = rm->GetDiffusionGrid(fields::kImmunePressureId);
  if (!grid) return 0;
  return grid->GetValue(qpos);
}

// Q10 temperature scaling factor. temp_normalized is the grid value
// (fraction of 37C), q10 is the scaling coefficient per 10C.
// Returns the multiplicative rate adjustment.
inline real_t Q10Factor(real_t temp_normalized, real_t q10) {
  real_t temp_c = temp_normalized * 37.0;
  return std::pow(q10, (temp_c - 37.0) / 10.0);
}

}  // namespace skibidy
}  // namespace bdm

#endif  // INFRA_UTIL_H_
