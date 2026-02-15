#ifndef CONTINUUM_HANDOFF_H_
#define CONTINUUM_HANDOFF_H_

// UWYN Continuum Handoff -- central dispatcher for agent-to-continuum
// state transfer. Call before RemoveFromSimulation() at every removal site.
//
// The UWYN paradigm: agents and continuum are the same cells in different
// representations. When an agent is removed, its state must be faithfully
// written to the continuum. An agent can also create NEW continuum types
// distinct from its surroundings (scar tissue, tumor mass).
//
// Each cell type branch documents what it writes, or explicitly documents
// why it is a no-op. Adding a new cell type means adding a branch here.

#include "tissue/keratinocyte.h"
#include "tumor/tumor_cell.h"
#include "fibroblast/fibroblast.h"
#include "immune/immune_cell.h"
#include "scar/scar_op.h"
#include "core/field_names.h"
#include "infra/sim_param.h"
#include "infra/util.h"

namespace bdm {
namespace skibidy {

inline void ContinuumHandoff(Agent* agent) {
  auto* sim = Simulation::GetActive();
  auto* sp = sim->GetParam()->Get<SimParam>();

  // --- Keratinocyte ---
  // Writes stratum value to continuum. Can create NEW continuum (scar)
  // when local collagen exceeds threshold (stratum+5 encoding).
  // Context-dependent handoffs (basal demotion, per-step stratum writes)
  // remain inline in differentiation.h. This is the fallback for any
  // removal path (safety timeout, coverage resolution).
  if (auto* kc = dynamic_cast<Keratinocyte*>(agent)) {
    if (sp->wound_enabled) {
      Real3 pos = kc->GetPosition();
      real_t dx = pos[0] - sp->wound_center_x;
      real_t dy = pos[1] - sp->wound_center_y;
      if (dx * dx + dy * dy <= sp->wound_radius * sp->wound_radius) {
        WriteScarValue(kc);
      }
    }
    return;
  }

  // --- TumorCell ---
  // Tumor mass persists after individual cell death. Writes binary marker
  // to kTumorId (new continuum type: tumor tissue distinct from normal
  // tissue) and stamps stratum for visualization.
  if (auto* tc = dynamic_cast<TumorCell*>(agent)) {
    if (!sp->tumor_enabled) return;
    auto* rm = sim->GetResourceManager();
    Real3 pos = ClampToBounds(tc->GetPosition(), sim->GetParam());

    auto* tumor_grid = rm->GetDiffusionGrid(fields::kTumorId);
    if (tumor_grid) {
      size_t tidx = tumor_grid->GetBoxIndex(pos);
      real_t t_cur = tumor_grid->GetConcentration(tidx);
      if (t_cur < 0.5) {
        tumor_grid->ChangeConcentrationBy(tidx, 1.0 - t_cur);
      }
    }

    if (sp->tumor_stratum_value > 0) {
      auto* stratum_grid = rm->GetDiffusionGrid(fields::kStratumId);
      if (stratum_grid) {
        size_t sidx = stratum_grid->GetBoxIndex(pos);
        real_t s_cur = stratum_grid->GetConcentration(sidx);
        if (sp->tumor_stratum_value > s_cur) {
          stratum_grid->ChangeConcentrationBy(
              sidx, sp->tumor_stratum_value - s_cur);
        }
      }
    }
    return;
  }

  // --- Fibroblast ---
  // Collagen, TGF-beta, MMP, fibronectin, elastin, and hyaluronan are
  // deposited every step during the fibroblast's life via fused_source.h.
  // No "fibroblast density" continuum exists. The ECM they deposited IS
  // their persistent continuum legacy. No additional writes needed.
  if (dynamic_cast<Fibroblast*>(agent)) {
    return;
  }

  // --- ImmuneCell (neutrophil + macrophage) ---
  // Cytokines (inflammation, anti-inflammation, immune pressure, VEGF,
  // MMP, NO) are deposited every step during the cell's life. No "immune
  // density" continuum exists. The cytokine fields they wrote to persist
  // and decay via their own PDE dynamics. No additional writes needed.
  if (dynamic_cast<ImmuneCell*>(agent)) {
    return;
  }
}

}  // namespace skibidy
}  // namespace bdm

#endif  // CONTINUUM_HANDOFF_H_
