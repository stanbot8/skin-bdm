#ifndef SCAR_H_
#define SCAR_H_

#include <cmath>

#include "tissue/keratinocyte.h"
#include "core/field_names.h"
#include "core/pde.h"
#include "infra/sim_param.h"
#include "infra/util.h"
#include "core/operation/operation.h"

namespace bdm {
namespace skibidy {

// Emergent scar-vs-normal decision when a keratinocyte dissolves.
// Checks local collagen: high collagen -> scar tissue (stratum+5),
// low collagen -> normal stratum (healthy re-epithelialization).
// This makes scar formation an emergent outcome of the fibroblast/collagen
// cascade rather than an unconditional stamp.
inline void WriteScarValue(Keratinocyte* cell) {
  auto* sim = Simulation::GetActive();
  auto* sp = sim->GetParam()->Get<SimParam>();
  if (!sp->scar_enabled) return;

  auto* rm = sim->GetResourceManager();
  Real3 pos = ClampToBounds(cell->GetPosition(), sim->GetParam());

  // Determine if this voxel is scar or normal tissue.
  // When fibroblasts are enabled, local collagen drives the decision.
  // When disabled, all wound-area cells become scar (legacy behavior).
  bool is_scar = true;
  if (sp->fibroblast_enabled) {
    auto* col_grid = rm->GetDiffusionGrid(fields::kCollagen);
    // Check collagen at cell position AND at dermal depth below --
    // fibroblasts deposit collagen in the dermis with zero diffusion,
    // so epidermal cells need to sample the dermal layer.
    real_t collagen = col_grid->GetValue(pos);
    Real3 dermal_pos = {pos[0], pos[1], sp->dermal_fibroblast_depth};
    real_t dermal_col = col_grid->GetValue(dermal_pos);
    is_scar = (std::max(collagen, dermal_col) > sp->scar_collagen_threshold);
  }

  // Stratum field: scar voxels get stratum+5, normal get plain stratum
  auto* stratum_grid = rm->GetDiffusionGrid(fields::kStratum);
  size_t st_idx = stratum_grid->GetBoxIndex(pos);
  real_t new_val = static_cast<real_t>(cell->GetStratum());
  if (is_scar) new_val += 5.0;
  real_t st_cur = stratum_grid->GetConcentration(st_idx);
  if (new_val > st_cur) {
    stratum_grid->ChangeConcentrationBy(st_idx, new_val - st_cur);
  }

  // Scar field: binary marker for data/metrics
  // Skip when proportional mode is on (ScarAccumulationOp handles it)
  if (!sp->scar_proportional_enabled && is_scar) {
    auto* scar_grid = rm->GetDiffusionGrid(fields::kScar);
    size_t sc_idx = scar_grid->GetBoxIndex(pos);
    real_t sc_cur = scar_grid->GetConcentration(sc_idx);
    if (sc_cur < 0.5) {
      scar_grid->ChangeConcentrationBy(sc_idx, 1.0 - sc_cur);
    }
  }
}

// ---------------------------------------------------------------------------
// ScarAccumulationOp -- proportional scar accumulation.
// When fibroblast module is enabled, scar mirrors collagen density (the
// mechanistic scar output: myofibroblast-deposited ECM). When fibroblasts
// are disabled, falls back to inflammation integral (Ogawa 2017).
// ---------------------------------------------------------------------------
struct ScarAccumulationOp : public StandaloneOperationImpl {
  BDM_OP_HEADER(ScarAccumulationOp);

  void operator()() override {
    auto* sim = Simulation::GetActive();
    auto* sp = sim->GetParam()->Get<SimParam>();

    if (!sp->scar_proportional_enabled || !sp->wound_enabled) return;

    auto* rm = sim->GetResourceManager();
    auto* scar_grid = rm->GetDiffusionGrid(fields::kScar);

    // Collagen-driven scar (mechanistic): scar = collagen snapshot each step.
    // Collagen IS the scar -- deposited by myofibroblasts, degraded by MMPs.
    // Also stamps scar coloring into the Stratum field for visualization:
    // voxels where collagen > threshold get stratum 5 (scar basale).
    if (sp->fibroblast_enabled) {
      auto* col_grid = rm->GetDiffusionGrid(fields::kCollagen);
      auto* stratum_grid = rm->GetDiffusionGrid(fields::kStratum);
      GridContext ctx(scar_grid, sp);
      real_t z_max = sp->volume_z_cornified + ctx.box_len;
      real_t col_thresh = sp->scar_collagen_threshold;

      // Dermal z-level where fibroblasts deposit collagen
      real_t dermal_z = sp->dermal_fibroblast_depth;

      for (size_t idx = 0; idx < ctx.n; idx++) {
        real_t z = ctx.Z(idx);
        if (z < -5.0 || z > z_max) continue;
        real_t x = ctx.X(idx), y = ctx.Y(idx);
        if (!ctx.InWound(x, y)) continue;

        real_t collagen = col_grid->GetConcentration(idx);

        // Scar field: mirror collagen value
        real_t current = scar_grid->GetConcentration(idx);
        real_t delta = collagen - current;
        if (std::abs(delta) > 1e-10) {
          scar_grid->ChangeConcentrationBy(idx, delta);
        }

        // Stratum field: stamp scar coloring in epidermal wound voxels.
        // Collagen is deposited by dermal fibroblasts at z < 0 with no
        // diffusion, so for epidermal voxels we sample the dermal collagen
        // at the same (x,y) -- dermal ECM drives epidermal scar appearance.
        if (z >= 0) {
          real_t dermal_col = col_grid->GetValue({x, y, dermal_z});
          if (dermal_col > col_thresh) {
            real_t st_cur = stratum_grid->GetConcentration(idx);
            if (st_cur >= 1.0 && st_cur < 4.5) {
              // Normal healed stratum -> convert to scar variant (+5)
              real_t scar_val = std::round(st_cur) + 5.0;
              if (scar_val > st_cur) {
                stratum_grid->ChangeConcentrationBy(idx, scar_val - st_cur);
              }
            }
          }
        }
      }
      return;
    }

    // Fallback: inflammation integral when fibroblasts disabled
    DiffusionGrid* infl_grid = nullptr;
    if (sp->split_inflammation_enabled) {
      infl_grid = rm->GetDiffusionGrid(fields::kProInflammatory);
    } else {
      infl_grid = rm->GetDiffusionGrid(fields::kInflammation);
    }

    GridContext ctx(scar_grid, sp);
    real_t rate = sp->scar_accumulation_rate;
    real_t z_max = sp->volume_z_cornified + ctx.box_len;

    for (size_t idx = 0; idx < ctx.n; idx++) {
      real_t z = ctx.Z(idx);
      if (z < 0 || z > z_max) continue;
      real_t x = ctx.X(idx), y = ctx.Y(idx);
      if (!ctx.InWound(x, y)) continue;

      real_t inflammation = infl_grid->GetConcentration(idx);
      if (inflammation > 1e-10) {
        scar_grid->ChangeConcentrationBy(idx, inflammation * rate);
      }
    }
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // SCAR_H_
