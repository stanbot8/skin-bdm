#ifndef SCAR_POST_HOOK_H_
#define SCAR_POST_HOOK_H_

#include "core/grid_registry.h"
#include "core/voxel_snapshot.h"
#include "core/signal_board.h"

namespace bdm {
namespace skibidy {

// Post hook for scar module.
// Scar accumulation mirrors dermal collagen deposition (mechanistic) or
// falls back to inflammation integral (when fibroblasts disabled).
// Mechanotransduction amplification: stiff ECM triggers YAP/TAZ nuclear
// translocation, increasing alpha-SMA and collagen deposition.
// Aarabi et al. 2007 (doi:10.1096/fj.07-8218com)
struct ScarPostHook {
  DiffusionGrid* scar_grid = nullptr;
  DiffusionGrid* stiffness_grid = nullptr;
  DiffusionGrid* stratum_grid = nullptr;
  DiffusionGrid* col_grid = nullptr;
  DiffusionGrid* infl_grid = nullptr;
  const SimParam* sp_ = nullptr;
  bool active = false;
  bool do_scar_collagen = false;
  bool do_mechano = false;
  bool coarse_ = false;
  size_t dermal_z_offset = 0;
  size_t res2 = 0;
  std::vector<size_t> coarse_map_;

  inline void Init(const GridRegistry& reg, SignalBoard& sig) {
    sp_ = reg.Params();
    active = sp_->scar.proportional_enabled;
    if (!active) return;
    scar_grid = reg.Get(fields::kScarId);
    if (!scar_grid) { active = false; return; }

    stratum_grid = reg.Get(fields::kStratumId);

    // Compute grid geometry from stratum grid metadata
    auto* ref = reg.Get(fields::kStratumId);
    size_t res = ref->GetResolution();
    res2 = res * res;
    auto dims = ref->GetDimensions();
    real_t lo = static_cast<real_t>(dims[0]);
    real_t box_len = ref->GetBoxLength();

    do_scar_collagen = sp_->fibroblast.enabled;
    if (do_scar_collagen) {
      col_grid = reg.Get(fields::kCollagenId);
      // Precompute dermal z-slice offset for stratum stamping
      int dz_i = static_cast<int>(
          (sp_->dermal_fibroblast_depth - lo) / box_len);
      dz_i = std::max(0, std::min(dz_i, static_cast<int>(res) - 1));
      dermal_z_offset = static_cast<size_t>(dz_i) * res2;
      // Build coarse map for dermal_fine -> col_grid index mapping
      coarse_ = sp_->grid_resolution_structural > 0 &&
                sp_->grid_resolution_structural != sp_->grid_resolution;
      if (coarse_ && col_grid) {
        GridContext gctx(ref, sp_);
        coarse_map_.resize(gctx.n);
        for (size_t i = 0; i < gctx.n; i++) {
          Real3 pos = {gctx.X(i), gctx.Y(i), gctx.Z(i)};
          coarse_map_[i] = col_grid->GetBoxIndex(pos);
        }
      }
    } else {
      // Inflammation integral fallback
      infl_grid = reg.InflammationGrid();
    }

    do_mechano = sp_->mechanotransduction.enabled;
    if (do_mechano) {
      stiffness_grid = reg.Get(fields::kStiffnessId);
      if (!stiffness_grid) do_mechano = false;
    }
  }

  // Called per epidermal wound voxel. Fully self-contained.
  inline void ApplyEpiWound(const VoxelSnapshot& snap, SignalBoard& sig) {
    if (do_scar_collagen && col_grid) {
      // Mirror DERMAL collagen to epidermal scar field
      size_t dermal_fine = (snap.idx % res2) + dermal_z_offset;
      size_t dermal_ci = coarse_ ? coarse_map_[dermal_fine] : dermal_fine;
      real_t dermal_col = col_grid->GetConcentration(dermal_ci);
      size_t sc_si = snap.coarse_si;
      real_t sc_cur = scar_grid->GetConcentration(sc_si);
      real_t sc_delta = dermal_col - sc_cur;
      if (std::abs(sc_delta) > 1e-10) {
        scar_grid->ChangeConcentrationBy(sc_si, sc_delta);
      }
      // Stamp stratum +5 for scar visualization
      if (dermal_col > sp_->scar.collagen_threshold) {
        real_t st_cur = stratum_grid->GetConcentration(snap.idx);
        if (st_cur >= 1.0 && st_cur < 4.5) {
          real_t scar_val = std::round(st_cur) + 5.0;
          if (scar_val > st_cur) {
            stratum_grid->ChangeConcentrationBy(snap.idx, scar_val - st_cur);
          }
        }
      }
    } else if (infl_grid) {
      // Inflammation integral fallback (fibroblasts disabled)
      real_t inflammation = infl_grid->GetConcentration(snap.idx);
      if (inflammation > 1e-10) {
        scar_grid->ChangeConcentrationBy(
            snap.coarse_si,
            inflammation * sp_->scar.accumulation_rate * snap.coarse_w);
      }
    }

    // Mechanotransduction scar amplification
    if (do_mechano && stiffness_grid) {
      real_t stiff_val = stiffness_grid->GetConcentration(snap.idx);
      if (stiff_val > sp_->stiffness_yap_threshold) {
        size_t sc_si = snap.coarse_si;
        real_t sc_cur = scar_grid->GetConcentration(sc_si);
        if (sc_cur > 1e-10) {
          real_t amplify = sp_->stiffness_scar_factor *
              (stiff_val - sp_->stiffness_yap_threshold) * sc_cur;
          scar_grid->ChangeConcentrationBy(sc_si, amplify * snap.coarse_w);
        }
      }
    }
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // SCAR_POST_HOOK_H_
