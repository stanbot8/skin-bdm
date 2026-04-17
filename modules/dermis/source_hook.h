#ifndef DERMIS_SOURCE_HOOK_H_
#define DERMIS_SOURCE_HOOK_H_

#include "core/hook_api.h"

namespace bdm {
namespace skibidy {

// Source hook for dermis module.
// Dermis tissue recovery driven by collagen deposition (collagen is the
// structural scaffold; dermis density follows). MMP degrades dermis.
// Damaged tissue emits DAMPs (damage-associated molecular patterns)
// that sustain wound inflammation.
struct DermisSourceHook {
  DiffusionGrid* dermis_grid = nullptr;
  DiffusionGrid* col_grid = nullptr;
  DiffusionGrid* mmp_grid = nullptr;
  DiffusionGrid* infl_grid = nullptr;
  const SimParam* sp_ = nullptr;
  bool active = false;
  bool do_wound_infl = false;
  real_t wound_infl_rate = 0;

  inline void Init(const GridRegistry& reg, SignalBoard& sig) {
    sp_ = reg.Params();
    active = sp_->dermis.enabled;
    if (!active) return;
    dermis_grid = reg.Get(fields::kDermisId);
    if (!dermis_grid) { active = false; return; }
    if (sp_->fibroblast.enabled)
      col_grid = reg.Get(fields::kCollagenId);
    if (sp_->mmp.enabled)
      mmp_grid = reg.Get(fields::kMMPId);

    wound_infl_rate = reg.WoundInflRate();
    do_wound_infl = wound_infl_rate > 0;
    if (do_wound_infl) {
      infl_grid = reg.InflammationGrid();
      if (!infl_grid) do_wound_infl = false;
    }
  }

  // Dermal: tissue recovery + MMP degradation + DAMP source
  inline void ApplyDermal(const VoxelSnapshot& snap, SignalBoard& sig) {
    if (!snap.in_wound || !snap.post_wound) return;

    size_t d_idx = snap.coarse_si;
    real_t dermis_cur = dermis_grid->GetConcentration(d_idx);

    // Per-layer target and rate
    real_t d_target, d_factor;
    if (snap.z >= sp_->dermal_z_papillary) {
      d_target = sp_->dermis.papillary_density;
      d_factor = sp_->dermis.papillary_rate_factor;
    } else if (snap.z >= sp_->dermal_z_reticular) {
      d_target = sp_->dermis.reticular_density;
      d_factor = sp_->dermis.reticular_rate_factor;
    } else {
      d_target = sp_->dermis.hypodermis_density;
      d_factor = sp_->dermis.hypodermis_rate_factor;
    }

    // Collagen-driven recovery (collagen is structural)
    if (dermis_cur < d_target && col_grid) {
      real_t col_val = col_grid->GetConcentration(snap.coarse_si);
      if (col_val > sp_->dermis.collagen_threshold) {
        real_t gain = sp_->dermis.collagen_recovery_rate * d_factor *
                      (col_val - sp_->dermis.collagen_threshold);
        if (dermis_cur + gain > d_target) gain = d_target - dermis_cur;
        if (gain > 1e-10) {
          dermis_grid->ChangeConcentrationBy(d_idx, gain * snap.coarse_w);
          dermis_cur += gain * snap.coarse_w;
        }
      }
    }

    // MMP degradation (MMP is fine grid, dermis is structural)
    if (mmp_grid && dermis_cur > 0) {
      real_t mmp_val = mmp_grid->GetConcentration(snap.idx);
      if (mmp_val > 1e-10) {
        real_t loss = std::min(dermis_cur,
                               mmp_val * sp_->dermis.mmp_degradation);
        if (loss > 1e-10) {
          dermis_grid->ChangeConcentrationBy(d_idx, -loss * snap.coarse_w);
          dermis_cur -= loss * snap.coarse_w;
        }
      }
    }

    // Dermal DAMP source: inflammation from damaged tissue (fine grid)
    if (do_wound_infl && dermis_cur < d_target) {
      real_t damage_frac = 1.0 - dermis_cur / d_target;
      real_t infl_gain = wound_infl_rate * damage_frac;
      if (infl_gain > 1e-10) {
        infl_grid->ChangeConcentrationBy(snap.idx, infl_gain);
      }
    }
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // DERMIS_SOURCE_HOOK_H_
