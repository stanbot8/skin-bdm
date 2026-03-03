#ifndef LYMPHATIC_SOURCE_HOOK_H_
#define LYMPHATIC_SOURCE_HOOK_H_

#include "core/grid_registry.h"
#include "core/voxel_snapshot.h"
#include "core/signal_board.h"

namespace bdm {
namespace skibidy {

// Source hook for lymphatic module.
// Lymphatic endothelial cells sprout from intact margins via VEGF-C/D
// signaling through VEGFR-3. Slower than angiogenesis (~7 day lag).
// Also handles edema source (inflammatory vascular leak) and sink
// (lymphatic drainage).
// Kataru et al. 2009 (doi:10.1182/blood-2008-09-176776)
struct LymphaticSourceHook {
  DiffusionGrid* lymph_grid = nullptr;
  DiffusionGrid* edema_grid = nullptr;
  DiffusionGrid* vegf_grid = nullptr;
  DiffusionGrid* infl_grid = nullptr;
  DiffusionGrid* o2_grid = nullptr;
  const SimParam* sp_ = nullptr;
  bool active = false;

  inline void Init(const GridRegistry& reg, SignalBoard& sig) {
    sp_ = reg.Params();
    active = sp_->lymphatic.enabled;
    if (!active) return;
    lymph_grid = reg.Get(fields::kLymphaticId);
    edema_grid = reg.Get(fields::kEdemaId);
    if (!lymph_grid || !edema_grid) { active = false; return; }
    o2_grid = reg.Get(fields::kOxygenId);
    if (sp_->angiogenesis.enabled)
      vegf_grid = reg.Get(fields::kVEGFId);
    infl_grid = reg.InflammationGrid();
  }

  // Dermal: lymphatic regen + edema source/sink
  inline void ApplyDermal(const VoxelSnapshot& snap, SignalBoard& sig) {
    if (!snap.in_wound || !snap.post_wound) return;

    // Lymphatic regeneration
    size_t l_si = snap.coarse_si;
    real_t lymph_val = lymph_grid->GetConcentration(l_si);
    real_t lymph_target = sp_->lymphatic.basal_density;
    if (sp_->diabetic.mode) {
      lymph_target *= sp_->diabetic.lymphatic_factor;
    }
    if (lymph_val < lymph_target - 1e-10) {
      real_t rate = sp_->lymphatic.regen_rate;
      if (sp_->diabetic.mode) {
        rate *= sp_->diabetic.lymphatic_factor;
      }
      // VEGF promotes lymphangiogenesis
      if (vegf_grid) {
        real_t local_vegf = vegf_grid->GetConcentration(snap.idx);
        rate *= (1.0 + sp_->lymphatic.vegf_boost * local_vegf);
      }
      real_t gain = std::min(lymph_target - lymph_val, rate);
      if (gain > 1e-10) {
        lymph_grid->ChangeConcentrationBy(l_si, gain * snap.coarse_w);
      }
    }

    // Edema source: inflammatory vascular leak (Starling forces)
    real_t infl_val = infl_grid ? infl_grid->GetConcentration(snap.idx) : 0;
    real_t leak = sp_->edema_leak_rate * infl_val * snap.vasc;
    if (leak > 1e-10) {
      edema_grid->ChangeConcentrationBy(snap.idx, leak);
    }

    // Edema sink: lymphatic drainage
    real_t edema_cur = edema_grid->GetConcentration(snap.idx);
    if (edema_cur > 1e-10 && lymph_val > 1e-10) {
      real_t drain = std::min(edema_cur,
          sp_->edema_drainage_rate * lymph_val * edema_cur);
      edema_grid->ChangeConcentrationBy(snap.idx, -drain);
    }
  }

  // Epidermal wound: edema O2 impairment
  inline void ApplyEpiWound(const VoxelSnapshot& snap, SignalBoard& sig) {
    real_t edema_val = edema_grid->GetConcentration(snap.idx);
    if (edema_val > 1e-10) {
      real_t o2_cur = o2_grid->GetConcentration(snap.idx);
      if (o2_cur > 1e-10) {
        real_t impair = std::min(o2_cur,
            sp_->edema_o2_impairment * edema_val * o2_cur);
        o2_grid->ChangeConcentrationBy(snap.idx, -impair);
      }
    }
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // LYMPHATIC_SOURCE_HOOK_H_
