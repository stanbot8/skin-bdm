#ifndef LYMPHATIC_POST_HOOK_H_
#define LYMPHATIC_POST_HOOK_H_

#include "core/hook_api.h"

namespace bdm {
namespace skibidy {

// Post hook for lymphatic TGF-beta drainage.
// Local lymphatic vessel density determines drainage capacity.
// Kataru et al. 2009 (doi:10.1182/blood-2008-09-176776)
struct LymphaticPostHook {
  DiffusionGrid* tgfb_grid = nullptr;
  DiffusionGrid* lymph_grid = nullptr;
  const SimParam* sp_ = nullptr;
  bool active = false;

  inline void Init(const GridRegistry& reg, SignalBoard& sig) {
    sp_ = reg.Params();
    active = sp_->lymphatic.enabled;
    if (!active) return;
    tgfb_grid = reg.Get(fields::kTGFBetaId);
    lymph_grid = reg.Get(fields::kLymphaticId);
    if (!tgfb_grid || !lymph_grid) active = false;
  }

  inline void ApplyEpiWound(const VoxelSnapshot& snap, SignalBoard& sig) {
    real_t tgfb_val = tgfb_grid->GetConcentration(snap.idx);
    if (tgfb_val <= 1e-10) return;
    real_t lymph_val = lymph_grid->GetConcentration(snap.coarse_si);
    if (lymph_val <= 1e-10) return;
    real_t sink = sp_->edema_drainage_rate * lymph_val * tgfb_val;
    tgfb_grid->ChangeConcentrationBy(snap.idx, -std::min(sink, tgfb_val));
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // LYMPHATIC_POST_HOOK_H_
