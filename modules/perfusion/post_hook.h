#ifndef PERFUSION_POST_HOOK_H_
#define PERFUSION_POST_HOOK_H_

#include "core/grid_registry.h"
#include "core/voxel_snapshot.h"
#include "core/signal_board.h"

namespace bdm {
namespace skibidy {

// Post hook for venous perfusion clearance of TGF-beta.
// Blood flow through granulation tissue vasculature carries away
// soluble TGF-beta1 via venous drainage. Clearance scales with
// local perfusion (vascular density).
struct PerfusionPostHook {
  DiffusionGrid* tgfb_grid = nullptr;
  DiffusionGrid* perf_grid = nullptr;
  const SimParam* sp_ = nullptr;
  bool active = false;

  inline void Init(const GridRegistry& reg, SignalBoard& sig) {
    sp_ = reg.Params();
    active = sp_->perfusion.clearance_rate > 0;
    if (!active) return;
    tgfb_grid = reg.Get(fields::kTGFBetaId);
    perf_grid = reg.Get(fields::kVascularId);
    if (!tgfb_grid || !perf_grid) active = false;
  }

  inline void ApplyEpiWound(const VoxelSnapshot& snap, SignalBoard& sig) {
    real_t tgfb_val = tgfb_grid->GetConcentration(snap.idx);
    if (tgfb_val <= 1e-10) return;
    real_t perf = perf_grid->GetConcentration(snap.idx);
    if (perf <= 1e-10) return;
    real_t sink = sp_->perfusion.clearance_rate * perf * tgfb_val;
    tgfb_grid->ChangeConcentrationBy(snap.idx, -std::min(sink, tgfb_val));
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // PERFUSION_POST_HOOK_H_
