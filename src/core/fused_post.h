#ifndef FUSED_POST_H_
#define FUSED_POST_H_

#include "infra/util.h"
#include <cmath>
#include <tuple>
#include <vector>

#include "core/field_names.h"
#include "core/pde.h"
#include "core/grid_registry.h"
#include "core/voxel_snapshot.h"
#include "core/signal_board.h"
#include "core/hook_traits.h"
#include "core/hook_registry.h"
#include "infra/sim_param.h"
#include "core/operation/operation.h"

namespace bdm {
namespace skibidy {

// ---------------------------------------------------------------------------
// FusedWoundPostOp -- pure iteration skeleton.
// All domain-specific biology is delegated to self-contained module hooks
// registered in hook_registry.h. Ordering within each voxel preserves the
// dependency chain (set by PostHooks tuple order in hook_registry.h).
//
// Adding a new post module: create post_hook.h, add to PostHooks tuple,
// and (if needed) update the need_post check in registration.h.
// ---------------------------------------------------------------------------
struct FusedWoundPostOp : public StandaloneOperationImpl {
  BDM_OP_HEADER(FusedWoundPostOp);

  void operator()() override {
    auto* sim = Simulation::GetActive();
    auto* sp = sim->GetParam()->Get<SimParam>();
    PerfTimer timer(sp->debug_perf);

    // Pre-wound: nothing to do
    uint64_t step = GetGlobalStep(sim);
    uint64_t wound_step = static_cast<uint64_t>(sp->wound.trigger_step);
    if (step <= wound_step) return;

    // --- Registry + SignalBoard ---
    GridRegistry reg;
    reg.Fill(sim);
    SignalBoard sig;

    // --- Hook construction + Init (engine + study) ---
    PostHooks post_hooks;
    InitHooks(post_hooks, reg, sig);
    StudyPostHooks study_post_hooks;
    InitHooks(study_post_hooks, reg, sig);

    // Early exit if nothing is active
    if (!AnyActive(post_hooks) && !AnyActive(study_post_hooks)) return;

    // --- Grid context + wound mask (lazy init) ---
    GridContext ctx(reg.Get(fields::kStratumId), sp);
    if (mask_.mask.empty()) {
      real_t z_max = sp->volume_z_cornified + ctx.box_len;
      mask_ = GridContext::ComputeWoundMask(ctx, z_max, sp);
    }

    // --- Coarse map (lazy init) ---
    bool coarse = sp->grid_resolution_structural > 0 &&
                  sp->grid_resolution_structural != sp->grid_resolution;
    real_t coarse_w = 1.0;
    if (coarse) {
      real_t rf = static_cast<real_t>(sp->grid_resolution);
      real_t rc = static_cast<real_t>(sp->grid_resolution_structural);
      coarse_w = (rc * rc * rc) / (rf * rf * rf);
      if (coarse_map_.empty()) {
        DiffusionGrid* ref = reg.Get(fields::kCollagenId);
        if (!ref) ref = reg.Get(fields::kFibronectinId);
        if (!ref) ref = reg.Get(fields::kScarId);
        if (ref) {
          coarse_map_.resize(ctx.n);
          for (size_t i = 0; i < ctx.n; i++) {
            Real3 pos = {ctx.X(i), ctx.Y(i), ctx.Z(i)};
            coarse_map_[i] = ref->GetBoxIndex(pos);
          }
        }
      }
    }

    // --- SnapshotFiller ---
    SnapshotFiller filler;
    filler.Init(reg, ctx, mask_, coarse_map_, coarse, coarse_w);

    // --- Dermal wound loop (erosion, tissue effects) ---
    VoxelSnapshot snap;
    for (size_t idx : mask_.dermal_all) {
      filler.FillDermal(snap, idx);
      sig.ResetPerVoxel();
      std::apply([&](auto&... h) {
        (DispatchDermal(h, snap, sig), ...);
      }, post_hooks);
      std::apply([&](auto&... h) {
        (DispatchDermal(h, snap, sig), ...);
      }, study_post_hooks);
    }

    // --- Epidermal wound loop ---
    for (size_t idx : mask_.epi_wound) {
      filler.FillEpi(snap, idx);
      sig.ResetPerVoxel();
      std::apply([&](auto&... h) {
        (DispatchEpiWound(h, snap, sig), ...);
      }, post_hooks);
      std::apply([&](auto&... h) {
        (DispatchEpiWound(h, snap, sig), ...);
      }, study_post_hooks);
    }

    timer.Print("fused_post");
  }

 private:
  GridContext::WoundMaskData mask_;
  std::vector<size_t> coarse_map_;
};

}  // namespace skibidy
}  // namespace bdm

#endif  // FUSED_POST_H_
