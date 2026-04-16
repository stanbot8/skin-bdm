#ifndef FUSED_POST_H_
#define FUSED_POST_H_

#include "infra/util.h"
#include <tuple>

#include "core/hook_registry.h"
#include "core/hook_traits.h"
#include "core/operation/operation.h"
#include "core/wound_iteration.h"
#include "infra/sim_param.h"

namespace bdm {
namespace skibidy {

// FusedWoundPostOp: pure iteration skeleton. All domain-specific biology is
// delegated to self-contained module hooks registered in hook_registry.h.
// Ordering within each voxel preserves the dependency chain set by the
// PostHooks tuple order in hook_registry.h.
//
// Adding a new post module: create post_hook.h, add to PostHooks tuple, and
// (if needed) update the need_post check in registration.h.
struct FusedWoundPostOp : public StandaloneOperationImpl {
  BDM_OP_HEADER(FusedWoundPostOp);

  void operator()() override {
    auto* sim = Simulation::GetActive();
    auto* sp = sim->GetParam()->Get<SimParam>();
    PerfTimer timer(sp->debug_perf);

    ctx_.Prepare(sim, sp);

    // Pre-wound: nothing to do
    if (!ctx_.reg.PostWound()) return;

    PostHooks post_hooks;
    InitHooks(post_hooks, ctx_.reg, ctx_.sig);
    StudyPostHooks study_post_hooks;
    InitHooks(study_post_hooks, ctx_.reg, ctx_.sig);

    if (!AnyActive(post_hooks) && !AnyActive(study_post_hooks)) return;

    ctx_.IterateDermal([&](VoxelSnapshot& snap, SignalBoard& sig) {
      sig.ResetPerVoxel();
      std::apply([&](auto&... h) {
        (DispatchDermal(h, snap, sig), ...);
      }, post_hooks);
      std::apply([&](auto&... h) {
        (DispatchDermal(h, snap, sig), ...);
      }, study_post_hooks);
    });

    ctx_.IterateEpiWound([&](VoxelSnapshot& snap, SignalBoard& sig) {
      sig.ResetPerVoxel();
      std::apply([&](auto&... h) {
        (DispatchEpiWound(h, snap, sig), ...);
      }, post_hooks);
      std::apply([&](auto&... h) {
        (DispatchEpiWound(h, snap, sig), ...);
      }, study_post_hooks);
    });

    timer.Print("fused_post");
  }

 private:
  WoundIterationContext ctx_;
};

}  // namespace skibidy
}  // namespace bdm

#endif  // FUSED_POST_H_
