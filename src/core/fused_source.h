#ifndef FUSED_SOURCE_H_
#define FUSED_SOURCE_H_

#include "infra/util.h"
#include <tuple>

#include "core/hook_registry.h"
#include "core/hook_traits.h"
#include "core/operation/operation.h"
#include "core/wound_iteration.h"
#include "infra/sim_param.h"

namespace bdm {
namespace skibidy {

// FusedWoundSourceOp: pure iteration skeleton. All domain-specific biology is
// delegated to self-contained module hooks registered in hook_registry.h.
// Adding a new source module requires only creating a source_hook.h and
// adding the type to the SourceHooks tuple.
//
// Architecture (inspired by Unreal Engine Subsystems):
//   GridRegistry  : single-point grid cache (filled once per step)
//   VoxelSnapshot : per-voxel read-once context (L1-friendly)
//   SignalBoard   : cross-hook coupling channel (Delegate pattern)
//   hook_traits   : compile-time dispatch (if constexpr eliminates dead code)
struct FusedWoundSourceOp : public StandaloneOperationImpl {
  BDM_OP_HEADER(FusedWoundSourceOp);

  void operator()() override {
    auto* sim = Simulation::GetActive();
    auto* sp = sim->GetParam()->Get<SimParam>();
    PerfTimer timer(sp->debug_perf);

    ctx_.Prepare(sim, sp);

    SourceHooks src_hooks;
    InitHooks(src_hooks, ctx_.reg, ctx_.sig);
    StudySourceHooks study_src_hooks;
    InitHooks(study_src_hooks, ctx_.reg, ctx_.sig);

    DecayHooks decay_hooks;
    InitHooks(decay_hooks, ctx_.reg, ctx_.sig);
    StudyDecayHooks study_decay_hooks;
    InitHooks(study_decay_hooks, ctx_.reg, ctx_.sig);

    ctx_.IterateDermal([&](VoxelSnapshot& snap, SignalBoard& sig) {
      std::apply([&](auto&... h) {
        (DispatchDermal(h, snap, sig), ...);
      }, src_hooks);
      std::apply([&](auto&... h) {
        (DispatchDermal(h, snap, sig), ...);
      }, study_src_hooks);
    });

    ctx_.IterateEpiWound([&](VoxelSnapshot& snap, SignalBoard& sig) {
      std::apply([&](auto&... h) {
        (DispatchEpiWound(h, snap, sig), ...);
      }, src_hooks);
      std::apply([&](auto&... h) {
        (DispatchEpiWound(h, snap, sig), ...);
      }, study_src_hooks);
    });

    real_t dt = sim->GetParam()->simulation_time_step;
    std::apply([&](auto&... h) {
      (DispatchDecay(h, ctx_.mask, dt), ...);
    }, decay_hooks);
    std::apply([&](auto&... h) {
      (DispatchDecay(h, ctx_.mask, dt), ...);
    }, study_decay_hooks);

    timer.Print("fused_source");
  }

 private:
  WoundIterationContext ctx_;
};

}  // namespace skibidy
}  // namespace bdm

#endif  // FUSED_SOURCE_H_
