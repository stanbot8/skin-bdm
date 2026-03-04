#ifndef FUSED_SOURCE_H_
#define FUSED_SOURCE_H_

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
// FusedWoundSourceOp -- pure iteration skeleton.
// All domain-specific biology is delegated to self-contained module hooks
// registered in hook_registry.h. Adding a new source module requires only
// creating a source_hook.h and adding the type to the SourceHooks tuple.
//
// Architecture (inspired by Unreal Engine Subsystems):
//   GridRegistry  -- single-point grid cache (filled once per step)
//   VoxelSnapshot -- per-voxel read-once context (L1-friendly)
//   SignalBoard   -- cross-hook coupling channel (Delegate pattern)
//   hook_traits   -- compile-time dispatch (if constexpr eliminates dead code)
// ---------------------------------------------------------------------------
struct FusedWoundSourceOp : public StandaloneOperationImpl {
  BDM_OP_HEADER(FusedWoundSourceOp);

  void operator()() override {
    auto* sim = Simulation::GetActive();
    auto* sp = sim->GetParam()->Get<SimParam>();
    PerfTimer timer(sp->debug_perf);

    // --- Registry + SignalBoard ---
    GridRegistry reg;
    reg.Fill(sim);
    SignalBoard sig;

    // --- Hook construction + Init (engine + study) ---
    SourceHooks src_hooks;
    InitHooks(src_hooks, reg, sig);
    StudySourceHooks study_src_hooks;
    InitHooks(study_src_hooks, reg, sig);

    DecayHooks decay_hooks;
    InitHooks(decay_hooks, reg, sig);
    StudyDecayHooks study_decay_hooks;
    InitHooks(study_decay_hooks, reg, sig);

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

    // --- Loop 1: ALL dermal voxels ---
    VoxelSnapshot snap;
    for (size_t idx : mask_.dermal_all) {
      filler.FillDermal(snap, idx);
      std::apply([&](auto&... h) {
        (DispatchDermal(h, snap, sig), ...);
      }, src_hooks);
      std::apply([&](auto&... h) {
        (DispatchDermal(h, snap, sig), ...);
      }, study_src_hooks);
    }

    // --- Loop 2: epidermal wound voxels ---
    for (size_t idx : mask_.epi_wound) {
      filler.FillEpi(snap, idx);
      std::apply([&](auto&... h) {
        (DispatchEpiWound(h, snap, sig), ...);
      }, src_hooks);
      std::apply([&](auto&... h) {
        (DispatchEpiWound(h, snap, sig), ...);
      }, study_src_hooks);
    }

    // --- ECM decay hooks (self-contained iteration) ---
    real_t dt = sim->GetParam()->simulation_time_step;
    std::apply([&](auto&... h) {
      (DispatchDecay(h, mask_, dt), ...);
    }, decay_hooks);
    std::apply([&](auto&... h) {
      (DispatchDecay(h, mask_, dt), ...);
    }, study_decay_hooks);

    timer.Print("fused_source");
  }

 private:
  GridContext::WoundMaskData mask_;
  std::vector<size_t> coarse_map_;
};

}  // namespace skibidy
}  // namespace bdm

#endif  // FUSED_SOURCE_H_
