#ifndef WOUND_ITERATION_H_
#define WOUND_ITERATION_H_

#include <cstddef>
#include <vector>

#include "biodynamo.h"
#include "core/field_names.h"
#include "core/grid_registry.h"
#include "core/pde.h"
#include "core/signal_board.h"
#include "core/voxel_snapshot.h"
#include "infra/sim_param.h"

namespace bdm {
namespace skibidy {

// Shared per-step iteration scaffold for the fused source and post ops.
// Owns the per-invocation registry/signal board plus the cached wound mask
// and snapshot filler. Coarse map now lives in GridRegistry.
struct WoundIterationContext {
  GridRegistry reg;
  SignalBoard sig;

  GridContext::WoundMaskData mask;
  GridContext ctx;
  SnapshotFiller filler;

  void Prepare(Simulation* sim, const SimParam* sp) {
    reg.Fill(sim);

    ctx = GridContext(reg.Get(fields::kStratumId), sp);
    if (mask.mask.empty()) {
      real_t z_max = sp->volume_z_cornified + ctx.box_len;
      mask = GridContext::ComputeWoundMask(ctx, z_max, sp);
    }

    // SnapshotFiller still takes coarse params for its own CoarseSI(); pass
    // from registry so there's one source of truth.
    static const std::vector<size_t> empty_map;
    const auto& cm = reg.HasCoarseGrid() ? reg.coarse_map_ : empty_map;
    filler.Init(reg, ctx, mask, cm, reg.HasCoarseGrid(), reg.CoarseWeight());
  }

  template <class F>
  void IterateDermal(F&& fn) {
    VoxelSnapshot snap;
    for (size_t idx : mask.dermal_all) {
      filler.FillDermal(snap, idx);
      fn(snap, sig);
    }
  }

  template <class F>
  void IterateEpiWound(F&& fn) {
    VoxelSnapshot snap;
    for (size_t idx : mask.epi_wound) {
      filler.FillEpi(snap, idx);
      fn(snap, sig);
    }
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // WOUND_ITERATION_H_
