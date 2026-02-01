#ifndef BASELINE_INFLAMMATION_H_
#define BASELINE_INFLAMMATION_H_

#include "core/field_names.h"
#include "core/pde.h"
#include "infra/sim_param.h"
#include "core/operation/operation.h"

namespace bdm {
namespace skibidy {

// ---------------------------------------------------------------------------
// BaselineInflammationOp -- diabetic chronic low-grade inflammation.
// Adds a small amount of pro-inflammatory signal to wound epidermal voxels
// each step, modeling the AGE/RAGE-driven sterile inflammation present in
// diabetic tissue even before wounding.
// ---------------------------------------------------------------------------
struct BaselineInflammationOp : public StandaloneOperationImpl {
  BDM_OP_HEADER(BaselineInflammationOp);

  void operator()() override {
    auto* sim = Simulation::GetActive();
    auto* sp = sim->GetParam()->Get<SimParam>();

    if (!sp->diabetic_mode || !sp->wound_enabled) return;
    if (sp->diabetic_baseline_inflammation <= 0) return;

    auto* rm = sim->GetResourceManager();

    DiffusionGrid* infl_grid = nullptr;
    if (sp->split_inflammation_enabled) {
      infl_grid = rm->GetDiffusionGrid(fields::kProInflammatory);
    } else {
      infl_grid = rm->GetDiffusionGrid(fields::kInflammation);
    }

    GridContext ctx(infl_grid, sp);
    real_t z_max = sp->volume_z_cornified + ctx.box_len;
    real_t rate = sp->diabetic_baseline_inflammation;

    auto* stratum_grid = rm->GetDiffusionGrid(fields::kStratum);

    for (size_t idx = 0; idx < ctx.n; idx++) {
      real_t z = ctx.Z(idx);
      if (z < 0 || z > z_max) continue;
      if (!ctx.InWound(ctx.X(idx), ctx.Y(idx))) continue;

      // Gated by stratum: re-epithelialized tissue seals the AGE/RAGE source.
      real_t sv = stratum_grid->GetConcentration(idx);
      real_t gate = std::max(static_cast<real_t>(0),
                             static_cast<real_t>(1) - sv);
      if (gate > 1e-10) {
        infl_grid->ChangeConcentrationBy(idx, gate * rate);
      }
    }
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // BASELINE_INFLAMMATION_H_
