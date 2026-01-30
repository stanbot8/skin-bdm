#ifndef VEGF_H_
#define VEGF_H_

#include <cmath>

#include "core/field_names.h"
#include "core/pde.h"
#include "infra/sim_param.h"
#include "core/operation/operation.h"

namespace bdm {
namespace skibidy {

// ---------------------------------------------------------------------------
// VEGFSourceOp -- hypoxia-driven VEGF production in wound tissue.
// Reads O2 field; where O2 < threshold, produces VEGF proportional to
// hypoxia severity. Diabetic mode reduces production (HIF-1alpha impairment).
// VEGF consumption by angiogenesis is handled in VascularPDE::ApplySource().
// ---------------------------------------------------------------------------
struct VEGFSourceOp : public StandaloneOperationImpl {
  BDM_OP_HEADER(VEGFSourceOp);

  void operator()() override {
    auto* sim = Simulation::GetActive();
    auto* sp = sim->GetParam()->Get<SimParam>();

    if (!sp->angiogenesis_enabled || !sp->wound_enabled) return;

    auto* scheduler = sim->GetScheduler();
    uint64_t step = scheduler->GetSimulatedSteps();
    uint64_t wound_step = static_cast<uint64_t>(sp->wound_trigger_step);
    if (step <= wound_step) return;

    auto* rm = sim->GetResourceManager();
    auto* o2_grid = rm->GetDiffusionGrid(fields::kOxygen);
    auto* vegf_grid = rm->GetDiffusionGrid(fields::kVEGF);
    GridContext ctx(vegf_grid, sp);

    real_t threshold = sp->vegf_hypoxia_threshold;
    real_t prod_rate = sp->vegf_production_rate;
    if (sp->diabetic_mode) {
      prod_rate *= sp->diabetic_vegf_factor;
    }
    real_t z_max = sp->volume_z_cornified + ctx.box_len;

    for (size_t idx = 0; idx < ctx.n; idx++) {
      real_t z = ctx.Z(idx);
      // Both dermal and epidermal wound voxels produce VEGF under hypoxia
      if (z > z_max) continue;
      real_t x = ctx.X(idx), y = ctx.Y(idx);
      if (!ctx.InWound(x, y)) continue;

      real_t local_o2 = o2_grid->GetConcentration(idx);
      if (local_o2 < threshold) {
        // Linear ramp: max production at O2=0, zero at O2=threshold
        real_t vegf_prod = prod_rate * (threshold - local_o2) / threshold;
        if (vegf_prod > 1e-10) {
          vegf_grid->ChangeConcentrationBy(idx, vegf_prod);
        }
      }
    }
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // VEGF_H_
