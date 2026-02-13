#ifndef DERIVED_FIELDS_OP_H_
#define DERIVED_FIELDS_OP_H_

#include "core/derived_field.h"
#include "core/field_names.h"
#include "core/pde.h"
#include "infra/sim_param.h"
#include "core/operation/operation.h"

namespace bdm {
namespace skibidy {

// Recomputes all derived composite fields once per step.
// Runs as PostSchedule after FusedWoundPostOp but before MetricsExporter.
// Single fused loop: reads all source grids, writes all composites in one pass.
struct DerivedFieldsOp : public StandaloneOperationImpl {
  BDM_OP_HEADER(DerivedFieldsOp);

  DerivedFieldsOp(DerivedField* ecm_quality,
                  DerivedField* tissue_viability,
                  DerivedField* wound_microenv)
      : ecm_quality_(ecm_quality),
        tissue_viability_(tissue_viability),
        wound_microenv_(wound_microenv) {}

  void operator()() override {
    auto* sim = Simulation::GetActive();
    auto* sp = sim->GetParam()->Get<SimParam>();
    PerfTimer timer(sp->debug_perf);
    auto* rm = sim->GetResourceManager();

    // Source grids (null when module disabled)
    DiffusionGrid* col_grid = sp->fibroblast_enabled
        ? rm->GetDiffusionGrid(fields::kCollagen) : nullptr;
    DiffusionGrid* fn_grid = sp->fibronectin_enabled
        ? rm->GetDiffusionGrid(fields::kFibronectin) : nullptr;
    DiffusionGrid* el_grid = sp->elastin_enabled
        ? rm->GetDiffusionGrid(fields::kElastin) : nullptr;
    DiffusionGrid* fb_grid = sp->hemostasis_enabled
        ? rm->GetDiffusionGrid(fields::kFibrin) : nullptr;
    DiffusionGrid* o2_grid = rm->GetDiffusionGrid(fields::kOxygen);
    DiffusionGrid* water_grid = rm->GetDiffusionGrid(fields::kWater);
    DiffusionGrid* vasc_grid = rm->GetDiffusionGrid(fields::kVascular);

    // Weights from config
    real_t w_col = sp->ecm_weight_collagen;
    real_t w_fn = sp->ecm_weight_fibronectin;
    real_t w_el = sp->ecm_weight_elastin;
    real_t w_fb = sp->ecm_weight_fibrin;

    size_t n = ecm_quality_->GetNumBoxes();

    // Handle multi-resolution: source grids may differ in resolution.
    // If structural grids are coarser, use GetValue via world coordinates
    // instead of direct index access.
    bool coarse = sp->grid_resolution_structural > 0 &&
                  sp->grid_resolution_structural != sp->grid_resolution;

    // Build a fine GridContext for coordinate conversion
    GridContext ctx(o2_grid, sp);

    // Single fused loop over all voxels
    for (size_t idx = 0; idx < n; ++idx) {
      // ECM scaffold quality
      real_t ecm = 0;
      if (coarse) {
        // Structural grids at coarse resolution: use coordinate lookup
        real_t x = ctx.X(idx), y = ctx.Y(idx), z = ctx.Z(idx);
        Real3 pos = {x, y, z};
        if (col_grid) ecm += w_col * col_grid->GetValue(pos);
        if (fn_grid) ecm += w_fn * fn_grid->GetValue(pos);
        if (el_grid) ecm += w_el * el_grid->GetValue(pos);
        if (fb_grid) ecm += w_fb * fb_grid->GetValue(pos);
      } else {
        // Same resolution: direct index access
        if (col_grid) ecm += w_col * col_grid->GetConcentration(idx);
        if (fn_grid) ecm += w_fn * fn_grid->GetConcentration(idx);
        if (el_grid) ecm += w_el * el_grid->GetConcentration(idx);
        if (fb_grid) ecm += w_fb * fb_grid->GetConcentration(idx);
      }
      ecm_quality_->SetConcentration(idx, ecm);

      // Tissue viability: O2 * water * perfusion (all fine-grid)
      real_t o2 = o2_grid->GetConcentration(idx);
      real_t water = water_grid->GetConcentration(idx);
      real_t vasc = vasc_grid->GetConcentration(idx);
      real_t viability = o2 * water * vasc;
      tissue_viability_->SetConcentration(idx, viability);

      // Wound microenvironment: viability * ecm (extensible)
      wound_microenv_->SetConcentration(idx, viability * ecm);
    }
    timer.Print("derived_fields");
  }

 private:
  DerivedField* ecm_quality_;
  DerivedField* tissue_viability_;
  DerivedField* wound_microenv_;
};

}  // namespace skibidy
}  // namespace bdm

#endif  // DERIVED_FIELDS_OP_H_
