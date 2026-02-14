#ifndef VOXEL_ENV_H_
#define VOXEL_ENV_H_

// Fused per-voxel environment snapshot for agent behaviors.
//
// Instead of each behavior independently calling GetDiffusionGrid(id)->
// GetValue(pos) (redundant pos-to-index conversion, scattered cache misses
// across 28 separate grid arrays), behaviors read from this compact
// array-of-structs buffer. One cache-line fetch brings all fields together.
//
// Filled once per step by VoxelEnvFillOp (kPreSchedule, before behaviors).
// Read-only during behavior dispatch. Behaviors that WRITE to grids still
// use the BDM DiffusionGrid API directly.

#include <vector>

#include "biodynamo.h"
#include "core/field_names.h"
#include "infra/sim_param.h"
#include "infra/util.h"

namespace bdm {
namespace skibidy {

// Packed per-voxel snapshot of all fields that agent behaviors read.
// 13 doubles x 8 bytes = 104 bytes/voxel. At 1000 voxels = 101.6 KB (fits L2).
struct VoxelEnv {
  real_t o2;
  real_t kgf;
  real_t water;
  real_t calcium;
  real_t immune_pressure;
  real_t stratum;
  real_t glucose;
  real_t ecm_quality;
  real_t ph;
  real_t temperature;
  real_t inflammation;
  real_t tgf_beta;
  real_t biofilm;
};

// Global fused environment cache.
struct VoxelEnvCache {
  std::vector<VoxelEnv> data{1};       // always >= 1 element (zero fallback)
  DiffusionGrid* ref_grid = nullptr;   // any grid, for GetBoxIndex
  Simulation* filled_sim = nullptr;    // tracks which sim owns ref_grid
  size_t num_boxes = 1;

  // Position to flat voxel index (one call replaces 3-5 GetValue calls).
  // Lazy fill when called from a new sim (tests that skip VoxelEnvFillOp).
  size_t Index(const Real3& pos) {
    if (filled_sim != Simulation::GetActive()) {
      Fill(Simulation::GetActive());
    }
    if (!ref_grid) return 0;
    return ref_grid->GetBoxIndex(pos);
  }

  const VoxelEnv& operator[](size_t idx) const { return data[idx]; }

  void Fill(Simulation* sim) {
    auto* rm = sim->GetResourceManager();
    auto* sp = sim->GetParam()->Get<SimParam>();

    // Reference grid: use O2 (always registered in real sims).
    // Fall back to any available grid for tests that register only a subset.
    ref_grid = rm->GetDiffusionGrid(fields::kOxygenId);
    if (!ref_grid) {
      for (int id = 0; id <= fields::kBasalDensityId && !ref_grid; ++id) {
        ref_grid = rm->GetDiffusionGrid(id);
      }
    }
    if (!ref_grid) {
      data.assign(1, {});
      num_boxes = 1;
      filled_sim = sim;
      return;
    }
    filled_sim = sim;

    num_boxes = ref_grid->GetNumBoxes();
    data.resize(num_boxes);

    // Grid pointers (null if module disabled)
    auto* kgf_grid   = rm->GetDiffusionGrid(fields::kKGFId);
    auto* water_grid  = rm->GetDiffusionGrid(fields::kWaterId);
    auto* ca_grid     = rm->GetDiffusionGrid(fields::kCalciumId);
    auto* ip_grid     = rm->GetDiffusionGrid(fields::kImmunePressureId);
    auto* strat_grid  = rm->GetDiffusionGrid(fields::kStratumId);
    auto* gluc_grid   = (sp->glucose_enabled && sp->diabetic_mode)
                          ? rm->GetDiffusionGrid(fields::kGlucoseId)
                          : nullptr;
    auto* ph_grid     = rm->GetDiffusionGrid(fields::kPHId);
    auto* temp_grid   = sp->temperature_enabled
                          ? rm->GetDiffusionGrid(fields::kTemperatureId)
                          : nullptr;
    auto* infl_grid   = sp->split_inflammation_enabled
                          ? rm->GetDiffusionGrid(fields::kProInflammatoryId)
                          : rm->GetDiffusionGrid(fields::kInflammationId);
    auto* tgf_grid    = sp->fibroblast_enabled
                          ? rm->GetDiffusionGrid(fields::kTGFBetaId)
                          : nullptr;
    auto* bio_grid    = sp->biofilm_enabled
                          ? rm->GetDiffusionGrid(fields::kBiofilmId)
                          : nullptr;

    // Multi-resolution: structural grids (pH, biofilm) may use coarser
    // resolution. Detect and handle by mapping fine index -> world pos ->
    // coarse index for those grids.
    bool coarse = sp->grid_resolution_structural > 0 &&
                  sp->grid_resolution_structural != sp->grid_resolution;
    // Helper: fine-grid voxel center coordinates for cross-resolution lookup
    size_t res = ref_grid->GetResolution();
    real_t lo = static_cast<real_t>(ref_grid->GetDimensions()[0]);
    real_t bl = ref_grid->GetBoxLength();
    auto fine_pos = [res, lo, bl](size_t i) -> Real3 {
      real_t x = lo + (i % res) * bl + bl / 2.0;
      real_t y = lo + ((i / res) % res) * bl + bl / 2.0;
      real_t z = lo + (i / (res * res)) * bl + bl / 2.0;
      return {x, y, z};
    };
    // Resolve fine index to coarse-grid concentration
    auto coarse_get = [&](DiffusionGrid* g, size_t i) -> real_t {
      if (!coarse || g->GetNumBoxes() == num_boxes) {
        return g->GetConcentration(i);
      }
      return g->GetValue(fine_pos(i));
    };

    // Sequential scan: each source grid is 8 KB (1000 doubles). The tight
    // loop ensures both source and destination stay in cache.
    for (size_t i = 0; i < num_boxes; ++i) {
      auto& v = data[i];
      v.o2              = ref_grid->GetConcentration(i);
      v.kgf             = kgf_grid   ? kgf_grid->GetConcentration(i)   : 0;
      v.water           = water_grid ? water_grid->GetConcentration(i) : 0;
      v.calcium         = ca_grid    ? ca_grid->GetConcentration(i)    : 0;
      v.immune_pressure = ip_grid    ? ip_grid->GetConcentration(i)    : 0;
      v.stratum         = strat_grid ? strat_grid->GetConcentration(i) : 0;
      v.glucose         = gluc_grid  ? gluc_grid->GetConcentration(i)  : 0;
      v.ph              = ph_grid    ? coarse_get(ph_grid, i)          : 0;
      v.temperature     = temp_grid  ? temp_grid->GetConcentration(i)  : 1.0;
      v.inflammation    = infl_grid  ? infl_grid->GetConcentration(i)  : 0;
      v.tgf_beta        = tgf_grid   ? tgf_grid->GetConcentration(i)  : 0;
      v.biofilm         = bio_grid   ? coarse_get(bio_grid, i)         : 0;
      v.ecm_quality     = g_ecm_quality
                            ? g_ecm_quality->GetConcentration(i) : 0;
    }
  }
};

// Single global instance.
inline VoxelEnvCache g_voxel_env;

// Pre-schedule operation: fills g_voxel_env once per step before behaviors.
struct VoxelEnvFillOp : public StandaloneOperationImpl {
  BDM_OP_HEADER(VoxelEnvFillOp);

  void operator()() override {
    g_voxel_env.Fill(Simulation::GetActive());
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // VOXEL_ENV_H_
