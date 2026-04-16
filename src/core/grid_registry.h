#ifndef GRID_REGISTRY_H_
#define GRID_REGISTRY_H_

#include <array>
#include <cmath>
#include <vector>
#include "core/field_names.h"
#include "core/material.h"
#include "core/pde.h"
#include "infra/sim_param.h"
#include "infra/util.h"

namespace bdm {
namespace skibidy {

// ---------------------------------------------------------------------------
// GridRegistry -- single-point grid pointer cache (Unreal Subsystem pattern).
// Filled once per step. All hooks read from this instead of calling
// rm->GetDiffusionGrid() independently. Eliminates 100+ redundant lookups.
// ---------------------------------------------------------------------------
struct GridRegistry {
  static constexpr int kMaxFields = 48;

  std::array<DiffusionGrid*, kMaxFields> grids_{};
  const SimParam* sp_ = nullptr;
  uint64_t step_ = 0;
  real_t dt_ = 0.1;
  MaterialRegistry mat_reg_;
  bool mat_initialized_ = false;

  bool post_wound_ = false;
  uint64_t wound_age_ = 0;

  bool coarse_ = false;
  real_t coarse_w_ = 1.0;
  std::vector<size_t> coarse_map_;

  void Fill(Simulation* sim) {
    auto* rm = sim->GetResourceManager();
    sp_ = sim->GetParam()->Get<SimParam>();
    step_ = GetGlobalStep(sim);
    dt_ = sim->GetParam()->simulation_time_step;
    grids_.fill(nullptr);
    rm->ForEachContinuum([&](Continuum* cm) {
      auto* dg = dynamic_cast<DiffusionGrid*>(cm);
      if (dg) {
        int id = static_cast<int>(dg->GetContinuumId());
        if (id >= 0 && id < kMaxFields) grids_[id] = dg;
      }
    });

    // Wound state (derived once per step, consumed by many hooks)
    uint64_t wound_step = static_cast<uint64_t>(sp_->wound.trigger_step);
    post_wound_ = (step_ > wound_step);
    wound_age_ = post_wound_ ? step_ - wound_step : 0;

    // Coarse structural grid map (lazy, computed once)
    coarse_ = sp_->grid_resolution_structural > 0 &&
              sp_->grid_resolution_structural != sp_->grid_resolution;
    if (coarse_) {
      real_t rf = static_cast<real_t>(sp_->grid_resolution);
      real_t rc = static_cast<real_t>(sp_->grid_resolution_structural);
      coarse_w_ = (rc * rc * rc) / (rf * rf * rf);
      if (coarse_map_.empty()) {
        DiffusionGrid* ref = Get(fields::kCollagenId);
        if (!ref) ref = Get(fields::kFibronectinId);
        if (!ref) ref = Get(fields::kScarId);
        if (ref) {
          auto* strat = Get(fields::kStratumId);
          GridContext ctx(strat, sp_);
          coarse_map_.resize(ctx.n);
          for (size_t i = 0; i < ctx.n; i++) {
            Real3 pos = {ctx.X(i), ctx.Y(i), ctx.Z(i)};
            coarse_map_[i] = ref->GetBoxIndex(pos);
          }
        }
      }
    } else {
      coarse_w_ = 1.0;
    }

    // Lazy-init material registry (once, not per step)
    if (!mat_initialized_) {
      mat_reg_.RegisterSkinDefaults();
      if (sp_->photon.enabled) {
        mat_reg_.RegisterBrainDefaults();
      }
      if (sp_->ra.enabled) {
        mat_reg_.RegisterJointDefaults();
      }
      mat_initialized_ = true;
    }
  }

  // Field access (returns nullptr if field not registered).
  DiffusionGrid* Get(int field_id) const { return grids_[field_id]; }

  const SimParam* Params() const { return sp_; }
  uint64_t Step() const { return step_; }
  real_t Dt() const { return dt_; }
  const MaterialRegistry& Materials() const { return mat_reg_; }

  bool PostWound() const { return post_wound_; }
  uint64_t WoundAge() const { return wound_age_; }

  // Pre-tapered wound inflammation rate (0 if pre-wound or decayed below 1e-10).
  real_t WoundInflRate() const {
    if (!post_wound_ || sp_->wound.inflammation_source_rate <= 0) return 0;
    real_t rate = sp_->wound.inflammation_source_rate;
    real_t taper = sp_->wound.inflammation_source_taper;
    if (taper > 0) {
      rate *= std::exp(-taper * static_cast<real_t>(wound_age_));
      if (rate < 1e-10) return 0;
    }
    return rate;
  }

  // Coarse-to-fine structural grid index mapping. Returns fine_idx unchanged
  // if no structural resolution is configured.
  size_t CoarseIndex(size_t fine_idx) const {
    if (!coarse_ || coarse_map_.empty()) return fine_idx;
    return coarse_map_[fine_idx];
  }
  bool HasCoarseGrid() const { return coarse_; }
  real_t CoarseWeight() const { return coarse_w_; }

  // Resolve split vs unified inflammation grid.
  DiffusionGrid* InflammationGrid() const {
    if (sp_->inflammation.split_inflammation_enabled)
      return grids_[fields::kProInflammatoryId];
    return grids_[fields::kInflammationId];
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // GRID_REGISTRY_H_
