#ifndef HOT_RELOAD_H_
#define HOT_RELOAD_H_

#include <sys/stat.h>

#include <fstream>
#include <iostream>
#include <string>

#include "biodynamo.h"
#include "core/field_names.h"
#include "infra/sim_param.h"
#include "core/operation/operation.h"

namespace bdm {
namespace skibidy {

// ---------------------------------------------------------------------------
// HotReloadOp -- checks bdm.toml for modifications and re-applies parameter
// values at runtime. Enables interactive tuning without restarting the sim.
//
// Checked every metrics_interval steps (default: every 100 steps = 10h sim).
// Only safe-to-change parameters are updated (rates, thresholds, flags).
// Grid topology and domain bounds are NOT changed (would corrupt state).
// ---------------------------------------------------------------------------
struct HotReloadOp : public StandaloneOperationImpl {
  BDM_OP_HEADER(HotReloadOp);

  void operator()() override {
    auto* sim = Simulation::GetActive();
    auto* sp = sim->GetParam()->Get<SimParam>();
    if (!sp->hot_reload) return;

    uint64_t step = sim->GetScheduler()->GetSimulatedSteps();
    int interval = sp->metrics_interval > 0 ? sp->metrics_interval : 100;
    if (step % static_cast<uint64_t>(interval) != 0) return;

    // Check file modification time
    struct stat st;
    if (stat("bdm.toml", &st) != 0) return;
    if (last_mtime_ != 0 && st.st_mtime <= last_mtime_) {
      return;
    }
    if (last_mtime_ == 0) {
      // First call: record mtime without reloading
      last_mtime_ = st.st_mtime;
      return;
    }
    last_mtime_ = st.st_mtime;

    try {
      auto config = toml::parse_file("bdm.toml");
      auto* sp_mut = const_cast<SimParam*>(sp);
      sp_mut->LoadConfig(config);
      sp_mut->ValidateConfig();
      UpdateGridDecayRates(sim, sp_mut);
      std::cout << "[hot-reload] parameters updated at step " << step
                << " (day " << step * sim->GetParam()->simulation_time_step / 24.0
                << ")" << std::endl;
    } catch (const std::exception& e) {
      std::cout << "[hot-reload] parse error: " << e.what() << std::endl;
    }
  }

 private:
  time_t last_mtime_ = 0;

  void UpdateGridDecayRates(Simulation* sim, const SimParam* sp) {
    auto* rm = sim->GetResourceManager();
    if (sp->mmp_enabled) {
      auto* grid = rm->GetDiffusionGrid(fields::kMMP);
      if (grid) {
        real_t decay = sp->diabetic_mode
                           ? sp->mmp_decay * sp->diabetic_timp_factor
                           : sp->mmp_decay;
        grid->SetDecayConstant(decay);
      }
    }
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // HOT_RELOAD_H_
