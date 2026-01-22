#ifndef FIBROBLAST_RECRUITMENT_H_
#define FIBROBLAST_RECRUITMENT_H_

#include <cmath>

#include "fibroblast/fibroblast.h"
#include "fibroblast/fibroblast_behavior.h"
#include "infra/sim_param.h"
#include "core/operation/operation.h"

namespace bdm {
namespace skibidy {

// ---------------------------------------------------------------------------
// FibroblastRecruitment -- standalone operation that spawns resident dermal
// fibroblasts in a ring around the wound across multiple recruitment waves.
// Staggered arrival produces gradual myofibroblast differentiation instead
// of a single spike, matching literature kinetics (Desmouliere 1995).
// ---------------------------------------------------------------------------
struct FibroblastRecruitment : public StandaloneOperationImpl {
  BDM_OP_HEADER(FibroblastRecruitment);

  void operator()() override {
    auto* sim = Simulation::GetActive();
    auto* scheduler = sim->GetScheduler();
    auto* sp = sim->GetParam()->Get<SimParam>();

    if (!sp->wound_enabled || !sp->fibroblast_enabled) return;

    uint64_t step = scheduler->GetSimulatedSteps();
    uint64_t wound_step = static_cast<uint64_t>(sp->wound_trigger_step);
    if (step <= wound_step) return;

    uint64_t wound_age = step - wound_step;
    int waves = std::max(1, sp->fibroblast_spawn_waves);

    if (waves_spawned_ >= waves) return;  // all waves done

    // Compute wave trigger time
    uint64_t base_delay = static_cast<uint64_t>(sp->fibroblast_spawn_delay);
    uint64_t window = static_cast<uint64_t>(sp->fibroblast_spawn_window);
    uint64_t wave_interval = (waves > 1) ? window / (waves - 1) : 0;
    uint64_t wave_trigger = base_delay + waves_spawned_ * wave_interval;

    if (wound_age >= wave_trigger) {
      // Compute total cell count (same formula as before)
      if (total_cells_ < 0) {
        real_t outer_r = sp->wound_radius + sp->dermal_fibroblast_margin;
        total_cells_ = static_cast<int>(
            std::ceil(2.0 * M_PI * outer_r / sp->fibroblast_diameter *
                      sp->fibroblast_density_factor));
        if (total_cells_ < 6) total_cells_ = 6;
      }

      int cells_per_wave = total_cells_ / waves;
      int this_wave = cells_per_wave;
      // Last wave gets remainder
      if (waves_spawned_ == waves - 1) {
        this_wave = total_cells_ - cells_per_wave * (waves - 1);
      }

      SpawnFibroblasts(sim, sp, this_wave, waves_spawned_);
      waves_spawned_++;
    }
  }

 private:
  int waves_spawned_ = 0;
  int total_cells_ = -1;  // computed once on first wave

  void SpawnFibroblasts(Simulation* sim, const SimParam* sp,
                        int count, int wave_idx) {
    auto* ctxt = sim->GetExecutionContext();
    auto* random = sim->GetRandom();

    real_t cx = sp->wound_center_x;
    real_t cy = sp->wound_center_y;
    real_t r = sp->wound_radius;
    real_t d = sp->fibroblast_diameter;
    real_t z = sp->dermal_fibroblast_depth;  // in dermis (z < 0)
    real_t outer_r = r + sp->dermal_fibroblast_margin;

    // Offset angular start per wave so cells don't overlap
    real_t d_theta = 2.0 * M_PI / count;
    real_t theta_offset = wave_idx * (d_theta * 0.5);

    for (int i = 0; i < count; i++) {
      real_t theta = theta_offset + i * d_theta;
      real_t frac = random->Uniform(0, 1);
      real_t ring_r = r + frac * sp->dermal_fibroblast_margin;
      real_t x = cx + ring_r * std::cos(theta);
      real_t y = cy + ring_r * std::sin(theta);

      x = std::max(sp->tissue_min + d, std::min(x, sp->tissue_max - d));
      y = std::max(sp->tissue_min + d, std::min(y, sp->tissue_max - d));

      auto* cell = new Fibroblast({x, y, z});
      cell->SetDiameter(d);
      cell->SetFibroblastState(kFibroQuiescent);
      cell->AddBehavior(new FibroblastBehavior());
      ctxt->AddAgent(cell);
    }

    if (sp->debug_fibroblast) {
      std::cout << "[fibroblast] wave " << (wave_idx + 1)
                << ": spawned " << count
                << " (z=" << z << ", r=" << r << "-" << (r + sp->dermal_fibroblast_margin) << ")"
                << std::endl;
    }
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // FIBROBLAST_RECRUITMENT_H_
