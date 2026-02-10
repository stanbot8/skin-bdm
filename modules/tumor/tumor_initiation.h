#ifndef TUMOR_INITIATION_H_
#define TUMOR_INITIATION_H_

#include <cmath>
#include <iostream>

#include "tumor/tumor_cell.h"
#include "tumor/tumor_behavior.h"
#include "infra/sim_param.h"
#include "core/operation/operation.h"

namespace bdm {
namespace skibidy {

// ---------------------------------------------------------------------------
// TumorInitiation -- standalone operation that seeds a cluster of tumor cells
// at a configurable location and step.  Independent of wound_enabled: composes
// with other events purely through bdm.toml config.
// ---------------------------------------------------------------------------
struct TumorInitiation : public StandaloneOperationImpl {
  BDM_OP_HEADER(TumorInitiation);

  void operator()() override {
    auto* sim = Simulation::GetActive();
    auto* scheduler = sim->GetScheduler();
    auto* sp = sim->GetParam()->Get<SimParam>();

    if (!sp->tumor_enabled || fired_) return;

    uint64_t step = scheduler->GetSimulatedSteps();
    if (step < static_cast<uint64_t>(sp->tumor_seed_time)) return;

    SpawnTumor(sim, sp);
    fired_ = true;
  }

 private:
  bool fired_ = false;

  void SpawnTumor(Simulation* sim, const SimParam* sp) {
    auto* ctxt = sim->GetExecutionContext();
    auto* random = sim->GetRandom();

    real_t cx = sp->tumor_seed_x;
    real_t cy = sp->tumor_seed_y;
    real_t cz = sp->tumor_seed_z;
    real_t d = sp->tumor_diameter;
    int n = sp->tumor_seed_count;

    // Sphere radius: scales with cube root of cell count (loose packing ~50%)
    real_t R = (d / 2.0) * std::cbrt(static_cast<real_t>(n) / 0.5);

    for (int i = 0; i < n; i++) {
      // Uniform random point in sphere of radius R
      real_t u = std::cbrt(random->Uniform(0, 1));
      real_t theta = random->Uniform(0, 2.0 * M_PI);
      real_t cos_phi = random->Uniform(-1.0, 1.0);
      real_t sin_phi = std::sqrt(1.0 - cos_phi * cos_phi);
      real_t r = R * u;

      real_t x = cx + r * sin_phi * std::cos(theta);
      real_t y = cy + r * sin_phi * std::sin(theta);
      real_t z = cz + r * cos_phi;

      // Clamp to domain bounds (tumor can invade dermis, so no z floor)
      auto* param = sim->GetParam();
      real_t margin = d;
      x = std::max(param->min_bound + margin, std::min(x, param->max_bound - margin));
      y = std::max(param->min_bound + margin, std::min(y, param->max_bound - margin));
      z = std::max(param->min_bound + margin, std::min(z, param->max_bound - margin));

      auto* cell = new TumorCell({x, y, z});
      cell->SetDiameter(d);
      cell->AddBehavior(new TumorBehavior());
      ctxt->AddAgent(cell);
    }

    std::cout << "[TumorInitiation] Seeded " << n
              << " tumor cells at (" << cx << ", " << cy << ", " << cz << ")"
              << std::endl;
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // TUMOR_INITIATION_H_
