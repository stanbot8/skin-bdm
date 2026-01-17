#ifndef IMMUNE_RESPONSE_H_
#define IMMUNE_RESPONSE_H_

#include <cmath>

#include "immune/immune_cell.h"
#include "immune/macrophage_behavior.h"
#include "immune/neutrophil_behavior.h"
#include "infra/sim_param.h"
#include "infra/util.h"
#include "core/operation/operation.h"

namespace bdm {
namespace skibidy {

// ---------------------------------------------------------------------------
// ImmuneResponse -- standalone operation that spawns immune cells.
// Neutrophils arrive in discrete waves (~2h, 3 waves over 24h).
// Macrophages are recruited continuously: each step after the spawn delay,
// probability of spawning = rate * max(0, infl - threshold) * saturation.
// Saturation = max(0, 1 - n_macrophages / capacity) models limited ICAM-1
// binding sites at postcapillary venules. This prevents runaway feedback
// (M1 inflammation -> recruitment -> more inflammation).
// ---------------------------------------------------------------------------
struct ImmuneResponse : public StandaloneOperationImpl {
  BDM_OP_HEADER(ImmuneResponse);

  void operator()() override {
    auto* sim = Simulation::GetActive();
    auto* scheduler = sim->GetScheduler();
    auto* sp = sim->GetParam()->Get<SimParam>();

    if (!sp->wound_enabled) return;

    uint64_t step = scheduler->GetSimulatedSteps();
    uint64_t wound_step = static_cast<uint64_t>(sp->wound_trigger_step);
    if (step <= wound_step) return;

    uint64_t wound_age = step - wound_step;

    // Neutrophil waves (discrete, fast arrival)
    int eff_waves = sp->neutrophil_spawn_waves;
    int eff_window = sp->neutrophil_spawn_window;
    if (sp->diabetic_mode) {
      eff_waves = static_cast<int>(
          std::ceil(eff_waves * sp->diabetic_neutrophil_waves_factor));
      eff_window = static_cast<int>(
          std::ceil(eff_window * sp->diabetic_neutrophil_window_factor));
    }
    if (neut_waves_spawned_ < eff_waves) {
      uint64_t wave_step = sp->neutrophil_spawn_delay;
      if (eff_waves > 1) {
        wave_step += neut_waves_spawned_ *
                     (eff_window / (eff_waves - 1));
      }
      if (wound_age >= wave_step) {
        int total = ComputePerimeterSlots(sp);
        if (sp->diabetic_mode) {
          total = static_cast<int>(
              std::ceil(total * sp->diabetic_neutrophil_factor));
        }
        int this_wave = total / eff_waves;
        if (neut_waves_spawned_ == eff_waves - 1)
          this_wave = total - this_wave * (eff_waves - 1);
        SpawnCells(sim, sp, kNeutrophil, this_wave);
        neut_waves_spawned_++;
        if (sp->debug_immune) {
          std::cout << "[immune] spawned " << this_wave
                    << " neutrophils (wave " << neut_waves_spawned_
                    << ")" << std::endl;
        }
      }
    }

    // Macrophage continuous recruitment (inflammation-driven)
    if (wound_age >= static_cast<uint64_t>(sp->macrophage_spawn_delay)) {
      real_t query_z = sp->immune_cell_diameter / 2.0;
      Real3 center = {sp->wound_center_x, sp->wound_center_y, query_z};
      Real3 qpos = ClampToBounds(center, sim->GetParam());
      real_t infl = GetNetInflammation(sim, qpos);
      real_t excess = std::max(static_cast<real_t>(0),
                               infl - sp->macrophage_spawn_threshold);

      // Carrying capacity: limited adhesion sites at wound margin
      int n_mac = CountMacrophages(sim);
      int capacity = ComputePerimeterSlots(sp);
      real_t saturation = std::max(static_cast<real_t>(0),
                                   1.0 - static_cast<real_t>(n_mac) / capacity);

      // Chemokine gradient decay: monocyte recruitment diminishes as wound matures
      real_t recruit_taper = 1.0;
      if (sp->macrophage_spawn_taper > 0) {
        real_t wound_age_h = wound_age * sim->GetParam()->simulation_time_step;
        real_t eff_taper = sp->macrophage_spawn_taper;
        if (sp->diabetic_mode) {
          eff_taper *= sp->diabetic_macrophage_taper_factor;
        }
        recruit_taper = std::exp(-eff_taper * wound_age_h);
      }

      real_t prob = sp->macrophage_spawn_rate * excess * saturation * recruit_taper;
      if (prob > 0 && sim->GetRandom()->Uniform(0, 1) < prob) {
        SpawnCells(sim, sp, kMacrophage, 1);
        mac_total_spawned_++;
        if (sp->debug_immune) {
          std::cout << "[immune] macrophage #" << mac_total_spawned_
                    << " infl=" << infl << " sat=" << saturation << std::endl;
        }
      }
    }
  }

 private:
  int neut_waves_spawned_ = 0;
  int mac_total_spawned_ = 0;

  // Perimeter slots: wound margin positions available for immune cell arrival
  int ComputePerimeterSlots(const SimParam* sp) {
    real_t r = sp->wound_radius;
    real_t diameter = sp->immune_cell_diameter;
    int n = static_cast<int>(std::round(2.0 * M_PI * r / diameter));
    if (n < 6) n = 6;
    return n;
  }

  int CountMacrophages(Simulation* sim) {
    int count = 0;
    sim->GetResourceManager()->ForEachAgent([&](Agent* a) {
      auto* ic = dynamic_cast<ImmuneCell*>(a);
      if (ic && ic->GetImmuneCellType() == kMacrophage) count++;
    });
    return count;
  }

  void SpawnCells(Simulation* sim, const SimParam* sp,
                  ImmuneCellType type, int n_cells) {
    auto* ctxt = sim->GetExecutionContext();
    auto* random = sim->GetRandom();

    real_t cx = sp->wound_center_x;
    real_t cy = sp->wound_center_y;
    real_t r = sp->wound_radius;
    real_t diameter = sp->immune_cell_diameter;
    real_t spawn_z = diameter / 2.0;

    real_t angle_step = 2.0 * M_PI / n_cells;
    real_t theta_offset = random->Uniform(0, 2.0 * M_PI);

    for (int i = 0; i < n_cells; i++) {
      real_t theta = theta_offset + i * angle_step;
      real_t x = cx + r * std::cos(theta);
      real_t y = cy + r * std::sin(theta);

      x = std::max(sp->tissue_min + diameter,
                    std::min(x, sp->tissue_max - diameter));
      y = std::max(sp->tissue_min + diameter,
                    std::min(y, sp->tissue_max - diameter));

      auto* cell = new ImmuneCell({x, y, spawn_z});
      cell->SetDiameter(diameter);
      cell->SetImmuneCellType(type);
      cell->SetState(kM1Active);

      if (type == kNeutrophil) {
        cell->AddBehavior(new NeutrophilBehavior());
      } else {
        cell->AddBehavior(new MacrophageBehavior());
      }
      ctxt->AddAgent(cell);
    }
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // IMMUNE_RESPONSE_H_
