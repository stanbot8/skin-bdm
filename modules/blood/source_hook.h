#ifndef BLOOD_SOURCE_HOOK_H_
#define BLOOD_SOURCE_HOOK_H_

#include "core/hook_api.h"

namespace bdm {
namespace skibidy {

// Blood module source hook: systemic modifier for perfusion and hemostasis.
// Models three coupled systems:
//   1. Hemorrhage: wound bleeds, reducing blood volume over time
//   2. Coagulation cascade: intrinsic/extrinsic pathways generate thrombin,
//      which enhances fibrin deposition and platelet aggregation
//   3. Blood volume: volume loss reduces perfusion (hemorrhagic shock)
// Versteeg et al. 2013 (doi:10.1152/physrev.00016.2011)
// Guo & DiPietro 2010 (doi:10.1177/0022034509359125)
struct BloodSourceHook {
  DiffusionGrid* perf_grid = nullptr;
  DiffusionGrid* o2_grid = nullptr;
  DiffusionGrid* fibrin_grid = nullptr;
  DiffusionGrid* infl_grid = nullptr;
  const SimParam* sp_ = nullptr;
  bool active = false;

  // Runtime state (evolved each step)
  real_t blood_volume_ = 1.0;     // current blood volume (normalized)
  real_t coag_state_ = 0.0;      // coagulation cascade activation (0 to 1)
  real_t thrombin_ = 0.0;        // thrombin concentration (0 to 1)
  int wound_age_ = 0;            // steps since wound creation

  inline void Init(const GridRegistry& reg, SignalBoard& sig) {
    sp_ = reg.Params();
    active = sp_->blood.enabled;
    if (!active) return;
    perf_grid = reg.Get(fields::kVascularId);
    o2_grid = reg.Get(fields::kOxygenId);
    if (sp_->hemostasis.enabled)
      fibrin_grid = reg.Get(fields::kFibrinId);
    infl_grid = reg.Get(fields::kInflammationId);
    blood_volume_ = sp_->blood.initial_volume;
    coag_state_ = 0.0;
    thrombin_ = 0.0;
    wound_age_ = 0;
  }

  // Dermal: apply systemic perfusion modifier based on blood volume and hemoglobin.
  // Blood volume loss reduces effective perfusion throughout wound.
  inline void ApplyDermal(const VoxelSnapshot& snap, SignalBoard& sig) {
    // Perfusion penalty from blood volume loss (hemorrhagic shock)
    real_t perf_factor = 1.0;
    if (blood_volume_ < sp_->blood.shock_threshold) {
      // Below shock threshold: sharp perfusion drop
      real_t shock_depth = (sp_->blood.shock_threshold - blood_volume_) /
                            sp_->blood.shock_threshold;
      perf_factor = 1.0 - shock_depth * (1.0 - sp_->blood.shock_perfusion_penalty);
      perf_factor = std::max(perf_factor, sp_->blood.shock_perfusion_penalty);
    }

    // Anemia penalty on O2 delivery
    real_t hb = sp_->blood.hemoglobin;
    if (hb < sp_->blood.anemia_threshold && o2_grid) {
      real_t anemia_depth = (sp_->blood.anemia_threshold - hb) /
                             sp_->blood.anemia_threshold;
      real_t o2_penalty = anemia_depth * sp_->blood.anemia_o2_penalty;
      real_t o2_val = o2_grid->GetConcentration(snap.idx);
      if (o2_val > 0) {
        o2_grid->ChangeConcentrationBy(snap.idx, -o2_val * o2_penalty * 0.01);
      }
    }

    // Apply perfusion modifier (only in wound where blood loss matters)
    if (snap.in_wound && perf_factor < 1.0 && perf_grid) {
      real_t current = perf_grid->GetConcentration(snap.idx);
      real_t delta = current * (perf_factor - 1.0) * 0.01;
      if (std::abs(delta) > 1e-10) {
        perf_grid->ChangeConcentrationBy(snap.idx, delta);
      }
    }
  }

  // Epidermal wound: hemorrhage, coagulation cascade, platelet cytokines.
  inline void ApplyEpiWound(const VoxelSnapshot& snap, SignalBoard& sig) {
    wound_age_++;

    // Hemorrhage: blood loss proportional to wound exposure
    real_t bleed = sp_->blood.bleed_rate * (1.0 - snap.barrier);
    // Anticoagulant therapy increases bleeding
    if (sp_->blood.anticoagulant_level > 0) {
      bleed *= (1.0 + (sp_->blood.anticoag_bleed_boost - 1.0) *
                sp_->blood.anticoagulant_level);
    }
    // Vascularity coupling
    if (sp_->blood.vascularity_coupling) {
      bleed *= snap.vasc;
    }
    // Volume loss from hemorrhage
    real_t dt = snap.dt;
    blood_volume_ -= bleed * sp_->blood.volume_loss_rate * dt;
    // Volume recovery (physiological compensation)
    if (blood_volume_ < sp_->blood.initial_volume) {
      blood_volume_ += sp_->blood.volume_recovery_rate * dt;
      blood_volume_ = std::min(blood_volume_, sp_->blood.initial_volume);
    }
    blood_volume_ = std::max(blood_volume_, static_cast<real_t>(0.1));

    // Coagulation cascade: intrinsic + extrinsic -> thrombin
    real_t coag_rate = sp_->blood.intrinsic_rate + sp_->blood.extrinsic_rate;
    // Anticoagulant therapy slows coagulation
    if (sp_->blood.anticoagulant_level > 0) {
      coag_rate *= (1.0 - sp_->blood.anticoag_coag_penalty *
                    sp_->blood.anticoagulant_level);
    }
    // Thrombocytopenia impairs clotting
    if (sp_->blood.platelet_count < sp_->blood.thrombocytopenia_threshold) {
      coag_rate *= sp_->blood.platelet_count /
                   sp_->blood.thrombocytopenia_threshold;
    }
    coag_state_ += coag_rate * (1.0 - coag_state_) * dt;
    coag_state_ = std::min(coag_state_, static_cast<real_t>(1.0));

    // Thrombin generation from converged pathways
    thrombin_ += sp_->blood.thrombin_rate * coag_state_ * (1.0 - thrombin_) * dt;
    thrombin_ = std::min(thrombin_, static_cast<real_t>(1.0));

    // Thrombin -> fibrin coupling (enhances hemostasis module)
    if (fibrin_grid && thrombin_ > 0.1) {
      real_t fibrin_boost = sp_->blood.thrombin_fibrin_coupling * thrombin_ * dt;
      fibrin_grid->ChangeConcentrationBy(snap.idx, fibrin_boost);
    }

    // Platelet cytokine release (PDGF, TGF-beta) -> inflammation
    if (infl_grid && coag_state_ > 0.1) {
      real_t cytokine = sp_->blood.platelet_cytokine * coag_state_ *
                        sp_->blood.platelet_count * dt;
      if (cytokine > 1e-10) {
        infl_grid->ChangeConcentrationBy(snap.idx, cytokine);
      }
    }
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // BLOOD_SOURCE_HOOK_H_
