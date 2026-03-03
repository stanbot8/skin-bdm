#ifndef BLOOD_PARAMS_H_
#define BLOOD_PARAMS_H_

#include "biodynamo.h"

namespace bdm {
namespace skibidy {

struct BloodParams {

  // Blood module (hemorrhage, coagulation cascade, blood volume effects)
  // Versteeg et al. 2013 (doi:10.1152/physrev.00016.2011)
  bool enabled = false;
  // Hemorrhage
  real_t bleed_rate = 0.02;  // baseline blood loss rate (normalized)
  real_t depth_bleed_factor = 1.5;  // deeper wounds bleed more
  bool vascularity_coupling = true;  // bleed_rate *= local perfusion
  // Coagulation cascade
  real_t intrinsic_rate = 0.01;  // contact activation (factor XII)
  real_t extrinsic_rate = 0.05;  // tissue factor exposure
  real_t thrombin_rate = 0.03;  // thrombin generation from converged pathways
  real_t thrombin_fibrin_coupling = 0.02;  // thrombin to fibrin conversion
  real_t platelet_aggregation = 0.04;  // platelet aggregation at wound
  real_t platelet_cytokine = 0.001;  // PDGF/TGF-beta release to inflammation
  real_t clot_maturation_h = 4.0;  // clot stabilization time (hours)
  real_t fibrinolysis_delay_h = 48.0;  // fibrinolysis onset delay (hours)
  // Blood volume and systemic effects
  real_t initial_volume = 1.0;  // normalized blood volume
  real_t volume_loss_rate = 0.01;  // volume loss per hour at bleed_rate=1
  real_t volume_recovery_rate = 0.005;  // physiological compensation
  real_t shock_threshold = 0.6;  // below this: hemorrhagic shock
  real_t shock_perfusion_penalty = 0.3;  // perfusion multiplier in shock
  // Hemoglobin and oxygen delivery
  real_t hemoglobin = 1.0;  // normalized (1.0 = 14 g/dL)
  real_t anemia_threshold = 0.7;  // below this: impaired O2 delivery
  real_t anemia_o2_penalty = 0.5;  // O2 delivery reduction factor
  // Anticoagulant therapy
  real_t anticoagulant_level = 0.0;  // 0 = none, 1.0 = therapeutic
  real_t anticoag_coag_penalty = 0.4;  // coagulation rate reduction
  real_t anticoag_bleed_boost = 1.8;  // bleed rate increase
  // Platelet count
  real_t platelet_count = 1.0;  // normalized (1.0 = 250k/uL)
  real_t thrombocytopenia_threshold = 0.4;  // below this: impaired clotting
};

}  // namespace skibidy
}  // namespace bdm

#endif  // BLOOD_PARAMS_H_
