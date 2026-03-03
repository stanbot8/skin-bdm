#ifndef DIABETIC_PARAMS_H_
#define DIABETIC_PARAMS_H_

#include "biodynamo.h"

namespace bdm {
namespace skibidy {

struct DiabeticParams {

  // Diabetic / chronic wound model (Extension 5)
  bool mode = false;  // impair M1→M2 transition
  real_t m1_duration_factor = 3.0;  // M1 duration multiplier
  real_t resolution_factor = 0.3;  // M2 resolution rate multiplier
  real_t efferocytosis_factor = 0.5;  // efferocytosis probability multiplier
  real_t prolif_factor = 0.5;  // keratinocyte G1->S multiplier
  real_t migration_factor = 0.45;  // keratinocyte migration speed multiplier
  real_t fibroblast_activation_factor = 2.0;  // activation delay multiplier (>1 = slower)
  real_t collagen_factor = 0.4;  // collagen deposition multiplier
  real_t fibroblast_lifespan_factor = 0.6;  // fibroblast lifespan multiplier
  real_t baseline_inflammation = 0.001;  // per-step AGE/RAGE inflammation (Bierhaus et al. 2005)
  real_t neutrophil_factor = 1.8;  // neutrophil spawn count multiplier
  real_t neutrophil_lifespan_factor = 4.0;  // neutrophil lifespan multiplier (Brem & Tomic-Canic 2007)
  real_t neutrophil_apoptosis_factor = 0.15;  // severely impaired apoptosis (Alba-Loureiro et al. 2007)
  real_t neutrophil_waves_factor = 2.0;  // recruitment wave multiplier
  real_t neutrophil_window_factor = 3.0;  // spawn window multiplier
  real_t macrophage_taper_factor = 0.3;  // recruitment taper multiplier
  real_t macrophage_apoptosis_factor = 0.5;  // impaired macrophage apoptosis (NF-kB anti-apoptotic; Mirza & Koh 2011)
  real_t macrophage_emigration_factor = 0.5;  // M2 emigration rate multiplier (impaired lymphatic clearance)
  real_t inflammation_sensitivity = 2.0;  // keratinocyte inflammation sensitivity multiplier
  real_t migration_infl_K = 0.15;  // Hill half-max for migration inflammation gating (Brownlee 2005)
  real_t prolif_infl_K = 0.20;  // Hill half-max for proliferation inflammation gating (Rasik & Shukla 2000)
  real_t biofilm_clearance_factor = 0.5;  // impaired neutrophil oxidative burst
  real_t vegf_factor = 0.4;  // HIF-1alpha impairment
  real_t ph_recovery_factor = 0.5;  // impaired acid mantle restoration
  real_t mmp_factor = 3.0;  // MMP overexpression multiplier
  real_t timp_production_factor = 0.4;  // reduced TIMP expression (Lobmann et al. 2002)
  real_t no_factor = 0.6;  // impaired iNOS in diabetic tissue
  real_t senescence_factor = 2.5;  // hyperglycemia accelerates senescence
  real_t nerve_factor = 0.3;  // residual nerve density (70% loss)
  real_t nerve_regen_factor = 0.4;  // impaired regeneration
  real_t mito_ros_factor = 3.0;  // hyperglycemia superoxide overproduction
  real_t antioxidant_factor = 0.5;  // impaired SOD/catalase/GPx
  real_t lymphatic_factor = 0.4;  // impaired lymphangiogenesis
  real_t voltage_factor = 0.6;  // impaired ion transport
};

}  // namespace skibidy
}  // namespace bdm

#endif  // DIABETIC_PARAMS_H_
