#ifndef GLUCOSE_MOD_PARAMS_H_
#define GLUCOSE_MOD_PARAMS_H_

#include "biodynamo.h"

namespace bdm {
namespace skibidy {

struct GlucoseParams {

  // Glucose (metabolic substrate and diabetic hyperglycemia)
  // Brem et al. 2007 (doi:10.1016/j.jmb.2007.03.007)
  bool enabled = true;
  real_t diffusion = 0.005;  // interstitial glucose diffusion
  real_t decay = 0.0002;  // cellular uptake
  real_t basal_conc = 1.0;  // normalized healthy level
  real_t perfusion_supply = 0.008;  // perfusion-driven supply rate
  real_t age_rate = 0.0001;  // AGE formation rate from glucose (Brownlee 2005)
  real_t age_inflammation = 0.002;  // AGE-driven inflammation per step (legacy, now via AGE field)
  real_t bacterial_consumption = 0.01;  // biofilm glucose consumption
  real_t prolif_threshold = 0.05;  // min glucose for full proliferation
};

}  // namespace skibidy
}  // namespace bdm

#endif  // GLUCOSE_MOD_PARAMS_H_
