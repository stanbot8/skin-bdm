#ifndef TEMPERATURE_PARAMS_H_
#define TEMPERATURE_PARAMS_H_

#include "biodynamo.h"

namespace bdm {
namespace skibidy {

struct TemperatureParams {

  // Temperature (wound thermal regulation)
  // Fierheller & Sibbald 2010 (doi:10.1097/01.ASW.0000363527.76169.b2)
  bool enabled = true;
  real_t diffusion = 0.001;  // thermal conduction
  real_t decay = 0.0;  // no intrinsic decay
  real_t wound_surface = 0.946;  // ~35C/37C normalized
  real_t perfusion_warming = 0.025;  // rewarming rate per step
  real_t surface_cooling = 0.002;  // evaporative heat loss
  real_t q10_migration = 1.0;  // cell migration Q10
  real_t q10_proliferation = 1.0;  // cell cycle Q10 (disabled normal; active diabetic)
  real_t q10_mmp = 1.3;  // MMP enzymatic Q10
  real_t q10_biofilm = 2.0;  // bacterial growth Q10
};

}  // namespace skibidy
}  // namespace bdm

#endif  // TEMPERATURE_PARAMS_H_
