#ifndef PRESSURE_PARAMS_H_
#define PRESSURE_PARAMS_H_

#include "biodynamo.h"

namespace bdm {
namespace skibidy {

struct PressureParams {

  // Pressure ulcer module (ischemia-reperfusion injury)
  // Gefen 2007 (doi:10.12968/jowc.2007.16.8.27854)
  bool enabled = false;
  real_t compression_threshold = 0.3;  // compression above this occludes capillaries
  real_t ischemia_rate = 0.02;  // tissue damage accumulation under compression
  real_t reperfusion_ros_burst = 0.15;  // ROS burst on pressure release
  real_t shear_factor = 1.5;  // shear stress amplifies compression damage
  real_t tissue_damage_rate = 0.003;  // general tissue damage accumulation
  real_t reposition_interval_h = 0.0;  // hours between pressure relief (0 = constant)
  real_t moisture_damage_rate = 0.001;  // moisture-associated skin damage
  // Stage thresholds (tissue damage accumulator -> NPUAP staging)
  real_t stage_1 = 0.15;  // non-blanchable erythema
  real_t stage_2 = 0.35;  // partial thickness skin loss
  real_t stage_3 = 0.60;  // full thickness skin loss
  real_t stage_4 = 0.85;  // full thickness tissue loss (muscle/bone)
};

}  // namespace skibidy
}  // namespace bdm

#endif  // PRESSURE_PARAMS_H_
