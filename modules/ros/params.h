#ifndef ROS_PARAMS_H_
#define ROS_PARAMS_H_

#include "biodynamo.h"

namespace bdm {
namespace skibidy {

struct RosParams {

  // ROS (reactive oxygen species / oxidative stress)
  // Schafer & Werner 2008 (doi:10.1016/j.phrs.2008.06.004)
  bool enabled = true;
  real_t diffusion = 0.03;  // fast diffusion (small molecule)
  real_t decay = 0.04;  // enzymatic scavenging (SOD, catalase, GPx)
  real_t neutrophil_burst = 0.015;  // NADPH oxidase respiratory burst
  real_t m1_burst = 0.008;  // M1 macrophage oxidative burst
  real_t mitochondrial_rate = 0.001;  // basal mitochondrial ETC leak
  real_t hypoxia_threshold = 0.4;  // O2 below this amplifies mito ROS
  real_t hypoxia_amplification = 2.0;  // fold increase under hypoxia
  real_t senescence_rate = 0.0002;  // DNA damage driving senescence
  real_t mmp_activation = 0.005;  // oxidative pro-MMP cysteine switch
  real_t angiogenesis_impairment = 0.1;  // endothelial dysfunction
  real_t inflammation_amplification = 0.0005;  // NF-kB activation
  real_t collagen_damage = 0.001;  // oxidative cross-link damage
  real_t tissue_antioxidant = 0.02;  // enzymatic clearance by resident cells
  real_t perfusion_clearance = 0.005;  // vascular washout
};

}  // namespace skibidy
}  // namespace bdm

#endif  // ROS_PARAMS_H_
