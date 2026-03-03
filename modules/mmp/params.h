#ifndef MMP_PARAMS_H_
#define MMP_PARAMS_H_

#include "biodynamo.h"

namespace bdm {
namespace skibidy {

struct MmpParams {

  // MMP dynamics (matrix metalloproteinase -- ECM remodeling)
  // Lobmann et al. 2002 (doi:10.1007/s00125-002-0868-8): MMP-1 65x, MMP-9 14x in diabetic wounds
  bool enabled = true;
  real_t diffusion = 0.02;  // MMP diffusion coefficient
  real_t m1_rate = 0.003;  // MMP-9 from M1 macrophages per step
  real_t neutrophil_rate = 0.002;  // MMP-8 from neutrophils per step (Nwomeh 1998)
  real_t fibroblast_rate = 0.001;  // MMP-1/3 from activated fibroblasts per step
  real_t keratinocyte_rate = 0.0005;  // MMP-1 from wound-edge keratinocytes per step
  real_t collagen_degradation = 0.005;  // collagen loss per MMP unit per step
  real_t fibronectin_degradation = 0.002;  // fibronectin loss per MMP unit per step
  real_t residual_decay = 0.01;  // non-TIMP clearance (endocytosis, auto-degradation)
  // TIMP dynamics (tissue inhibitor of metalloproteinases)
  // Brew et al. 2000 (doi:10.1016/S0167-4838(99)00252-5)
  real_t timp_diffusion = 0.02;  // TIMP diffusion (similar MW to MMP)
  real_t timp_decay = 0.005;  // TIMP natural degradation (more stable than MMP)
  real_t timp_fibroblast_rate = 0.002;  // TIMP-1 from activated fibroblasts (primary source)
  real_t timp_macrophage_rate = 0.001;  // TIMP-1 from M2 macrophages (pro-resolution)
  real_t timp_keratinocyte_rate = 0.0005;  // TIMP-1 from wound-edge keratinocytes
  real_t timp_inhibition_rate = 0.08;  // second-order MMP*TIMP neutralization
  // Pro-MMP zymogen activation cascade
  // Visse & Nagase 2003 (doi:10.1161/01.res.0000070112.80711.3d): MMPs secreted as zymogens,
  // activated by plasmin cleavage and MMP-3 autocatalysis
  real_t prommp_decay = 0.002;  // pro-MMP natural degradation
  real_t prommp_activation_rate = 0.15;  // basal plasmin-mediated activation rate
  real_t prommp_autocatalytic_rate = 0.02;  // MMP-3 autocatalytic activation (rate * proMMP * MMP)
};

}  // namespace skibidy
}  // namespace bdm

#endif  // MMP_PARAMS_H_
