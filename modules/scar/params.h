#ifndef SCAR_PARAMS_H_
#define SCAR_PARAMS_H_

#include "biodynamo.h"

namespace bdm {
namespace skibidy {

struct ScarParams {

  // Inflammation-proportional scarring (Extension 1)
  bool proportional_enabled = false;  // cumulative inflammation integral
  real_t accumulation_rate = 0.001;  // scar += inflammation * rate per step
  real_t collagen_threshold = 0.003;  // local collagen above this drives scar

  // Scar maturation kinetics: 0 immature granulation (cellular, active
  // myofibroblasts, red raised tissue), 1 mature fibrous scar (hypocellular,
  // cross-linked collagen, pale flat tissue). Maturation proceeds via
  // MMP-mediated remodeling and time, slowed by persistent inflammation
  // and active myofibroblasts.
  // Bond JS 2008 (Plast Reconstr Surg, doi:10.1097/PRS.0b013e3181858fa3):
  //   scar maturation timeline spans months in normal healing.
  // Marshall CD 2018 (Adv Wound Care, doi:10.1089/wound.2016.0696):
  //   abnormal remodeling in hypertrophic and keloid lesions.
  bool maturity_enabled = false;          // opt-in; requires ScarMaturityPDE
  real_t maturation_rate = 0.00015;       // base rate per step
  real_t maturation_mmp_boost = 4.0;      // MMP accelerates cross-link turnover
  real_t maturation_infl_block = 3.0;     // inflammation slows maturation
  real_t maturation_myofib_block = 2.0;   // active myofib keep scar immature
};

}  // namespace skibidy
}  // namespace bdm

#endif  // SCAR_PARAMS_H_
