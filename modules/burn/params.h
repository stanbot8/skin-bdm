#ifndef BURN_PARAMS_H_
#define BURN_PARAMS_H_

#include "biodynamo.h"

namespace bdm {
namespace skibidy {

struct BurnParams {

  // Burn module (Jackson's burn wound model, three concentric zones)
  // Jackson 1953, Herndon 2012 (ISBN 978-1-4377-2786-9)
  bool enabled = false;
  real_t coagulation_necrosis = 0.95;  // irreversible cell death fraction
  real_t stasis_viability = 0.5;  // zone of stasis initial viability
  real_t stasis_deterioration = 0.02;  // ischemic deterioration rate (h^-1)
  real_t hyperemia_boost = 0.3;  // inflammatory enhancement in hyperemia zone
  real_t depth_fraction = 0.5;  // 0-0.3 superficial, 0.3-0.7 partial, 0.7+ full
  real_t eschar_rate = 0.01;  // eschar formation from necrotic tissue
  real_t tewl_multiplier = 5.0;  // transepidermal water loss multiplier
  real_t fluid_loss_rate = 0.05;  // fluid loss through burn surface
  real_t infection_susceptibility = 3.0;  // barrier breach infection multiplier
  real_t contracture_rate = 0.002;  // excessive scar contraction tendency
};

}  // namespace skibidy
}  // namespace bdm

#endif  // BURN_PARAMS_H_
