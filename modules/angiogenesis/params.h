#ifndef ANGIOGENESIS_PARAMS_H_
#define ANGIOGENESIS_PARAMS_H_

#include "biodynamo.h"

namespace bdm {
namespace skibidy {

struct AngiogenesisParams {
  real_t vegf_migration_weight = 0.5;  // relative weight of VEGF chemotaxis (vs O2=1.0)

  // VEGF-driven angiogenesis (Extension 8)
  bool enabled = true;
  real_t vegf_diffusion = 0.08;  // VEGF diffusion coefficient
  real_t vegf_decay = 0.01;  // VEGF natural decay rate
  real_t vegf_production_rate = 0.002;  // hypoxia-driven VEGF production
  real_t vegf_hypoxia_threshold = 0.5;  // O2 level below which VEGF is produced
  real_t vegf_consumption_rate = 0.1;  // VEGF consumed per unit perfusion recovery
  real_t vegf_production_taper = 0.001;  // PHD2/3 negative feedback: exp(-k*wound_age) on VEGF production
  real_t vegf_receptor_clearance = 0.005;  // VEGFR-2 internalization rate per vascular density (Ferrara 2003)
};

}  // namespace skibidy
}  // namespace bdm

#endif  // ANGIOGENESIS_PARAMS_H_
