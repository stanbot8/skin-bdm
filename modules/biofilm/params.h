#ifndef BIOFILM_PARAMS_H_
#define BIOFILM_PARAMS_H_

#include "biodynamo.h"

namespace bdm {
namespace skibidy {

struct BiofilmParams {

  // Biofilm dynamics (Extension 7)
  bool enabled = false;
  real_t growth_rate = 0.005;  // per-step logistic growth
  real_t carrying_capacity = 1.0;  // max biofilm density per voxel
  int seed_delay = 480;  // TOML: 48h post-wound
  real_t seed_amount = 0.01;  // initial inoculum density
  real_t neutrophil_clearance = 0.003;  // clearance per neutrophil per step
  real_t macrophage_clearance = 0.001;  // clearance per macrophage per step
  real_t m1_block_threshold = 0.3;  // biofilm level blocking M1->M2
  real_t inflammation_rate = 0.001;  // PAMP-driven inflammation per biofilm
};

}  // namespace skibidy
}  // namespace bdm

#endif  // BIOFILM_PARAMS_H_
