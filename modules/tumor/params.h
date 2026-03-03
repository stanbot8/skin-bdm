#ifndef TUMOR_PARAMS_H_
#define TUMOR_PARAMS_H_

#include "biodynamo.h"

namespace bdm {
namespace skibidy {

struct TumorParams {

  // Tumor module
  bool enabled = false;
  int seed_time = 0;  // TOML: hours; 0h = frame 0
  real_t seed_x = -5.0;  // cluster center X (offset from wound)
  real_t seed_y = -5.0;  // cluster center Y (offset from wound)
  real_t seed_z = 2.0;  // cluster center Z (basal layer)
  int seed_count = 5;  // initial cluster size
  real_t diameter = 4.0;  // cell diameter
  real_t cycle_factor = 3.85;  // all-phase duration multiplier (calibrated for Td~148d)
  real_t g1_factor = 3.85;  // G1 extension factor (equal to cycle_factor; Ki-67 set by contact inhibition)
  int max_neighbors = 12;  // contact inhibition threshold (sphere packing max ~12-14)
  real_t ci_steepness = 2.0;  // soft CI exponent (0 = hard threshold; 2 = quadratic YAP/TAZ)
  real_t growth_rate = 5;  // volume growth rate (same as normal)
  int max_cells = 5000;  // carrying capacity (0 = unlimited)
  int handoff_delay = 5000;  // TOML: 500h in G0 before field conversion (0 = disabled)
  real_t stratum_value = 10.0;  // value written to Stratum field on handoff
  real_t apoptosis_rate = 0.0003;  // per-step probability (calibrated with cycle for Td~148d)
  real_t o2_threshold = 0.05;  // BCC hypoxia tolerance via HIF-1alpha (Gruber et al. 2004)
};

}  // namespace skibidy
}  // namespace bdm

#endif  // TUMOR_PARAMS_H_
