#ifndef FIBROBLAST_PARAMS_H_
#define FIBROBLAST_PARAMS_H_

#include "biodynamo.h"

namespace bdm {
namespace skibidy {

struct FibroblastParams {

  // Fibroblast / TGF-beta / Collagen module (TOML: hours -> steps internally)
  bool enabled = true;
  int spawn_delay = 10;  // 1h post-wound
  int spawn_waves = 4;  // recruitment pulses
  int spawn_window = 960;  // 96h window
  real_t diameter = 4.0;
  real_t density_factor = 1.0;
  int activation_delay = 100;  // 10h before auto-activation
  real_t activation_threshold = 0.1;
  int myofibroblast_delay = 960;  // 96h min in activated
  real_t myofibroblast_threshold = 0.3;
  real_t apoptosis_threshold = 0.003;
  int apoptosis_onset = 1920;  // 192h as myofibroblast
  real_t apoptosis_rate = 0.0008;  // per-step probability
  int min_lifespan = 1680;  // 168h (7d)
  int lifespan = 6720;  // 672h (28d)
  real_t migration_speed = 1.0;
  real_t tgfb_diffusion = 0.03;
  real_t tgfb_decay = 0.005;
  real_t tgfb_wound_seed = 0.15;  // platelet alpha-granule TGF-beta bolus (Shah et al. 1995)
  real_t tgfb_rate = 0.005;
  real_t collagen_deposition_rate = 0.002;
  real_t collagen_decay = 0.0;
  real_t tgfb_receptor_consumption = 0.0;  // per-cell receptor endocytosis (Vilar 2006)
  real_t tgfb_tissue_clearance = 0.0;  // tissue receptor clearance scaled by local cell density
  real_t collagen_o2_half_max = 0.0;  // O2 for half-max prolyl hydroxylation (Myllyharju 2003)
  real_t collagen_fn_transition = 0.15;  // collagen level at which FN deposition ceases (ECM maturation)
  real_t tgfb_taper = 0.0012;  // myofibroblast TGF-beta exp(-k*state_age)

  // Activated fibroblast apoptosis
  int activated_apoptosis_onset = 1680;  // 168h as activated before stochastic death
  real_t activated_apoptosis_rate = 0.0003;  // per-step removal for activated fibroblasts
};

}  // namespace skibidy
}  // namespace bdm

#endif  // FIBROBLAST_PARAMS_H_
