#ifndef FIBROBLAST_PARAMS_H_
#define FIBROBLAST_PARAMS_H_

#include "biodynamo.h"

namespace bdm {
namespace skibidy {

struct FibroblastParams {

  // Fibroblast / TGF-beta / Collagen module (TOML: hours -> steps internally)
  bool enabled = true;
  int spawn_delay = 240;  // 24h post-wound (Darby et al. 2014)
  int spawn_waves = 6;  // recruitment pulses
  int spawn_window = 1200;  // 120h window (~5 days)
  real_t diameter = 4.0;
  real_t density_factor = 1.0;
  int activation_delay = 600;  // 60h before auto-activation (Darby et al. 2014)
  real_t activation_threshold = 0.005;  // TGF-beta for quiescent->activated (Tomasek et al. 2002)
  int myofibroblast_delay = 480;  // 48h min in activated (~2d; Van De Water et al. 2013)
  real_t myofibroblast_threshold = 0.015;  // TGF-beta for activated->myofibroblast (Tomasek et al. 2002)
  real_t apoptosis_threshold = 0.001;  // TGF-beta below this triggers removal (Hinz 2007)
  int apoptosis_onset = 1680;  // 168h as myofibroblast (~7d; Desmouliere et al. 1995)
  real_t apoptosis_rate = 0.0006;  // per-step probability (Desmouliere et al. 1995)
  int min_lifespan = 1680;  // 168h (7d)
  int lifespan = 12000;  // 1200h (50d; stochastic apoptosis handles natural decline)
  real_t migration_speed = 1.0;
  real_t tgfb_diffusion = 0.03;
  real_t tgfb_decay = 0.0;  // clearance now fully mechanistic (receptor + decorin)
  real_t tgfb_wound_seed = 0.04;  // platelet alpha-granule TGF-beta1 bolus (Shah et al. 1992)
  real_t tgfb_rate = 0.001;  // TGF-beta per myofibroblast per step
  real_t collagen_deposition_rate = 0.0005;  // collagen per myofibroblast per step (Murphy et al. 2012)
  real_t collagen_decay = 0.0;
  real_t tgfb_receptor_consumption = 0.001;  // per-cell TbRII/TbRI endocytosis rate (Vilar et al. 2006)
  real_t tgfb_tissue_clearance = 0.025;  // receptor clearance by resident tissue cells
  real_t collagen_o2_half_max = 0.15;  // O2 for half-max prolyl hydroxylation (Myllyharju 2003)
  real_t collagen_fn_transition = 0.15;  // collagen level at which FN deposition ceases (ECM maturation)
  real_t tgfb_taper = 0.0012;  // myofibroblast TGF-beta exp(-k*state_age)

  // Two-factor myofibroblast apoptosis gate (Jiang et al. 2025,
  // doi:10.1038/s41467-025-58906-z): mechanical unloading AND local
  // IL-1beta signaling jointly trigger apoptosis. Implementation multiplies
  // apoptosis_rate by (1 + scale * local_closure * local_inflammation).
  // Disabled by default: needs multi-replicate calibration. When enabled
  // along with fibromodulin/decorin (collagen proxy here), it replicates
  // the fetal-like rapid activation + expedited clearance pattern that
  // minimizes scarring.
  real_t closure_apoptosis_scale = 0.0;  // combined closure*infl multiplier

  // Collagen-TGF-beta sequestration: decorin/biglycan in mature collagen
  // bind and neutralize TGF-beta (Yamaguchi et al. 1990, doi:10.1038/346281a0).
  // Disabled by default: see above note.
  real_t tgfb_collagen_sequestration = 0.0;  // max TGF-beta sink rate
  real_t tgfb_sequestration_half_max = 0.15;  // collagen for half-max sequestration

  // Activated fibroblast apoptosis
  int activated_apoptosis_onset = 1680;  // 168h as activated before stochastic death
  real_t activated_apoptosis_rate = 0.0003;  // per-step removal for activated fibroblasts
};

}  // namespace skibidy
}  // namespace bdm

#endif  // FIBROBLAST_PARAMS_H_
