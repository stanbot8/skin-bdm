#ifndef IMMUNE_PARAMS_H_
#define IMMUNE_PARAMS_H_

#include "biodynamo.h"

namespace bdm {
namespace skibidy {

struct ImmuneParams {

  // Immune response -- neutrophils (TOML: hours -> steps internally)
  int neutrophil_spawn_delay = 20;  // 2h (Kim et al. 2008)
  int neutrophil_spawn_waves = 3;  // recruitment waves (gradual rise)
  int neutrophil_spawn_window = 240;  // 24h between first and last wave
  int neutrophil_lifespan = 480;  // 48h hard ceiling (Wilgus et al. 2013)
  int neutrophil_min_survival = 120;  // 12h before apoptosis possible
  real_t neutrophil_apoptosis_rate = 0.0029;  // half-life ~24h (Wilgus et al. 2013, doi:10.1089/wound.2012.0383)

  // Immune response -- macrophages (TOML: hours -> steps internally)
  int macrophage_spawn_delay = 240;  // 24h (Rodero & Khosrotehrani 2010)
  real_t macrophage_spawn_threshold = 0.1;  // min inflammation for monocyte extravasation
  real_t macrophage_spawn_rate = 0.1;  // per-step spawn probability per unit excess inflammation
  real_t macrophage_spawn_taper = 0.0;  // per-hour decay of recruitment probability (chemokine gradient decline)
  int macrophage_m1_duration = 480;  // 48h in M1 (Krzyszczyk et al. 2018)
  int macrophage_lifespan = 6720;  // 672h (28d) hard ceiling (Lucas et al. 2010)
  int macrophage_min_survival = 240;  // 24h before apoptosis possible
  real_t macrophage_apoptosis_rate = 0.00072;  // half-life ~4d (Lucas et al. 2010, doi:10.4049/jimmunol.0903356)
  real_t macrophage_emigration_rate = 0.02;  // M2 exit rate once stratum re-epithelialized
  real_t cytokine_rate = 0.002;  // inflammation added per cell per step
  real_t resolution_rate = 0.005;  // inflammation consumed by M2 per step
  real_t migration_speed = 1.5;  // tractor force magnitude

  // Cytokine taper rates (configurable exponential decay constants)
  real_t neutrophil_cytokine_taper = 0.002;  // neutrophil cytokine exp(-k*age)
};

}  // namespace skibidy
}  // namespace bdm

#endif  // IMMUNE_PARAMS_H_
