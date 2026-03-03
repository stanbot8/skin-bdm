#ifndef PERFUSION_PARAMS_H_
#define PERFUSION_PARAMS_H_

#include "biodynamo.h"

namespace bdm {
namespace skibidy {

struct PerfusionParams {

  // Vascular perfusion field
  real_t diffusion = 0.03;  // capillary sprouting / lateral spread
  real_t decay = 0.0;  // no decay (stable vessels)
  real_t basal = 1.0;  // healthy dermis = fully perfused
  real_t angio_rate = 0.005;  // recovery speed per step
  int angio_delay = 240;  // TOML: 24h before sprouting begins

  // Per-layer perfusion baselines (fraction of perfusion_basal)
  // Braverman 2000 (doi:10.1046/j.1087-0024.2000.00010.x)
  real_t papillary_fraction = 1.0;  // dense capillary loops
  real_t reticular_fraction = 0.7;  // arterioles and venules
  real_t hypodermis_fraction = 0.3;  // sparse large vessels
  real_t clearance_rate = 0.0;  // venous clearance of soluble TGF-beta via blood flow (Bramhall 1999)
};

}  // namespace skibidy
}  // namespace bdm

#endif  // PERFUSION_PARAMS_H_
