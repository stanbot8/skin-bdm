#ifndef SENESCENCE_PDE_H_
#define SENESCENCE_PDE_H_

#include "core/pde.h"

namespace bdm {
namespace skibidy {

// Senescence field: represents local senescent cell density in the wound bed.
// Accumulates from DNA damage (wound-induced), oxidative stress (inflammation),
// and glycative stress (AGE in diabetic mode). Cleared slowly by immune-mediated
// surveillance (NK cells, macrophage phagocytosis).
// Non-diffusing: senescent cells are immobile once growth-arrested.
// Demaria et al. 2014 (doi:10.1016/j.devcel.2014.11.012)
struct SenescencePDE : public SimpleStructuralPDE {
  SenescencePDE(const SimParam* sp)
      : SimpleStructuralPDE(fields::kSenescence, fields::kSenescenceId,
                            sp->senescence.diffusion, sp->senescence.decay) {}
};

}  // namespace skibidy
}  // namespace bdm

#endif  // SENESCENCE_PDE_H_
