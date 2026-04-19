#ifndef SCAR_MATURITY_PDE_H_
#define SCAR_MATURITY_PDE_H_

#include "core/pde.h"

namespace bdm {
namespace skibidy {

// Scar maturity: 0 = immature granulation (red, raised, cellular, active
// myofibroblasts), 1 = mature fibrous tissue (pale, flat, hypocellular,
// stable). Non-diffusing structural field; evolved mechanistically in
// the scar maturation hook.
struct ScarMaturityPDE : public SimpleStructuralPDE {
  ScarMaturityPDE()
      : SimpleStructuralPDE(fields::kScarMaturity,
                            fields::kScarMaturityId, 0, 0) {}
};

}  // namespace skibidy
}  // namespace bdm

#endif  // SCAR_MATURITY_PDE_H_
