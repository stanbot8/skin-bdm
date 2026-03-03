#ifndef MECHANOTRANSDUCTION_PARAMS_H_
#define MECHANOTRANSDUCTION_PARAMS_H_

#include "biodynamo.h"

namespace bdm {
namespace skibidy {

struct MechanotransductionParams {

  // Mechanotransduction (tissue stiffness and wound contraction)
  // Tomasek et al. 2002 (doi:10.1038/nrm809)
  bool enabled = true;
};

}  // namespace skibidy
}  // namespace bdm

#endif  // MECHANOTRANSDUCTION_PARAMS_H_
