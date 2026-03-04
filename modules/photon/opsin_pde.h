#ifndef OPSIN_PDE_H_
#define OPSIN_PDE_H_

#include "core/pde.h"
#include "infra/sim_param.h"

namespace bdm {
namespace skibidy {

// Opsin activation field: channelrhodopsin open-state fraction (0 to 1).
// Two-state model: closed <-> open, driven by local fluence.
// Non-diffusing (opsins are membrane-bound proteins), no intrinsic decay
// (dark recovery handled in source hook kinetics).
// Fenno et al. 2011 (doi:10.1146/annurev-neuro-061010-113817)
struct OpsinPDE : public SimpleStructuralPDE {
  explicit OpsinPDE(const SimParam* sp)
      : SimpleStructuralPDE(fields::kOpsin, fields::kOpsinId, 0, 0) {}
};

}  // namespace skibidy
}  // namespace bdm

#endif  // OPSIN_PDE_H_
