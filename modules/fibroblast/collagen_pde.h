#ifndef COLLAGEN_PDE_H_
#define COLLAGEN_PDE_H_

#include "core/pde.h"
#include "infra/sim_param.h"

namespace bdm {
namespace skibidy {

struct CollagenPDE : public SimpleStructuralPDE {
  explicit CollagenPDE(const SimParam* sp)
      : SimpleStructuralPDE(fields::kCollagen, fields::kCollagenId,
                            0, sp->collagen_decay) {}
};

}  // namespace skibidy
}  // namespace bdm

#endif  // COLLAGEN_PDE_H_
