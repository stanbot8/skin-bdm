#ifndef SCAR_PDE_H_
#define SCAR_PDE_H_

#include "core/pde.h"

namespace bdm {
namespace skibidy {

struct ScarPDE : public SimpleStructuralPDE {
  ScarPDE() : SimpleStructuralPDE(fields::kScar, fields::kScarId, 0, 0) {}
};

}  // namespace skibidy
}  // namespace bdm

#endif  // SCAR_PDE_H_
