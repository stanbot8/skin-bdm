#ifndef TUMOR_PDE_H_
#define TUMOR_PDE_H_

#include "core/pde.h"

namespace bdm {
namespace skibidy {

struct TumorPDE : public SimpleStructuralPDE {
  TumorPDE() : SimpleStructuralPDE(fields::kTumor, fields::kTumorId, 0, 0) {}
};

}  // namespace skibidy
}  // namespace bdm

#endif  // TUMOR_PDE_H_
