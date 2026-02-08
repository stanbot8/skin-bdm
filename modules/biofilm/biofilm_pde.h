#ifndef BIOFILM_PDE_H_
#define BIOFILM_PDE_H_

#include "core/pde.h"

namespace bdm {
namespace skibidy {

struct BiofilmPDE : public SimpleStructuralPDE {
  BiofilmPDE() : SimpleStructuralPDE(fields::kBiofilm, fields::kBiofilmId, 0, 0) {}
};

}  // namespace skibidy
}  // namespace bdm

#endif  // BIOFILM_PDE_H_
