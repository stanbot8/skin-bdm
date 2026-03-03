#ifndef BIOELECTRIC_PARAMS_H_
#define BIOELECTRIC_PARAMS_H_

#include "biodynamo.h"

namespace bdm {
namespace skibidy {

struct BioelectricParams {

  // Bioelectric wound currents (transepithelial potential and galvanotaxis)
  // Zhao et al. 2006 (doi:10.1038/nature04925)
  bool enabled = true;
};

}  // namespace skibidy
}  // namespace bdm

#endif  // BIOELECTRIC_PARAMS_H_
