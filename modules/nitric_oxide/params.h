#ifndef NITRIC_OXIDE_PARAMS_H_
#define NITRIC_OXIDE_PARAMS_H_

#include "biodynamo.h"

namespace bdm {
namespace skibidy {

struct NitricOxideParams {

  // Nitric oxide (immune antimicrobial and vasodilator)
  // Schaffer et al. 1996 (doi:10.1016/S0022-4804(96)80068-5)
  bool enabled = true;
};

}  // namespace skibidy
}  // namespace bdm

#endif  // NITRIC_OXIDE_PARAMS_H_
