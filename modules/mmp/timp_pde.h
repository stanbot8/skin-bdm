#ifndef TIMP_PDE_H_
#define TIMP_PDE_H_

#include "core/pde.h"
#include "infra/sim_param.h"

namespace bdm {
namespace skibidy {

// TIMP (tissue inhibitor of metalloproteinases) PDE channel.
// Explicit TIMP field enables second-order MMP*TIMP neutralization kinetics,
// replacing the implicit TIMP decay constant in the MMP PDE.
// Brew et al. 2000 (doi:10.1016/S0167-4838(99)00252-5)
struct TIMPPDE : public SimplePDE {
  explicit TIMPPDE(const SimParam* sp)
      : SimplePDE(fields::kTIMP, fields::kTIMPId,
                  sp->mmp.timp_diffusion, sp->mmp.timp_decay) {}
};

}  // namespace skibidy
}  // namespace bdm

#endif  // TIMP_PDE_H_
