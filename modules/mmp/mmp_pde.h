#ifndef MMP_PDE_H_
#define MMP_PDE_H_

#include "core/pde.h"
#include "infra/sim_param.h"

namespace bdm {
namespace skibidy {

struct MMPPDE : public SimplePDE {
  explicit MMPPDE(const SimParam* sp)
      : SimplePDE(fields::kMMP, fields::kMMPId,
                  sp->mmp_diffusion,
                  sp->diabetic_mode ? sp->mmp_decay * sp->diabetic_timp_factor
                                    : sp->mmp_decay) {}
};

}  // namespace skibidy
}  // namespace bdm

#endif  // MMP_PDE_H_
