#ifndef MMP_PDE_H_
#define MMP_PDE_H_

#include "core/pde.h"
#include "infra/sim_param.h"

namespace bdm {
namespace skibidy {

// MMP PDE with residual-only decay. Bulk MMP clearance is now handled by
// second-order MMP*TIMP neutralization in fused_post.h. The residual decay
// represents non-TIMP clearance (endocytosis, auto-degradation).
struct MMPPDE : public SimplePDE {
  explicit MMPPDE(const SimParam* sp)
      : SimplePDE(fields::kMMP, fields::kMMPId,
                  sp->mmp_diffusion, sp->mmp_residual_decay) {}
};

}  // namespace skibidy
}  // namespace bdm

#endif  // MMP_PDE_H_
