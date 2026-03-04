#ifndef PROMMP_PDE_H_
#define PROMMP_PDE_H_

#include "core/pde.h"
#include "infra/sim_param.h"

namespace bdm {
namespace skibidy {

// Pro-MMP (zymogen) PDE: secreted as inactive zymogens by immune cells and
// fibroblasts. Activated to MMP by plasmin and autocatalytic MMP-3 cleavage.
// Visse & Nagase 2003 (doi:10.1161/01.res.0000070112.80711.3d)
struct ProMMPPDE : public SimplePDE {
  explicit ProMMPPDE(const SimParam* sp)
      : SimplePDE(fields::kProMMP, fields::kProMMPId,
                  sp->mmp.diffusion, sp->mmp.prommp_decay) {}
};

}  // namespace skibidy
}  // namespace bdm

#endif  // PROMMP_PDE_H_
