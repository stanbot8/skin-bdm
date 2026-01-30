#ifndef VEGF_PDE_H_
#define VEGF_PDE_H_

#include "core/pde.h"
#include "infra/sim_param.h"

namespace bdm {
namespace skibidy {

struct VEGFPDE : public SimplePDE {
  explicit VEGFPDE(const SimParam* sp)
      : SimplePDE(fields::kVEGF, fields::kVEGFId,
                  sp->vegf_diffusion, sp->vegf_decay) {}
};

}  // namespace skibidy
}  // namespace bdm

#endif  // VEGF_PDE_H_
