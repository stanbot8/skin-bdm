#ifndef AGE_PDE_H_
#define AGE_PDE_H_

#include "core/pde.h"
#include "infra/sim_param.h"

namespace bdm {
namespace skibidy {

// AGE (Advanced Glycation End-product) PDE: non-diffusing accumulator field.
// AGEs form via non-enzymatic glycation of proteins (Maillard reaction) in
// hyperglycemic tissue. Once formed, AGEs are essentially irreversible and
// persist in the ECM, driving chronic inflammation via RAGE signaling.
// Brownlee 2005 (doi:10.2337/diabetes.54.6.1615)
struct AGEPDE : public SimplePDE {
  explicit AGEPDE(const SimParam* sp)
      : SimplePDE(fields::kAGE, fields::kAGEId,
                  0.0, sp->age_decay) {}
};

}  // namespace skibidy
}  // namespace bdm

#endif  // AGE_PDE_H_
