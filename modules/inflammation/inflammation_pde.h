#ifndef INFLAMMATION_PDE_H_
#define INFLAMMATION_PDE_H_

#include "core/pde.h"
#include "infra/sim_param.h"

namespace bdm {
namespace skibidy {

struct InflammationPDE : public SimplePDE {
  explicit InflammationPDE(const SimParam* sp)
      : SimplePDE(fields::kInflammation, fields::kInflammationId,
                  sp->inflammation_diffusion, sp->inflammation_decay) {}
};

struct ProInflammatoryPDE : public SimplePDE {
  explicit ProInflammatoryPDE(const SimParam* sp)
      : SimplePDE(fields::kProInflammatory, fields::kProInflammatoryId,
                  sp->inflammation_diffusion, sp->inflammation_decay) {}
};

struct AntiInflammatoryPDE : public SimplePDE {
  explicit AntiInflammatoryPDE(const SimParam* sp)
      : SimplePDE(fields::kAntiInflammatory, fields::kAntiInflammatoryId,
                  sp->anti_inflammation_diffusion, sp->anti_inflammation_decay) {}
};

// Immune pressure: mirrors inflammation dynamics (same D, mu) but excludes
// wound-derived DAMPs. Used by keratinocytes for proliferation/migration
// suppression so that tissue damage signals don't create a feedback loop.
struct ImmunePressurePDE : public SimplePDE {
  explicit ImmunePressurePDE(const SimParam* sp)
      : SimplePDE(fields::kImmunePressure, fields::kImmunePressureId,
                  sp->inflammation_diffusion, sp->inflammation_decay) {}
};

}  // namespace skibidy
}  // namespace bdm

#endif  // INFLAMMATION_PDE_H_
