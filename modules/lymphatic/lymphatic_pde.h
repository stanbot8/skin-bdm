#ifndef LYMPHATIC_PDE_H_
#define LYMPHATIC_PDE_H_

#include "core/pde.h"
#include "infra/sim_param.h"

namespace bdm {
namespace skibidy {

// Lymphatic vessel density field: regenerates via VEGF-driven
// lymphangiogenesis after wound disruption. Slow diffusion models
// lymphatic sprouting from wound margins. Provides drainage capacity
// for edema resolution and cytokine clearance.
// Kataru et al. 2009 (doi:10.1182/blood-2008-09-176776)
struct LymphaticPDE : public PDE {
  explicit LymphaticPDE(const SimParam* sp)
      : diffusion_(sp->lymphatic.diffusion) {}

  const char* GetName() const override { return fields::kLymphatic; }
  int GetId() const override { return fields::kLymphaticId; }

  void Init(Simulation* sim) override {
    // Slow diffusion models lymphatic sprouting; no passive decay
    // (lymphatic vessels are structural, cleared only by MMP degradation).
    DefineGrid(sim, diffusion_, 0.0);
  }

 private:
  real_t diffusion_;
};

}  // namespace skibidy
}  // namespace bdm

#endif  // LYMPHATIC_PDE_H_
