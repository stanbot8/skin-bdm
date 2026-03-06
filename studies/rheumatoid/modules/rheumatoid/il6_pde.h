#ifndef IL6_PDE_H_
#define IL6_PDE_H_

#include "core/pde.h"
#include "infra/sim_param.h"

namespace bdm {
namespace skibidy {

// IL-6 cytokine field: second major inflammatory axis in rheumatoid arthritis.
// Produced by macrophages, fibroblast-like synoviocytes (FLS), and T cells
// via NF-kB and TNF-alpha stimulation. Drives acute phase response (CRP, SAA),
// B cell differentiation, Th17 polarization, and RANKL-mediated bone erosion.
// Cleared by IL-6R/gp130 receptor endocytosis. Therapeutic target of
// tocilizumab (anti-IL-6R monoclonal antibody).
// Kishimoto 2005 (doi:10.1146/annurev.immunol.23.021704.115806)
struct IL6PDE : public PDE {
  explicit IL6PDE(const SimParam* sp)
      : diffusion_(sp->il6_diffusion),
        decay_(sp->il6_decay) {}

  const char* GetName() const override { return fields::kIL6; }
  int GetId() const override { return fields::kIL6Id; }

  void Init(Simulation* sim) override {
    DefineGrid(sim, diffusion_, decay_);
    // Starts at 0; produced by autoimmune source and TNF-alpha induction
  }

 private:
  real_t diffusion_;
  real_t decay_;
};

}  // namespace skibidy
}  // namespace bdm

#endif  // IL6_PDE_H_
