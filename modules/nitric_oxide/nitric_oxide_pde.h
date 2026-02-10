#ifndef NITRIC_OXIDE_PDE_H_
#define NITRIC_OXIDE_PDE_H_

#include "core/pde.h"
#include "infra/sim_param.h"

namespace bdm {
namespace skibidy {

// Nitric oxide field: produced by iNOS in M1 macrophages and neutrophils.
// Very short half-life (seconds to minutes in tissue, modeled as high decay).
// Functions: vasodilation (increases local perfusion), antimicrobial
// (kills bacteria in biofilm), and anti-fibrotic (suppresses excessive
// collagen deposition).
struct NitricOxidePDE : public PDE {
  explicit NitricOxidePDE(const SimParam* sp)
      : diffusion_(sp->no_diffusion),
        decay_(sp->no_decay) {}

  const char* GetName() const override { return fields::kNitricOxide; }
  int GetId() const override { return fields::kNitricOxideId; }

  void Init(Simulation* sim) override {
    DefineGrid(sim, diffusion_, decay_);
    // Starts at 0; produced by immune cells after wounding
  }

 private:
  real_t diffusion_;
  real_t decay_;
};

}  // namespace skibidy
}  // namespace bdm

#endif  // NITRIC_OXIDE_PDE_H_
