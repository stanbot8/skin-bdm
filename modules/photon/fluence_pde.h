#ifndef FLUENCE_PDE_H_
#define FLUENCE_PDE_H_

#include "core/pde.h"
#include "infra/sim_param.h"

namespace bdm {
namespace skibidy {

// Fluence field: photon density (light intensity) in tissue.
// Transport modeled via diffusion approximation with absorption sink.
// D = 1 / (3 * (mu_a + mu_s')), Jacques 2013 (doi:10.1088/0031-9155/58/11/R37)
struct FluencePDE : public PDE {
  explicit FluencePDE(const SimParam* sp)
      : diffusion_(sp->photon.diffusion),
        decay_(sp->photon.decay) {}

  const char* GetName() const override { return fields::kFluence; }
  int GetId() const override { return fields::kFluenceId; }

  void Init(Simulation* sim) override {
    DefineGrid(sim, diffusion_, decay_);
    // Starts at 0; photons injected by PhotonSourceHook beam profile
  }

 private:
  real_t diffusion_;
  real_t decay_;
};

}  // namespace skibidy
}  // namespace bdm

#endif  // FLUENCE_PDE_H_
