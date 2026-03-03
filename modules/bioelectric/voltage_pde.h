#ifndef VOLTAGE_PDE_H_
#define VOLTAGE_PDE_H_

#include "core/pde.h"
#include "infra/sim_param.h"

namespace bdm {
namespace skibidy {

// Transepithelial potential (TEP) field: maintained by Na+/K+ ATPase in
// intact epithelium, collapses at wound edge creating lateral voltage
// gradient that drives galvanotaxis. Fast diffusion models ionic current
// flow through extracellular fluid.
// Zhao et al. 2006 (doi:10.1038/nature04925)
struct VoltagePDE : public PDE {
  explicit VoltagePDE(const SimParam* sp)
      : diffusion_(sp->voltage_diffusion),
        decay_(sp->voltage_decay) {}

  const char* GetName() const override { return fields::kVoltage; }
  int GetId() const override { return fields::kVoltageId; }

  void Init(Simulation* sim) override {
    // Fast diffusion: ionic current spreads quickly through tissue fluid.
    // Decay: charge leakage through tissue resistance.
    DefineGrid(sim, diffusion_, decay_);
  }

 private:
  real_t diffusion_;
  real_t decay_;
};

}  // namespace skibidy
}  // namespace bdm

#endif  // VOLTAGE_PDE_H_
