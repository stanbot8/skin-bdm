#ifndef CARTILAGE_PDE_H_
#define CARTILAGE_PDE_H_

#include "core/pde.h"
#include "infra/sim_param.h"

namespace bdm {
namespace skibidy {

// Cartilage integrity field: structural tissue representing articular
// cartilage (or in skin context, any tissue substrate targeted by RA).
// Non-diffusing (D=0); initialized at full integrity in dermal voxels.
// Degraded by MMP-mediated proteolysis (aggrecanase/collagenase) and
// direct TNF-alpha-driven osteoclast activation.
// Goldring & Goldring 2007 (doi:10.1002/art.22783)
struct CartilagePDE : public PDE {
  explicit CartilagePDE(const SimParam* sp)
      : basal_(sp->ra.cartilage_basal) {}

  const char* GetName() const override { return fields::kCartilage; }
  int GetId() const override { return fields::kCartilageId; }

  void Init(Simulation* sim) override {
    DefineStructuralGrid(sim, 0.0, 0.0);  // D=0, no decay (structural)
    // Initialize cartilage in dermal voxels (maps to synovial/cartilage zone)
    real_t b = basal_;
    ModelInitializer::InitializeSubstance(GetId(),
        [b](real_t x, real_t y, real_t z) -> real_t {
          return (z < 0) ? b : 0;
        });
  }

 private:
  real_t basal_;
};

}  // namespace skibidy
}  // namespace bdm

#endif  // CARTILAGE_PDE_H_
