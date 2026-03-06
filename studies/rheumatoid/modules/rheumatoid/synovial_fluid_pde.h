#ifndef SYNOVIAL_FLUID_PDE_H_
#define SYNOVIAL_FLUID_PDE_H_

#include "core/pde.h"
#include "infra/sim_param.h"

namespace bdm {
namespace skibidy {

// Synovial fluid field: represents pannus tissue and synovial hyperplasia.
// In RA, the synovial membrane thickens from 1-2 cell layers to a dense
// pannus of fibroblast-like synoviocytes (FLS) and infiltrating immune cells.
// This field tracks the pannus density: 0 = healthy synovium, 1.0 = fully
// developed pannus. Pannus growth is driven by TNF-alpha and IL-6, and
// amplifies local cytokine production and cartilage erosion.
// Non-diffusing (structural tissue, not soluble factor).
// Smolen et al. 2018 (doi:10.1038/nrdp.2018.1)
// Firestein 2003 (doi:10.1038/nature01661)
struct SynovialFluidPDE : public PDE {
  explicit SynovialFluidPDE(const SimParam* sp)
      : basal_(sp->ra.synovial_basal) {}

  const char* GetName() const override { return fields::kSynovialFluid; }
  int GetId() const override { return fields::kSynovialFluidId; }

  void Init(Simulation* sim) override {
    DefineStructuralGrid(sim, 0.0, 0.0);  // D=0, no decay (structural)
    // Initialize with baseline synovial lining in dermal wound voxels
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

#endif  // SYNOVIAL_FLUID_PDE_H_
