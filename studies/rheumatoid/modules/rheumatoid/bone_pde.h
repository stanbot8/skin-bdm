#ifndef BONE_PDE_H_
#define BONE_PDE_H_

#include "core/pde.h"
#include "infra/sim_param.h"

namespace bdm {
namespace skibidy {

// Subchondral bone integrity field: structural tissue representing bone
// beneath articular cartilage. Eroded by RANKL/osteoclast activation
// (driven by IL-6) and direct TNF-alpha osteoclast stimulation.
// Distinct from cartilage degradation: bone erosion occurs deeper and
// via different cellular mechanisms (osteoclasts vs chondrocyte apoptosis).
// McInnes & Schett 2011 (doi:10.1056/NEJMra1004965)
// Schett & Gravallese 2012 (doi:10.1038/nrrheum.2012.153)
struct BonePDE : public PDE {
  explicit BonePDE(const SimParam* sp)
      : basal_(sp->ra.bone_basal),
        depth_z_(sp->ra.bone_depth_z) {}

  const char* GetName() const override { return fields::kBone; }
  int GetId() const override { return fields::kBoneId; }

  void Init(Simulation* sim) override {
    DefineStructuralGrid(sim, 0.0, 0.0);  // D=0, no decay (structural)
    // Initialize bone in deep dermal voxels (below reticular dermis)
    real_t b = basal_;
    real_t dz = depth_z_;
    ModelInitializer::InitializeSubstance(GetId(),
        [b, dz](real_t x, real_t y, real_t z) -> real_t {
          return (z < dz) ? b : 0;
        });
  }

 private:
  real_t basal_;
  real_t depth_z_;
};

}  // namespace skibidy
}  // namespace bdm

#endif  // BONE_PDE_H_
