#ifndef FIELDS_STRATUM_H_
#define FIELDS_STRATUM_H_

#include "core/pde.h"
#include "infra/volume.h"

namespace bdm {
namespace skibidy {

// Stratum PDE channel.
// Static scalar field encoding stratum identity at each z-height.
// ParaView renders this as colored volume slabs - the tissue lasagna.
// No source term (agents write directly). Wound zeroes epidermal voxels only.
struct StratumPDE : public PDE {
  const char* GetName() const override { return fields::kStratum; }
  int GetId() const override { return fields::kStratumId; }

  void Init(Simulation* sim) override {
    auto* vm = VolumeManager::Get();
    DefineGrid(sim, 0, 0);
    // GetDermalSubLayer returns kPapillary/kReticular/kHypodermis for z<0
    // (fine-grained) instead of flat kDermis from GetStratumAt.
    ModelInitializer::InitializeSubstance(GetId(),
        [vm](real_t x, real_t y, real_t z) {
          return static_cast<real_t>(vm->GetDermalSubLayer({x, y, z}));
        });
    MarkPrescribed(sim);  // agent-written, no diffusion
  }

  void ApplyWound(Simulation* sim, real_t, real_t, real_t) override {
    auto* grid = Grid(sim);
    auto* sp = sim->GetParam()->Get<SimParam>();
    GridContext ctx(grid, sp);
    GridContext::ZeroWoundEpidermis(grid, ctx);
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // FIELDS_STRATUM_H_
