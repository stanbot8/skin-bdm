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
    // Allow wound void value of -1 (BDM default lower_threshold_ is 0).
    Grid(sim)->SetLowerThreshold(-2.0);
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
    // Set wound void to -1 (not 0) so it is distinguishable from kBasal=0.
    // Agents writing kBasal=0 will overwrite -1 as tissue arrives.
    for (size_t idx = 0; idx < ctx.n; idx++) {
      if (ctx.Z(idx) < 0) continue;
      if (ctx.InWound(ctx.X(idx), ctx.Y(idx))) {
        real_t val = grid->GetConcentration(idx);
        real_t delta = -1.0 - val;
        if (std::abs(delta) > 1e-10) {
          grid->ChangeConcentrationBy(idx, delta);
        }
      }
    }
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // FIELDS_STRATUM_H_
