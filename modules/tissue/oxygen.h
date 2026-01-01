#ifndef FIELDS_OXYGEN_H_
#define FIELDS_OXYGEN_H_

#include "core/pde.h"
#include "core/composite_field.h"
#include "infra/sim_param.h"

namespace bdm {
namespace skibidy {

// Oxygen PDE channel.
// Dermal vasculature supplies O2 at basement membrane. O2 diffuses through
// tissue and is consumed by cells. Source term reads the Vascular field for
// perfusion state, coupling O2 supply to explicit vessel integrity.
struct OxygenPDE : public PDE {
  const char* GetName() const override { return fields::kOxygen; }
  int GetId() const override { return fields::kOxygenId; }

  void Init(Simulation* sim) override {
    auto* sp = sim->GetParam()->Get<SimParam>();
    DefineGrid(sim, sp->oxygen_diffusion, sp->oxygen_decay);
    ModelInitializer::InitializeSubstance(GetId(),
        [sp](real_t x, real_t y, real_t z) {
          return GridContext::ExpDecay(z, sp->oxygen_basal_conc,
                                      sp->oxygen_decay_length);
        });
  }

  void ApplyWound(Simulation* sim, real_t cx, real_t cy, real_t r) override {
    ZeroInWound(sim);
  }

  // Dermal O2 source: pins dermal voxels to target proportional to
  // local vascular perfusion. All damage/recovery logic lives in VascularPDE.
  void ApplySource(Simulation* sim, const CompositeField& fields) override {
    auto* sp = sim->GetParam()->Get<SimParam>();
    auto* o2_grid = Grid(sim);
    auto* vasc_grid = fields.Grid(fields::kVascular, sim);
    GridContext ctx(o2_grid, sp);
    GridContext::PinDermal(o2_grid, vasc_grid, ctx, sp->oxygen_basal_conc);
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // FIELDS_OXYGEN_H_
