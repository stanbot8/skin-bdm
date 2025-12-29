#ifndef FIELDS_KGF_H_
#define FIELDS_KGF_H_

#include "core/pde.h"
#include "infra/sim_param.h"

namespace bdm {
namespace skibidy {

// KGF PDE channel.
// Dermal growth factor: high at basement membrane, decays upward.
// Basal keratinocytes read this to modulate G1->S transition.
// Static prescribed profile; no source term, not disrupted by wound.
struct KgfPDE : public PDE {
  const char* GetName() const override { return fields::kKGF; }
  int GetId() const override { return fields::kKGFId; }

  void Init(Simulation* sim) override {
    auto* sp = sim->GetParam()->Get<SimParam>();
    DefineGrid(sim, sp->kgf_diffusion, sp->kgf_decay);
    ModelInitializer::InitializeSubstance(GetId(),
        [sp](real_t x, real_t y, real_t z) {
          return GridContext::ExpDecay(z, sp->kgf_basal_conc,
                                      sp->kgf_decay_length);
        });
    MarkPrescribed(sim);  // paracrine dermal signal, not free diffusion
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // FIELDS_KGF_H_
