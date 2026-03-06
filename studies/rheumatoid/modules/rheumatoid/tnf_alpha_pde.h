#ifndef TNF_ALPHA_PDE_H_
#define TNF_ALPHA_PDE_H_

#include "core/pde.h"
#include "infra/sim_param.h"

namespace bdm {
namespace skibidy {

// TNF-alpha cytokine field: primary driver of rheumatoid arthritis pathology.
// Produced by synovial macrophages (M1), neutrophils, and fibroblast-like
// synoviocytes (FLS). Drives NF-kB inflammatory cascade, MMP transcription
// for cartilage/bone erosion, and pannus angiogenesis via VEGF induction.
// Cleared by receptor endocytosis and anti-TNF biologics (treatment).
// Feldmann & Maini 2003 (doi:10.1038/nm939)
struct TNFAlphaPDE : public PDE {
  explicit TNFAlphaPDE(const SimParam* sp)
      : diffusion_(sp->tnf_alpha_diffusion),
        decay_(sp->tnf_alpha_decay) {}

  const char* GetName() const override { return fields::kTNFAlpha; }
  int GetId() const override { return fields::kTNFAlphaId; }

  void Init(Simulation* sim) override {
    DefineGrid(sim, diffusion_, decay_);
    // Starts at 0; produced by immune cells and autoimmune source after wounding
  }

 private:
  real_t diffusion_;
  real_t decay_;
};

}  // namespace skibidy
}  // namespace bdm

#endif  // TNF_ALPHA_PDE_H_
