#ifndef ROS_PDE_H_
#define ROS_PDE_H_

#include "core/pde.h"
#include "infra/sim_param.h"

namespace bdm {
namespace skibidy {

// Reactive oxygen species field: produced by neutrophil NADPH oxidase
// respiratory burst, M1 macrophage oxidative burst, and mitochondrial
// electron transport chain leak. Cleared enzymatically (SOD, catalase,
// GPx) and by vascular perfusion. Mediates oxidative tissue damage that
// connects hyperglycemia to impaired wound healing in diabetic wounds.
struct ROSPDE : public PDE {
  explicit ROSPDE(const SimParam* sp)
      : diffusion_(sp->ros.diffusion),
        decay_(sp->ros.decay) {}

  const char* GetName() const override { return fields::kROS; }
  int GetId() const override { return fields::kROSId; }

  void Init(Simulation* sim) override {
    DefineGrid(sim, diffusion_, decay_);
    // Starts at 0; produced by immune cells and mitochondria after wounding
  }

 private:
  real_t diffusion_;
  real_t decay_;
};

}  // namespace skibidy
}  // namespace bdm

#endif  // ROS_PDE_H_
