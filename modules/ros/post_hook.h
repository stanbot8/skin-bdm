#ifndef ROS_POST_HOOK_H_
#define ROS_POST_HOOK_H_

#include "core/grid_registry.h"
#include "core/voxel_snapshot.h"
#include "core/signal_board.h"

namespace bdm {
namespace skibidy {

// Post hook for ROS module.
// ROS oxidizes cellular and ECM targets: DNA (senescence), pro-MMP cysteine
// switch (activation), collagen cross-links (damage), NF-kB (inflammation).
// Schafer & Werner 2008 (doi:10.1016/j.phrs.2008.06.004)
struct ROSPostHook {
  DiffusionGrid* ros_grid = nullptr;
  DiffusionGrid* infl_grid = nullptr;
  DiffusionGrid* prommp_grid = nullptr;
  DiffusionGrid* mmp_grid = nullptr;
  DiffusionGrid* col_grid = nullptr;
  DiffusionGrid* sen_grid = nullptr;
  const SimParam* sp_ = nullptr;
  bool active = false;
  bool do_senescence = false;
  bool do_prommp = false;

  inline void Init(const GridRegistry& reg, SignalBoard& sig) {
    sp_ = reg.Params();
    active = sp_->ros.enabled;
    if (!active) return;
    ros_grid = reg.Get(fields::kROSId);
    if (!ros_grid) { active = false; return; }
    infl_grid = reg.InflammationGrid();
    do_prommp = sp_->mmp.enabled;
    if (do_prommp) {
      prommp_grid = reg.Get(fields::kProMMPId);
      mmp_grid = reg.Get(fields::kMMPId);
    }
    if (sp_->fibroblast.enabled)
      col_grid = reg.Get(fields::kCollagenId);
    do_senescence = sp_->senescence.enabled;
    if (do_senescence)
      sen_grid = reg.Get(fields::kSenescenceId);
  }

  // Called per epidermal wound voxel.
  inline void ApplyEpiWound(const VoxelSnapshot& snap, SignalBoard& sig) {
    real_t ros_val = ros_grid->GetConcentration(snap.idx);
    if (ros_val <= 1e-10) return;

    // DNA damage -> senescence accumulation
    if (do_senescence && sen_grid && sp_->ros.senescence_rate > 0) {
      size_t sen_si = snap.coarse_si;
      real_t sen_gain = sp_->ros.senescence_rate * ros_val;
      sen_grid->ChangeConcentrationBy(sen_si, sen_gain * snap.coarse_w);
    }

    // Oxidative pro-MMP activation (cysteine switch oxidation)
    // Rajagopalan et al. 1996 (doi:10.1074/jbc.271.24.14030)
    if (do_prommp && prommp_grid && mmp_grid &&
        sp_->ros.mmp_activation > 0) {
      real_t prommp_val = prommp_grid->GetConcentration(snap.idx);
      if (prommp_val > 1e-10) {
        real_t activation = sp_->ros.mmp_activation * ros_val * prommp_val;
        activation = std::min(activation, prommp_val);
        if (activation > 1e-10) {
          prommp_grid->ChangeConcentrationBy(snap.idx, -activation);
          mmp_grid->ChangeConcentrationBy(snap.idx, activation);
        }
      }
    }

    // NF-kB activation: ROS amplifies pro-inflammatory signaling
    if (infl_grid && sp_->ros.inflammation_amplification > 0) {
      real_t gate = std::max(real_t{0},
                             static_cast<real_t>(1) - snap.stratum);
      if (gate > 1e-10) {
        infl_grid->ChangeConcentrationBy(snap.idx,
            sp_->ros.inflammation_amplification * ros_val * gate);
      }
    }

    // Oxidative collagen damage (cross-link disruption)
    if (col_grid && sp_->ros.collagen_damage > 0) {
      size_t col_si = snap.coarse_si;
      real_t col_val = col_grid->GetConcentration(col_si);
      if (col_val > 1e-10) {
        real_t damage = std::min(col_val,
            sp_->ros.collagen_damage * ros_val * col_val);
        if (damage > 1e-10) {
          col_grid->ChangeConcentrationBy(col_si, -damage * snap.coarse_w);
        }
      }
    }
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // ROS_POST_HOOK_H_
