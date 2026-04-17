#ifndef ROS_SOURCE_HOOK_H_
#define ROS_SOURCE_HOOK_H_

#include "core/hook_api.h"

namespace bdm {
namespace skibidy {

// Source hook for ROS module.
// Mitochondrial ETC leak produces superoxide constitutively. Rate increases
// under hypoxia (reverse electron transport) and hyperglycemia.
// Tissue antioxidant clearance scales with vascular density.
// Nishikawa et al. 2000 (doi:10.1038/35008121)
struct ROSSourceHook {
  DiffusionGrid* ros_grid = nullptr;
  const SimParam* sp_ = nullptr;
  bool active = false;

  inline void Init(const GridRegistry& reg, SignalBoard& sig) {
    sp_ = reg.Params();
    active = sp_->ros.enabled;
    if (!active) return;
    ros_grid = reg.Get(fields::kROSId);
    if (!ros_grid) { active = false; return; }
  }

  // Called per dermal voxel. Reads O2 and vascular from snapshot.
  inline void ApplyDermal(const VoxelSnapshot& snap, SignalBoard& sig) {
    if (!snap.in_wound || !snap.post_wound) return;

    // Mitochondrial production
    real_t mito_rate = sp_->ros.mitochondrial_rate;
    // Hypoxia amplification: reduced O2 causes ETC backup
    if (snap.o2 < sp_->ros.hypoxia_threshold) {
      real_t hypoxia_frac = (sp_->ros.hypoxia_threshold - snap.o2) /
                             sp_->ros.hypoxia_threshold;
      mito_rate *= (1.0 + (sp_->ros.hypoxia_amplification - 1.0) *
                    hypoxia_frac);
    }
    // Diabetic: hyperglycemia drives mitochondrial superoxide
    if (sp_->diabetic.mode) {
      mito_rate *= sp_->diabetic.mito_ros_factor;
    }
    if (mito_rate > 1e-10) {
      ros_grid->ChangeConcentrationBy(snap.idx, mito_rate);
    }

    // Tissue antioxidant clearance (SOD, catalase, GPx)
    real_t ros_cur = ros_grid->GetConcentration(snap.idx);
    if (ros_cur > 1e-10) {
      real_t tissue_density = std::max(static_cast<real_t>(0.1), snap.vasc);
      real_t antioxidant = sp_->ros.tissue_antioxidant * tissue_density;
      if (sp_->diabetic.mode) {
        antioxidant *= sp_->diabetic.antioxidant_factor;
      }
      real_t clear = std::min(ros_cur, antioxidant * ros_cur);
      if (clear > 1e-10) {
        ros_grid->ChangeConcentrationBy(snap.idx, -clear);
      }
    }

    // Perfusion clears ROS byproducts (vascular washout)
    ros_cur = ros_grid->GetConcentration(snap.idx);
    if (ros_cur > 1e-10 && snap.vasc > 0.1) {
      real_t clear = std::min(ros_cur,
                               sp_->ros.perfusion_clearance * snap.vasc);
      ros_grid->ChangeConcentrationBy(snap.idx, -clear);
    }
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // ROS_SOURCE_HOOK_H_
