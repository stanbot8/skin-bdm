#ifndef ANGIOGENESIS_SOURCE_HOOK_H_
#define ANGIOGENESIS_SOURCE_HOOK_H_

#include "core/hook_api.h"

namespace bdm {
namespace skibidy {

// Source hook for angiogenesis module.
// Vascular recovery (VEGF-driven angiogenesis with per-layer geometry),
// VEGF production from hypoxic tissue, VEGF receptor clearance, VEGF
// consumption by sprouting endothelial cells, ROS impairment of VEGFR2,
// and NO vasodilation.
// Ferrara 2003 (doi:10.1038/nm0603-669)
// Berra et al. 2003 (doi:10.1093/emboj/cdg392)
struct AngioSourceHook {
  DiffusionGrid* vasc_grid = nullptr;
  DiffusionGrid* vegf_grid = nullptr;
  DiffusionGrid* stratum_grid = nullptr;
  DiffusionGrid* o2_grid = nullptr;
  DiffusionGrid* ros_grid = nullptr;
  DiffusionGrid* no_grid = nullptr;
  const SimParam* sp_ = nullptr;

  bool active = false;
  bool do_no = false;

  inline void Init(const GridRegistry& reg, SignalBoard& sig) {
    sp_ = reg.Params();
    active = true;  // always active within fused source op

    vasc_grid = reg.Get(fields::kVascularId);
    stratum_grid = reg.Get(fields::kStratumId);
    o2_grid = reg.Get(fields::kOxygenId);

    if (sp_->angiogenesis.enabled)
      vegf_grid = reg.Get(fields::kVEGFId);
    if (sp_->ros.enabled)
      ros_grid = reg.Get(fields::kROSId);
    do_no = sp_->nitric_oxide.enabled;
    if (do_no) {
      no_grid = reg.Get(fields::kNitricOxideId);
      if (!no_grid) do_no = false;
    }

    sig.vasc_eligible = false;
    if (reg.PostWound()) {
      sig.vasc_eligible =
          (reg.WoundAge() >= static_cast<uint64_t>(sp_->perfusion.angio_delay));
    }

    sig.do_vegf = sp_->angiogenesis.enabled && reg.PostWound() && vegf_grid;
    sig.vegf_prod_rate = sp_->angiogenesis.vegf_production_rate;
    if (sp_->diabetic.mode) sig.vegf_prod_rate *= sp_->diabetic.vegf_factor;
    if (sig.do_vegf && sp_->angiogenesis.vegf_production_taper > 0) {
      sig.vegf_prod_rate *= std::exp(-sp_->angiogenesis.vegf_production_taper *
                                      static_cast<real_t>(reg.WoundAge()));
      if (sig.vegf_prod_rate < 1e-10) sig.do_vegf = false;
    }
    sig.vegf_threshold = sp_->angiogenesis.vegf_hypoxia_threshold;
  }

  // Hypoxia-driven VEGF production (shared by dermal and epidermal wound).
  inline void ProduceVEGF(size_t idx, real_t local_o2, const SignalBoard& sig) {
    if (!sig.do_vegf) return;
    if (local_o2 < sp_->angiogenesis.vegf_hypoxia_threshold) {
      real_t vegf_prod = sig.vegf_prod_rate *
          (sp_->angiogenesis.vegf_hypoxia_threshold - local_o2) /
          sp_->angiogenesis.vegf_hypoxia_threshold;
      if (vegf_prod > 1e-10) {
        vegf_grid->ChangeConcentrationBy(idx, vegf_prod);
      }
    }
  }

  // Dermal: vascular recovery + VEGF clearance + NO vasodilation + VEGF source
  inline void ApplyDermal(const VoxelSnapshot& snap, SignalBoard& sig) {
    real_t vasc_val = vasc_grid->GetConcentration(snap.idx);

    // Per-layer geometry
    real_t layer_fraction, layer_angio;
    if (snap.z >= sp_->dermal_z_papillary) {
      layer_fraction = sp_->perfusion.papillary_fraction;
      layer_angio = sp_->angio_papillary_factor;
    } else if (snap.z >= sp_->dermal_z_reticular) {
      layer_fraction = sp_->perfusion.reticular_fraction;
      layer_angio = sp_->angio_reticular_factor;
    } else {
      layer_fraction = sp_->perfusion.hypodermis_fraction;
      layer_angio = sp_->angio_hypodermis_factor;
    }
    real_t local_vasc_target = sp_->perfusion.basal * layer_fraction;

    // Vascular recovery (wound dermal only)
    if (snap.in_wound && sig.vasc_eligible && vasc_val < local_vasc_target - 1e-10) {
      Real3 above = {snap.x, snap.y, 1.0};
      real_t stratum_val = stratum_grid->GetValue(above);
      real_t demand = std::max(static_cast<real_t>(0.2),
                               GridContext::StratumGate(stratum_val));

      real_t eff_rate = sp_->perfusion.angio_rate * layer_angio;
      if (vegf_grid) {
        real_t local_vegf = vegf_grid->GetConcentration(snap.idx);
        eff_rate = sp_->angio_vegf_rate * local_vegf * layer_angio;
      }
      // ROS impairs endothelial VEGFR2 signaling
      if (ros_grid) {
        real_t ros_val = ros_grid->GetConcentration(snap.idx);
        if (ros_val > 1e-10) {
          eff_rate *= std::max(real_t{0},
              1.0 - sp_->ros.angiogenesis_impairment * ros_val);
        }
      }

      real_t delta = eff_rate * demand;
      if (vasc_val + delta > local_vasc_target)
        delta = local_vasc_target - vasc_val;
      if (delta > 1e-10) {
        vasc_grid->ChangeConcentrationBy(snap.idx, delta);
        vasc_val += delta;

        // VEGF consumption proportional to vascular recovery
        if (vegf_grid) {
          real_t vegf_cur = vegf_grid->GetConcentration(snap.idx);
          real_t consume =
              std::min(vegf_cur, sp_->angiogenesis.vegf_consumption_rate * delta);
          if (consume > 1e-10) {
            vegf_grid->ChangeConcentrationBy(snap.idx, -consume);
          }
        }
      }
    }

    // VEGF receptor clearance (all dermal voxels with vasculature)
    if (vegf_grid && vasc_val > 0.1) {
      real_t vegf_cur = vegf_grid->GetConcentration(snap.idx);
      if (vegf_cur > 1e-10) {
        real_t clearance = vegf_cur * vasc_val *
                           sp_->angiogenesis.vegf_receptor_clearance;
        vegf_grid->ChangeConcentrationBy(snap.idx,
            -std::min(vegf_cur, clearance));
      }
    }

    // VEGF source (wound dermal hypoxia)
    if (snap.in_wound)
      ProduceVEGF(snap.idx, o2_grid->GetConcentration(snap.idx), sig);

    // NO vasodilation: boost effective perfusion
    if (do_no && snap.in_wound && sig.vasc_eligible) {
      real_t no_val = no_grid->GetConcentration(snap.idx);
      if (no_val > 1e-10 && vasc_val < local_vasc_target) {
        real_t no_boost = sp_->perfusion.angio_rate *
            sp_->no_vasodilation_factor * no_val;
        real_t headroom = local_vasc_target - vasc_val;
        real_t gain = std::min(headroom, no_boost);
        if (gain > 1e-10) {
          vasc_grid->ChangeConcentrationBy(snap.idx, gain);
        }
      }
    }
  }

  // Epidermal wound: VEGF source (hypoxia-driven production)
  inline void ApplyEpiWound(const VoxelSnapshot& snap, SignalBoard& sig) {
    ProduceVEGF(snap.idx, snap.o2, sig);
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // ANGIOGENESIS_SOURCE_HOOK_H_
