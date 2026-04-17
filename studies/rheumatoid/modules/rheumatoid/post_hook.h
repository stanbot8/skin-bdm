#ifndef RHEUMATOID_POST_HOOK_H_
#define RHEUMATOID_POST_HOOK_H_

#include "core/hook_api.h"

namespace bdm {
namespace skibidy {

// Post hook for rheumatoid arthritis module.
// Two dispatch paths:
//   ApplyDermal: cartilage + bone erosion (deep joint tissue, no stratum gate)
//   ApplyEpiWound: TNF/IL-6 signaling cascades (stratum-gated surface effects)
//
// Erosion pathways:
//   Cartilage: MMP proteolysis + TNF osteoclast + IL-6 RANKL
//   Bone: RANKL/osteoclast (IL-6) + direct TNF activation
//   Pannus proximity amplifies both erosion types
// Signaling pathways:
//   TNF -> NF-kB -> inflammation, MMP transcription, VEGF
//   IL-6 -> JAK/STAT3 -> inflammation amplification
// Firestein 2003 (doi:10.1038/nature01661)
// McInnes & Schett 2011 (doi:10.1056/NEJMra1004965)
// Schett & Gravallese 2012 (doi:10.1038/nrrheum.2012.153)
struct RAPostHook {
  DiffusionGrid* tnf_grid = nullptr;
  DiffusionGrid* il6_grid = nullptr;
  DiffusionGrid* cartilage_grid = nullptr;
  DiffusionGrid* bone_grid = nullptr;
  DiffusionGrid* synovial_grid = nullptr;
  DiffusionGrid* vegf_grid = nullptr;
  DiffusionGrid* infl_grid = nullptr;
  DiffusionGrid* prommp_grid = nullptr;
  DiffusionGrid* mmp_grid = nullptr;
  const SimParam* sp_ = nullptr;
  bool active = false;

  inline void Init(const GridRegistry& reg, SignalBoard& sig) {
    sp_ = reg.Params();
    active = sp_->ra.enabled;
    if (!active) return;
    tnf_grid = reg.Get(fields::kTNFAlphaId);
    il6_grid = reg.Get(fields::kIL6Id);
    cartilage_grid = reg.Get(fields::kCartilageId);
    if (!tnf_grid || !il6_grid || !cartilage_grid) { active = false; return; }
    bone_grid = reg.Get(fields::kBoneId);
    infl_grid = reg.InflammationGrid();
    synovial_grid = reg.Get(fields::kSynovialFluidId);
    if (sp_->mmp.enabled) {
      prommp_grid = reg.Get(fields::kProMMPId);
      mmp_grid = reg.Get(fields::kMMPId);
    }
    if (sp_->angiogenesis.enabled)
      vegf_grid = reg.Get(fields::kVEGFId);
  }

  // Dermal loop: cartilage and bone erosion (joint tissue is deep, not epidermal).
  // Runs on ALL dermal wound voxels, independent of stratum/re-epithelialization.
  inline void ApplyDermal(const VoxelSnapshot& snap, SignalBoard& sig) {
    if (!snap.post_wound) return;

    real_t tnf_val = tnf_grid->GetConcentration(snap.idx);
    real_t il6_val = il6_grid->GetConcentration(snap.idx);
    if (tnf_val <= 1e-10 && il6_val <= 1e-10) return;

    // Pannus factor: amplifies both cartilage and bone erosion
    real_t pannus_factor = 1.0;
    if (synovial_grid) {
      real_t syn_val = synovial_grid->GetConcentration(snap.coarse_si);
      pannus_factor = 1.0 + sp_->ra.synovial_erosion_boost * syn_val;
    }

    // --- Cartilage degradation (combined TNF + MMP + IL-6/RANKL) ---
    size_t cart_si = snap.coarse_si;
    real_t cart_val = cartilage_grid->GetConcentration(cart_si);

    if (cart_val > 1e-10) {
      real_t mmp_erosion = 0;
      if (mmp_grid) {
        real_t mmp_val = mmp_grid->GetConcentration(snap.idx);
        mmp_erosion = sp_->ra.cartilage_mmp_degradation * mmp_val;
      }
      real_t tnf_erosion = sp_->ra.cartilage_tnf_degradation * tnf_val;
      real_t il6_erosion = sp_->ra.il6_cartilage_boost * il6_val;

      real_t total_loss = std::min(cart_val,
          (mmp_erosion + tnf_erosion + il6_erosion) * cart_val * pannus_factor);
      if (total_loss > 1e-10) {
        cartilage_grid->ChangeConcentrationBy(cart_si,
            -total_loss * snap.coarse_w);
      }
    }

    // --- Bone erosion (RANKL/osteoclast + TNF direct) ---
    // Schett & Gravallese 2012: focal bone erosions at pannus-bone interface.
    if (bone_grid) {
      size_t bone_si = snap.coarse_si;
      real_t bone_val = bone_grid->GetConcentration(bone_si);
      if (bone_val > 1e-10) {
        real_t rankl = sp_->ra.bone_rankl_erosion * il6_val;
        real_t tnf_bone = sp_->ra.bone_tnf_erosion * tnf_val;
        real_t bone_loss = std::min(bone_val,
            (rankl + tnf_bone) * bone_val * pannus_factor);
        if (bone_loss > 1e-10) {
          bone_grid->ChangeConcentrationBy(bone_si,
              -bone_loss * snap.coarse_w);
        }
      }
    }
  }

  // Epi wound loop: TNF/IL-6 signaling effects on surface biology.
  // Stratum-gated: strongest where tissue is disrupted.
  inline void ApplyEpiWound(const VoxelSnapshot& snap, SignalBoard& sig) {
    real_t tnf_val = tnf_grid->GetConcentration(snap.idx);
    real_t il6_val = il6_grid->GetConcentration(snap.idx);
    if (tnf_val <= 1e-10 && il6_val <= 1e-10) return;

    real_t gate = std::max(static_cast<real_t>(0),
                           static_cast<real_t>(1) - snap.stratum);
    if (gate <= 1e-10) return;

    if (tnf_val > 1e-10) {
      if (infl_grid && sp_->ra.tnf_inflammation_coupling > 0) {
        infl_grid->ChangeConcentrationBy(snap.idx,
            sp_->ra.tnf_inflammation_coupling * tnf_val * gate);
      }
      if (prommp_grid && sp_->ra.tnf_mmp_boost > 0) {
        prommp_grid->ChangeConcentrationBy(snap.idx,
            sp_->ra.tnf_mmp_boost * tnf_val);
      }
      if (vegf_grid && sp_->ra.tnf_vegf_boost > 0) {
        vegf_grid->ChangeConcentrationBy(snap.idx,
            sp_->ra.tnf_vegf_boost * tnf_val);
      }
    }

    if (il6_val > 1e-10) {
      if (infl_grid && sp_->ra.il6_inflammation_coupling > 0) {
        infl_grid->ChangeConcentrationBy(snap.idx,
            sp_->ra.il6_inflammation_coupling * il6_val * gate);
      }
    }
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // RHEUMATOID_POST_HOOK_H_
