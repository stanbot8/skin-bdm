#ifndef MMP_POST_HOOK_H_
#define MMP_POST_HOOK_H_

#include "core/grid_registry.h"
#include "core/voxel_snapshot.h"
#include "core/signal_board.h"

namespace bdm {
namespace skibidy {

// Post hook for MMP module.
// Matrix metalloproteinase cascade: pro-MMP activation, MMP-TIMP
// neutralization, ECM degradation (collagen, fibronectin, elastin, fibrin),
// matrikine positive feedback, and keratinocyte MMP/TIMP production.
// Visse & Nagase 2003 (doi:10.1161/01.res.0000070112.80711.3d)
struct MMPPostHook {
  DiffusionGrid* mmp_grid = nullptr;
  DiffusionGrid* prommp_grid = nullptr;
  DiffusionGrid* timp_grid = nullptr;
  DiffusionGrid* col_grid = nullptr;
  DiffusionGrid* fn_grid = nullptr;
  DiffusionGrid* elastin_grid = nullptr;
  DiffusionGrid* fibrin_grid = nullptr;
  DiffusionGrid* age_grid = nullptr;
  DiffusionGrid* temp_grid = nullptr;
  const SimParam* sp_ = nullptr;
  bool active = false;
  bool do_ecm = false;
  bool do_fn = false;
  bool do_elastin = false;
  bool do_fibrin = false;
  bool do_prommp = false;
  bool do_age = false;

  inline void Init(const GridRegistry& reg, SignalBoard& sig) {
    sp_ = reg.Params();
    do_ecm = sp_->mmp.enabled && sp_->fibroblast.enabled;
    do_prommp = sp_->mmp.enabled;
    active = do_ecm || do_prommp;
    if (!active) return;

    mmp_grid = reg.Get(fields::kMMPId);
    if (!mmp_grid) { active = false; do_prommp = false; return; }

    if (sp_->fibroblast.enabled)
      col_grid = reg.Get(fields::kCollagenId);

    if (do_prommp)
      prommp_grid = reg.Get(fields::kProMMPId);
    timp_grid = reg.Get(fields::kTIMPId);

    do_fn = sp_->fibronectin.enabled;
    if (do_fn)
      fn_grid = reg.Get(fields::kFibronectinId);

    do_elastin = sp_->elastin.enabled && sp_->mmp.enabled;
    if (do_elastin)
      elastin_grid = reg.Get(fields::kElastinId);

    do_fibrin = sp_->hemostasis.enabled;
    if (do_fibrin)
      fibrin_grid = reg.Get(fields::kFibrinId);

    do_age = sp_->glucose_mod.enabled && sp_->diabetic.mode;
    if (do_age)
      age_grid = reg.Get(fields::kAGEId);

    if (sp_->temperature.enabled)
      temp_grid = reg.Get(fields::kTemperatureId);
  }

  // Called per epidermal wound voxel.
  inline void ApplyEpiWound(const VoxelSnapshot& snap, SignalBoard& sig) {
    // -- MMP-mediated ECM degradation --
    real_t total_degraded = 0;
    if (do_ecm) {
      real_t mmp_val = mmp_grid->GetConcentration(snap.idx);
      if (mmp_val > 1e-10) {
        real_t mmp_ph = 1.0 + sp_->ph.mmp_boost * snap.ph;
        // Temperature Q10: enzyme kinetics scale with temperature
        real_t mmp_temp = 1.0;
        if (temp_grid) {
          real_t temp_val = temp_grid->GetConcentration(snap.idx);
          real_t temp_c = temp_val * 37.0;
          mmp_temp = std::pow(sp_->temperature.q10_mmp,
                              (temp_c - 37.0) / 10.0);
        }
        real_t eff_mmp = mmp_val * mmp_ph * mmp_temp;

        // Collagen degradation (structural)
        // AGE cross-linking: glycated collagen resists MMP cleavage
        // (Mott et al. 1997, doi:10.1038/ki.1997.455)
        if (col_grid) {
          size_t col_si = snap.coarse_si;
          real_t col_val = col_grid->GetConcentration(col_si);
          if (col_val > 1e-10) {
            real_t col_degrade_rate = sp_->mmp.collagen_degradation;
            if (do_age && age_grid) {
              real_t age_val = age_grid->GetConcentration(snap.idx);
              col_degrade_rate *= std::max(real_t{0},
                  1.0 - sp_->age_collagen_crosslink * age_val);
            }
            real_t degrade = std::min(col_val,
                eff_mmp * col_degrade_rate);
            if (degrade > 1e-10) {
              col_grid->ChangeConcentrationBy(col_si,
                  -degrade * snap.coarse_w);
              total_degraded += degrade;
            }
          }
        }
        // Fibronectin degradation (structural)
        if (do_fn && fn_grid) {
          size_t fn_si = snap.coarse_si;
          real_t fn_val = fn_grid->GetConcentration(fn_si);
          if (fn_val > 1e-10) {
            real_t fn_loss = std::min(fn_val,
                eff_mmp * sp_->mmp.fibronectin_degradation);
            if (fn_loss > 1e-10) {
              fn_grid->ChangeConcentrationBy(fn_si,
                  -fn_loss * snap.coarse_w);
              total_degraded += fn_loss;
            }
          }
        }
        // Elastin degradation (structural)
        if (do_elastin && elastin_grid) {
          size_t el_si = snap.coarse_si;
          real_t el_val = elastin_grid->GetConcentration(el_si);
          if (el_val > 1e-10) {
            real_t el_loss = std::min(el_val,
                eff_mmp * sp_->elastin.mmp_degradation);
            if (el_loss > 1e-10) {
              elastin_grid->ChangeConcentrationBy(el_si,
                  -el_loss * snap.coarse_w);
              total_degraded += el_loss;
            }
          }
        }
        // Fibrin degradation (structural)
        if (do_fibrin && fibrin_grid) {
          size_t fb_si = snap.coarse_si;
          real_t fb_val = fibrin_grid->GetConcentration(fb_si);
          if (fb_val > 1e-10) {
            real_t fb_loss = std::min(fb_val,
                eff_mmp * sp_->hemostasis.mmp_degradation);
            if (fb_loss > 1e-10) {
              fibrin_grid->ChangeConcentrationBy(fb_si,
                  -fb_loss * snap.coarse_w);
              total_degraded += fb_loss;
            }
          }
        }

        // Matrikine positive feedback: ECM degradation fragments
        // stimulate MMP transcription via NF-kB.
        if (sp_->matrikine_mmp_boost > 0 && total_degraded > 1e-10) {
          mmp_grid->ChangeConcentrationBy(
              snap.idx, sp_->matrikine_mmp_boost * total_degraded);
        }
      }
    }

    // Write total_degraded to SignalBoard for matrikine feedback
    sig.total_ecm_degraded += total_degraded;

    // -- Keratinocyte pro-MMP-1 at wound migration front --
    // Pilcher et al. 1997
    if (sp_->mmp.keratinocyte_rate > 0 && prommp_grid) {
      real_t edge = snap.stratum * (1.0 - snap.stratum) * 4.0;  // peaks at sv=0.5
      if (edge > 1e-10) {
        prommp_grid->ChangeConcentrationBy(snap.idx,
            sp_->mmp.keratinocyte_rate * edge);
      }
    }

    // -- Keratinocyte TIMP-1 at wound migration front --
    // Saarialho-Kere et al. 1995 (doi:10.1172/JCI117883)
    if (timp_grid && sp_->mmp.timp_keratinocyte_rate > 0) {
      real_t edge = snap.stratum * (1.0 - snap.stratum) * 4.0;
      if (edge > 1e-10) {
        real_t rate = sp_->mmp.timp_keratinocyte_rate;
        if (sp_->diabetic.mode) rate *= sp_->diabetic.timp_production_factor;
        timp_grid->ChangeConcentrationBy(snap.idx, rate * edge);
      }
    }

    // -- MMP-TIMP second-order neutralization --
    // Brew et al. 2000 (doi:10.1016/S0167-4838(99)00252-5)
    if (timp_grid && sp_->mmp.timp_inhibition_rate > 0) {
      real_t mmp_val = mmp_grid->GetConcentration(snap.idx);
      real_t timp_val = timp_grid->GetConcentration(snap.idx);
      if (mmp_val > 1e-10 && timp_val > 1e-10) {
        real_t neutralize =
            sp_->mmp.timp_inhibition_rate * mmp_val * timp_val;
        real_t loss = std::min(neutralize,
            std::min(mmp_val, timp_val));
        if (loss > 1e-10) {
          mmp_grid->ChangeConcentrationBy(snap.idx, -loss);
          timp_grid->ChangeConcentrationBy(snap.idx, -loss);
        }
      }
    }

    // -- Pro-MMP zymogen activation cascade --
    if (do_prommp && prommp_grid) {
      real_t prommp_val = prommp_grid->GetConcentration(snap.idx);
      if (prommp_val > 1e-10) {
        real_t mmp_val = mmp_grid->GetConcentration(snap.idx);
        // Basal plasmin activation + autocatalytic (MMP-3 amplifies)
        real_t activation = sp_->mmp.prommp_activation_rate * prommp_val +
            sp_->mmp.prommp_autocatalytic_rate * prommp_val * mmp_val;
        activation = std::min(activation, prommp_val);
        if (activation > 1e-10) {
          prommp_grid->ChangeConcentrationBy(snap.idx, -activation);
          mmp_grid->ChangeConcentrationBy(snap.idx, activation);
        }
      }
    }
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // MMP_POST_HOOK_H_
