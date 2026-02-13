#ifndef FUSED_SOURCE_H_
#define FUSED_SOURCE_H_

#include <cmath>
#include <vector>

#include "tissue/calcium.h"
#include "core/field_names.h"
#include "core/pde.h"
#include "infra/sim_param.h"
#include "core/operation/operation.h"

namespace bdm {
namespace skibidy {

// ---------------------------------------------------------------------------
// FusedWoundSourceOp -- single-pass replacement for 5 separate source-term
// loops when wound_enabled=true.
//
// Fuses: VascularPDE::ApplySource, OxygenPDE::ApplySource (PinDermal),
//        CalciumPDE::ApplySource, WaterPDE::ApplySource (both sub-loops),
//        VEGFSourceOp.
//
// One iteration over the grid replaces 6-7 separate iterations, each of which
// was computing InWound() and Z() redundantly. A precomputed wound mask +
// layer map eliminates all per-voxel geometry checks.
//
// Reads/writes: Vascular, Oxygen, Calcium, Water, VEGF, Hyaluronan,
//               Dermis, Inflammation, Stratum (9 fields total)
// Tight coupling is intentional: ~7x speedup from single-pass iteration.
// ---------------------------------------------------------------------------
struct FusedWoundSourceOp : public StandaloneOperationImpl {
  BDM_OP_HEADER(FusedWoundSourceOp);

  void operator()() override {
    auto* sim = Simulation::GetActive();
    auto* sp = sim->GetParam()->Get<SimParam>();
    PerfTimer timer(sp->debug_perf);
    auto* scheduler = sim->GetScheduler();
    uint64_t step = scheduler->GetSimulatedSteps();
    uint64_t wound_step = static_cast<uint64_t>(sp->wound_trigger_step);

    auto* rm = sim->GetResourceManager();

    // --- Cache grid pointers ---
    auto* vasc_grid = rm->GetDiffusionGrid(fields::kVascular);
    auto* o2_grid = rm->GetDiffusionGrid(fields::kOxygen);
    auto* stratum_grid = rm->GetDiffusionGrid(fields::kStratum);
    auto* ca_grid = rm->GetDiffusionGrid(fields::kCalcium);
    auto* water_grid = rm->GetDiffusionGrid(fields::kWater);

    DiffusionGrid* vegf_grid = nullptr;
    if (sp->angiogenesis_enabled) {
      vegf_grid = rm->GetDiffusionGrid(fields::kVEGF);
    }

    DiffusionGrid* ha_grid = nullptr;
    bool do_ha_water = sp->hyaluronan_enabled;
    real_t ha_water_factor = sp->hyaluronan_water_retention_factor;
    if (do_ha_water) {
      ha_grid = rm->GetDiffusionGrid(fields::kHyaluronan);
    }

    DiffusionGrid* ph_grid = rm->GetDiffusionGrid(fields::kPH);
    real_t ph_rate = sp->ph_recovery_rate;

    DiffusionGrid* dermis_grid = nullptr;
    DiffusionGrid* col_grid_dermis = nullptr;
    DiffusionGrid* mmp_grid_dermis = nullptr;
    bool do_dermis = sp->dermis_enabled;
    if (do_dermis) {
      dermis_grid = rm->GetDiffusionGrid(fields::kDermis);
      if (sp->fibroblast_enabled)
        col_grid_dermis = rm->GetDiffusionGrid(fields::kCollagen);
      if (sp->mmp_enabled)
        mmp_grid_dermis = rm->GetDiffusionGrid(fields::kMMP);
    }

    GridContext ctx(vasc_grid, sp);

    // --- Lazy-init wound mask + layer map ---
    if (mask_.mask.empty()) {
      real_t z_max = sp->volume_z_cornified + ctx.box_len;
      mask_ = GridContext::ComputeWoundMask(ctx, z_max, sp);
    }

    // Multi-resolution: detect if structural grids use a coarser resolution.
    // When active, fine-loop reads from structural grids use CoarseIndex;
    // writes to structural grids are scaled by coarse_w to compensate for
    // multiple fine voxels mapping to one coarse voxel.
    bool coarse = sp->grid_resolution_structural > 0 &&
                  sp->grid_resolution_structural != sp->grid_resolution;
    real_t coarse_w = 1.0;
    if (coarse) {
      real_t rf = static_cast<real_t>(sp->grid_resolution);
      real_t rc = static_cast<real_t>(sp->grid_resolution_structural);
      coarse_w = (rc * rc * rc) / (rf * rf * rf);
    }
    // Helper: resolve fine index to structural grid index
    auto sidx = [&](DiffusionGrid* g, size_t fine) -> size_t {
      return coarse ? GridContext::CoarseIndex(fine, ctx, g) : fine;
    };

    // --- Pre-compute constants ---
    real_t o2_basal = sp->oxygen_basal_conc;
    real_t water_basal = sp->water_basal_conc;
    real_t water_recovery = sp->water_recovery_rate;
    real_t water_loss = sp->water_surface_loss_rate;
    real_t water_decay_len = sp->water_decay_length;
    real_t ca_rate = sp->calcium_recovery_rate;

    bool post_wound = (step > wound_step);

    // Vascular recovery eligibility
    bool vasc_eligible = false;
    if (post_wound) {
      uint64_t wound_age = step - wound_step;
      vasc_eligible =
          (wound_age >= static_cast<uint64_t>(sp->perfusion_angio_delay));
    }
    real_t vasc_basal = sp->perfusion_basal;
    real_t vasc_rate = sp->perfusion_angio_rate;

    // Per-layer perfusion targets and angiogenesis rate multipliers
    real_t z_papillary = sp->dermal_z_papillary;
    real_t z_reticular = sp->dermal_z_reticular;
    real_t papillary_fraction = sp->perfusion_papillary_fraction;
    real_t reticular_fraction = sp->perfusion_reticular_fraction;
    real_t hypodermis_fraction = sp->perfusion_hypodermis_fraction;
    real_t angio_papillary = sp->angio_papillary_factor;
    real_t angio_reticular = sp->angio_reticular_factor;
    real_t angio_hypodermis = sp->angio_hypodermis_factor;

    // Wound-derived inflammation source (DAMPs / keratinocyte IL-1alpha)
    DiffusionGrid* infl_grid = nullptr;
    real_t wound_infl_rate = sp->wound_inflammation_source_rate;
    bool do_wound_infl = (post_wound && wound_infl_rate > 0);
    if (do_wound_infl) {
      // Temporal taper: DAMPs cleared by macrophage phagocytosis
      real_t taper = sp->wound_inflammation_source_taper;
      if (taper > 0) {
        uint64_t wound_age = step - wound_step;
        wound_infl_rate *= std::exp(-taper * static_cast<real_t>(wound_age));
        if (wound_infl_rate < 1e-10) do_wound_infl = false;
      }
    }
    if (do_wound_infl) {
      infl_grid = sp->split_inflammation_enabled
                      ? rm->GetDiffusionGrid(fields::kProInflammatory)
                      : rm->GetDiffusionGrid(fields::kInflammation);
      if (!infl_grid) do_wound_infl = false;
    }

    // VEGF production constants
    bool do_vegf = (sp->angiogenesis_enabled && post_wound && vegf_grid);
    real_t vegf_threshold = sp->vegf_hypoxia_threshold;
    real_t vegf_prod_rate = sp->vegf_production_rate;
    if (sp->diabetic_mode) vegf_prod_rate *= sp->diabetic_vegf_factor;
    // PHD2/3 negative feedback: prolyl hydroxylases desensitize HIF-1alpha
    // over time, reducing VEGF transcription even under sustained hypoxia
    // (Berra et al. 2003, doi:10.1093/emboj/cdg291)
    if (do_vegf && sp->vegf_production_taper > 0) {
      uint64_t wound_age = step - wound_step;
      vegf_prod_rate *= std::exp(-sp->vegf_production_taper *
                                  static_cast<real_t>(wound_age));
      if (vegf_prod_rate < 1e-10) do_vegf = false;
    }

    // Dermis recovery constants
    real_t dermis_col_threshold = do_dermis ? sp->dermis_collagen_threshold : 0;
    real_t dermis_col_rate = do_dermis ? sp->dermis_collagen_recovery_rate : 0;
    real_t dermis_mmp_rate = do_dermis ? sp->dermis_mmp_degradation : 0;
    real_t dermis_pap_target = sp->dermis_papillary_density;
    real_t dermis_ret_target = sp->dermis_reticular_density;
    real_t dermis_hyp_target = sp->dermis_hypodermis_density;
    real_t dermis_pap_factor = sp->dermis_papillary_rate_factor;
    real_t dermis_ret_factor = sp->dermis_reticular_rate_factor;
    real_t dermis_hyp_factor = sp->dermis_hypodermis_rate_factor;

    // --- Compact loops: iterate only relevant voxels ---
    // Loop 1: ALL dermal voxels (O2 pin, water pin, optional wound logic)
    for (size_t idx : mask_.dermal_all) {
      bool in_wound = mask_.mask[idx];
      real_t vasc_val = vasc_grid->GetConcentration(idx);

      // -- Vascular recovery (wound dermal only) --
      // Per-layer target and rate based on dermal sub-layer
      real_t layer_fraction, layer_angio;
      real_t z = ctx.Z(idx);
      if (z >= z_papillary) {
        layer_fraction = papillary_fraction;
        layer_angio = angio_papillary;
      } else if (z >= z_reticular) {
        layer_fraction = reticular_fraction;
        layer_angio = angio_reticular;
      } else {
        layer_fraction = hypodermis_fraction;
        layer_angio = angio_hypodermis;
      }
      real_t local_vasc_target = vasc_basal * layer_fraction;

      if (in_wound && vasc_eligible && vasc_val < local_vasc_target - 1e-10) {
        real_t x = ctx.X(idx), y = ctx.Y(idx);
        Real3 above = {x, y, 1.0};
        real_t stratum_val = stratum_grid->GetValue(above);
        real_t demand = std::max(static_cast<real_t>(0.2),
                                 GridContext::StratumGate(stratum_val));

        real_t eff_rate = vasc_rate * layer_angio;
        if (vegf_grid) {
          real_t local_vegf = vegf_grid->GetConcentration(idx);
          eff_rate = sp->angio_vegf_rate * local_vegf * layer_angio;
        }

        real_t delta = eff_rate * demand;
        if (vasc_val + delta > local_vasc_target)
          delta = local_vasc_target - vasc_val;
        if (delta > 1e-10) {
          vasc_grid->ChangeConcentrationBy(idx, delta);
          vasc_val += delta;

          // VEGF consumption (proportional to vascular recovery)
          if (vegf_grid) {
            real_t vegf_cur = vegf_grid->GetConcentration(idx);
            real_t consume =
                std::min(vegf_cur, sp->vegf_consumption_rate * delta);
            if (consume > 1e-10) {
              vegf_grid->ChangeConcentrationBy(idx, -consume);
            }
          }
        }
      }

      // -- VEGF receptor clearance (ALL dermal voxels with vasculature) --
      // VEGFR-2 on endothelial cells internalizes VEGF proportional to
      // vascular density (Ferrara 2003). Creates negative feedback: more
      // vessels = more receptors = faster VEGF clearance.
      if (vegf_grid && vasc_val > 0.1) {
        real_t vegf_cur = vegf_grid->GetConcentration(idx);
        if (vegf_cur > 1e-10) {
          real_t clearance = vegf_cur * vasc_val *
                             sp->vegf_receptor_clearance;
          vegf_grid->ChangeConcentrationBy(idx,
              -std::min(vegf_cur, clearance));
        }
      }

      // -- O2 pin (ALL dermal voxels) --
      // Bohr effect: alkaline pH reduces O2 release from hemoglobin
      real_t ph_alkalinity = ph_grid ?
          ph_grid->GetConcentration(sidx(ph_grid, idx)) : 0;
      real_t o2_target = o2_basal * vasc_val *
                         (1.0 - sp->ph_bohr_factor * ph_alkalinity);
      real_t o2_delta = o2_target - o2_grid->GetConcentration(idx);
      if (std::abs(o2_delta) > 1e-10) {
        o2_grid->ChangeConcentrationBy(idx, o2_delta);
      }

      // -- Water dermal --
      real_t water_target = water_basal * vasc_val;
      if (do_ha_water && ha_grid) {
        real_t local_ha = ha_grid->GetConcentration(sidx(ha_grid, idx));
        water_target *= (1.0 + local_ha * ha_water_factor);
      }
      real_t water_current = water_grid->GetConcentration(idx);
      if (in_wound) {
        if (water_current < water_target) {
          real_t gain = std::min(water_target - water_current, water_recovery);
          water_grid->ChangeConcentrationBy(idx, gain);
        }
      } else {
        real_t delta = water_target - water_current;
        if (std::abs(delta) > 1e-10) {
          water_grid->ChangeConcentrationBy(idx, delta);
        }
      }

      // -- VEGF source (wound dermal hypoxia) --
      if (in_wound && do_vegf) {
        real_t local_o2 = o2_grid->GetConcentration(idx);
        if (local_o2 < vegf_threshold) {
          real_t vegf_prod = vegf_prod_rate *
                             (vegf_threshold - local_o2) / vegf_threshold;
          if (vegf_prod > 1e-10) {
            vegf_grid->ChangeConcentrationBy(idx, vegf_prod);
          }
        }
      }

      // -- Dermis tissue recovery / MMP degradation (wound dermal only) --
      if (do_dermis && in_wound && post_wound) {
        size_t d_idx = sidx(dermis_grid, idx);
        real_t dermis_cur = dermis_grid->GetConcentration(d_idx);

        // Per-layer target and rate (reuses z already computed above)
        real_t d_target, d_factor;
        if (z >= z_papillary) {
          d_target = dermis_pap_target; d_factor = dermis_pap_factor;
        } else if (z >= z_reticular) {
          d_target = dermis_ret_target; d_factor = dermis_ret_factor;
        } else {
          d_target = dermis_hyp_target; d_factor = dermis_hyp_factor;
        }

        // Collagen-driven recovery (collagen is structural)
        if (dermis_cur < d_target && col_grid_dermis) {
          real_t col_val = col_grid_dermis->GetConcentration(
              sidx(col_grid_dermis, idx));
          if (col_val > dermis_col_threshold) {
            real_t gain = dermis_col_rate * d_factor *
                          (col_val - dermis_col_threshold);
            if (dermis_cur + gain > d_target) gain = d_target - dermis_cur;
            if (gain > 1e-10) {
              dermis_grid->ChangeConcentrationBy(d_idx, gain * coarse_w);
              dermis_cur += gain * coarse_w;
            }
          }
        }

        // MMP degradation (MMP is fine grid, dermis is structural)
        if (mmp_grid_dermis && dermis_cur > 0) {
          real_t mmp_val = mmp_grid_dermis->GetConcentration(idx);
          if (mmp_val > 1e-10) {
            real_t loss = std::min(dermis_cur, mmp_val * dermis_mmp_rate);
            if (loss > 1e-10) {
              dermis_grid->ChangeConcentrationBy(d_idx, -loss * coarse_w);
            }
          }
        }

        // Dermal DAMP source: inflammation from damaged tissue (fine grid)
        if (do_wound_infl && dermis_cur < d_target) {
          real_t damage_frac = 1.0 - dermis_cur / d_target;
          real_t infl_gain = wound_infl_rate * damage_frac;
          if (infl_gain > 1e-10) {
            infl_grid->ChangeConcentrationBy(idx, infl_gain);
          }
        }
      }
    }

    // Loop 2: epidermal wound voxels only
    for (size_t idx : mask_.epi_wound) {
      real_t stratum_val = stratum_grid->GetConcentration(idx);

      // -- VEGF source (epidermal wound hypoxia) --
      if (do_vegf) {
        real_t local_o2 = o2_grid->GetConcentration(idx);
        if (local_o2 < vegf_threshold) {
          real_t vegf_prod = vegf_prod_rate *
                             (vegf_threshold - local_o2) / vegf_threshold;
          if (vegf_prod > 1e-10) {
            vegf_grid->ChangeConcentrationBy(idx, vegf_prod);
          }
        }
      }

      // -- Surface DAMPs (keratinocyte IL-1a from open epidermis) --
      if (do_wound_infl && stratum_val < 1.0) {
        real_t infl_gain = wound_infl_rate * 0.25 * (1.0 - stratum_val);
        if (infl_gain > 1e-10) {
          infl_grid->ChangeConcentrationBy(idx, infl_gain);
        }
      }

      // -- Calcium restoration (proportional to re-stratification) --
      // Tight junctions reform progressively as the epithelium stratifies;
      // calcium gradient recovery scales with tissue maturity.
      if (stratum_val >= 0.5) {
        real_t barrier = GridContext::StratumGate(stratum_val);
        real_t z = ctx.Z(idx);
        real_t ca_target = CalciumPDE::Profile(sp, z);
        real_t ca_current = ca_grid->GetConcentration(idx);
        real_t ca_delta = (ca_target - ca_current) * ca_rate * barrier;
        if (std::abs(ca_delta) > 1e-10) {
          ca_grid->ChangeConcentrationBy(idx, ca_delta);
        }
      }

      // -- Water epidermal: evaporation + serum recovery --
      real_t barrier = GridContext::StratumGate(stratum_val);
      real_t evap = water_loss * (1.0 - barrier);
      real_t water_current = water_grid->GetConcentration(idx);
      if (evap > 0 && water_current > 0) {
        real_t loss = std::min(water_current, evap);
        water_grid->ChangeConcentrationBy(idx, -loss);
      }
      real_t z = ctx.Z(idx);
      real_t water_tgt =
          GridContext::ExpDecay(z, water_basal, water_decay_len);
      water_current = water_grid->GetConcentration(idx);
      if (water_current < water_tgt) {
        real_t gain = std::min(water_tgt - water_current, water_recovery);
        water_grid->ChangeConcentrationBy(idx, gain);
      }

      // -- pH recovery (wound acidification) --
      // Tissue metabolism produces lactic acid, lowering wound alkalinity.
      // Rate scales with perfusion (O2 delivery drives aerobic/anaerobic balance)
      // and tissue coverage (intact tissue buffers more effectively).
      if (ph_grid) {
        size_t ph_idx = sidx(ph_grid, idx);
        real_t ph_val = ph_grid->GetConcentration(ph_idx);
        if (ph_val > 1e-6) {
          Real3 below = {ctx.X(idx), ctx.Y(idx), -1.0};
          real_t local_perf = vasc_grid->GetValue(below);
          real_t eff_rate = ph_rate * std::max(static_cast<real_t>(0.2),
                                               local_perf) * barrier;
          if (sp->diabetic_mode) {
            eff_rate *= sp->diabetic_ph_recovery_factor;
          }
          real_t drop = std::min(ph_val, eff_rate);
          ph_grid->ChangeConcentrationBy(ph_idx, -drop * coarse_w);
        }
      }
    }

    // --- Manual decay for prescribed ECM fields (D=0, decay>0) ---
    // These fields skip the FTCS solver (SetTimeStep(1e30) in skibidy.h).
    // Apply simple exponential decay: c -= c * mu * dt  (first-order Euler).
    // Structural grids at coarser resolution use CoarseIndex + write scaling.
    real_t dt = sim->GetParam()->simulation_time_step;
    if (sp->fibronectin_enabled && sp->fibronectin_decay > 0) {
      auto* fn_grid = rm->GetDiffusionGrid(fields::kFibronectin);
      real_t fn_mu_dt = sp->fibronectin_decay * dt;
      for (size_t idx : mask_.epi_wound) {
        size_t si = sidx(fn_grid, idx);
        real_t v = fn_grid->GetConcentration(si);
        if (v > 1e-10) fn_grid->ChangeConcentrationBy(si, -v * fn_mu_dt * coarse_w);
      }
      for (size_t idx : mask_.dermal_wound) {
        size_t si = sidx(fn_grid, idx);
        real_t v = fn_grid->GetConcentration(si);
        if (v > 1e-10) fn_grid->ChangeConcentrationBy(si, -v * fn_mu_dt * coarse_w);
      }
    }
    if (sp->elastin_enabled && sp->elastin_decay > 0) {
      auto* el_grid = rm->GetDiffusionGrid(fields::kElastin);
      real_t el_mu_dt = sp->elastin_decay * dt;
      for (size_t idx : mask_.dermal_all) {
        size_t si = sidx(el_grid, idx);
        real_t v = el_grid->GetConcentration(si);
        if (v > 1e-10) el_grid->ChangeConcentrationBy(si, -v * el_mu_dt * coarse_w);
      }
    }
    // Fibrin: D=0, manual decay + coupling to TGF-beta and fibronectin
    if (sp->hemostasis_enabled) {
      auto* fb_grid = rm->GetDiffusionGrid(fields::kFibrin);
      // Manual decay
      if (sp->hemostasis_decay > 0) {
        real_t fb_mu_dt = sp->hemostasis_decay * dt;
        for (size_t idx : mask_.epi_wound) {
          size_t si = sidx(fb_grid, idx);
          real_t v = fb_grid->GetConcentration(si);
          if (v > 1e-10) fb_grid->ChangeConcentrationBy(si, -v * fb_mu_dt * coarse_w);
        }
      }
      // Fibrin -> TGF-beta and fibronectin coupling
      if (post_wound) {
        DiffusionGrid* tgfb_target = sp->fibroblast_enabled ?
            rm->GetDiffusionGrid(fields::kTGFBeta) : nullptr;
        DiffusionGrid* fn_target = sp->fibronectin_enabled ?
            rm->GetDiffusionGrid(fields::kFibronectin) : nullptr;
        real_t tgfb_rate = sp->hemostasis_tgfb_coupling;
        real_t fn_rate = sp->hemostasis_fibronectin_coupling;
        if (tgfb_target || fn_target) {
          for (size_t idx : mask_.epi_wound) {
            size_t fb_si = sidx(fb_grid, idx);
            real_t fb = fb_grid->GetConcentration(fb_si);
            if (fb < 1e-10) continue;
            // TGF-beta is fine grid: write at fine idx, no correction
            if (tgfb_target && tgfb_rate > 0) {
              tgfb_target->ChangeConcentrationBy(idx, fb * tgfb_rate);
            }
            // Fibronectin is structural: write at coarse idx, scaled
            if (fn_target && fn_rate > 0) {
              fn_target->ChangeConcentrationBy(
                  sidx(fn_target, idx), fb * fn_rate * coarse_w);
            }
          }
        }
      }
    }
    timer.Print("fused_source");
  }

 private:
  GridContext::WoundMaskData mask_;
};

}  // namespace skibidy
}  // namespace bdm

#endif  // FUSED_SOURCE_H_
