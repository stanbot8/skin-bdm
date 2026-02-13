#ifndef FUSED_POST_H_
#define FUSED_POST_H_

#include <cmath>
#include <iostream>
#include <vector>

#include "core/field_names.h"
#include "core/pde.h"
#include "infra/sim_param.h"
#include "core/operation/operation.h"

namespace bdm {
namespace skibidy {

// ---------------------------------------------------------------------------
// FusedWoundPostOp -- single-pass replacement for separate post-schedule
// voxel loops when wound_enabled=true.
//
// Fuses: BiofilmGrowthOp, BaselineInflammationOp, ScarAccumulationOp,
//        MMP-collagen degradation, fibronectin serum source.
//
// Only iterates epidermal wound voxels. The ordering within each voxel
// preserves dependencies: biofilm writes inflammation -> baseline writes
// inflammation -> MMP degrades ECM -> fibronectin serum -> scar reads
// inflammation.
//
// Reads/writes: Biofilm, Inflammation, Scar, MMP, Collagen, Fibronectin,
//               Elastin (7 fields total)
// Tight coupling is intentional: ~7x speedup from single-pass iteration.
// ---------------------------------------------------------------------------
struct FusedWoundPostOp : public StandaloneOperationImpl {
  BDM_OP_HEADER(FusedWoundPostOp);

  void operator()() override {
    auto* sim = Simulation::GetActive();
    auto* sp = sim->GetParam()->Get<SimParam>();
    PerfTimer timer(sp->debug_perf);
    auto* scheduler = sim->GetScheduler();
    uint64_t step = scheduler->GetSimulatedSteps();
    uint64_t wound_step = static_cast<uint64_t>(sp->wound_trigger_step);
    if (step <= wound_step) return;
    uint64_t wound_age = step - wound_step;

    bool do_biofilm = sp->biofilm_enabled;
    bool do_baseline = sp->diabetic_mode &&
                       sp->diabetic_baseline_inflammation > 0;
    bool do_scar = sp->scar_proportional_enabled;
    bool do_mmp = sp->mmp_enabled && sp->fibroblast_enabled;
    bool do_fn = sp->fibronectin_enabled;
    bool do_elastin_degrade = sp->elastin_enabled && sp->mmp_enabled;

    if (!do_biofilm && !do_baseline && !do_scar && !do_mmp && !do_fn &&
        !do_elastin_degrade) return;

    auto* rm = sim->GetResourceManager();

    // --- Grid pointers ---
    DiffusionGrid* stratum_grid = rm->GetDiffusionGrid(fields::kStratum);
    GridContext ctx(stratum_grid, sp);

    // Lazy-init wound mask
    if (mask_.mask.empty()) {
      real_t z_max = sp->volume_z_cornified + ctx.box_len;
      mask_ = GridContext::ComputeWoundMask(ctx, z_max, sp);
    }

    // pH grid for MMP and biofilm modulation
    DiffusionGrid* ph_grid = rm->GetDiffusionGrid(fields::kPH);

    DiffusionGrid* biofilm_grid = nullptr;
    if (do_biofilm) {
      biofilm_grid = rm->GetDiffusionGrid(fields::kBiofilm);
    }

    DiffusionGrid* infl_grid = nullptr;
    if (do_biofilm || do_baseline || do_scar) {
      if (sp->split_inflammation_enabled) {
        infl_grid = rm->GetDiffusionGrid(fields::kProInflammatory);
      } else {
        infl_grid = rm->GetDiffusionGrid(fields::kInflammation);
      }
    }

    DiffusionGrid* scar_grid = nullptr;
    if (do_scar) {
      scar_grid = rm->GetDiffusionGrid(fields::kScar);
    }

    DiffusionGrid* mmp_grid = nullptr;
    DiffusionGrid* col_grid = nullptr;
    if (do_mmp) {
      mmp_grid = rm->GetDiffusionGrid(fields::kMMP);
      col_grid = rm->GetDiffusionGrid(fields::kCollagen);
    }

    DiffusionGrid* fn_grid = nullptr;
    if (do_fn) {
      fn_grid = rm->GetDiffusionGrid(fields::kFibronectin);
    }

    DiffusionGrid* elastin_grid = nullptr;
    if (do_elastin_degrade) {
      elastin_grid = rm->GetDiffusionGrid(fields::kElastin);
    }

    bool do_fibrin_degrade = sp->hemostasis_enabled && do_mmp;
    DiffusionGrid* fibrin_grid = nullptr;
    if (do_fibrin_degrade) {
      fibrin_grid = rm->GetDiffusionGrid(fields::kFibrin);
    }

    // Multi-resolution: structural grids may use coarser resolution
    bool coarse = sp->grid_resolution_structural > 0 &&
                  sp->grid_resolution_structural != sp->grid_resolution;
    real_t coarse_w = 1.0;
    if (coarse) {
      real_t rf = static_cast<real_t>(sp->grid_resolution);
      real_t rc = static_cast<real_t>(sp->grid_resolution_structural);
      coarse_w = (rc * rc * rc) / (rf * rf * rf);
    }
    auto sidx = [&](DiffusionGrid* g, size_t fine) -> size_t {
      return coarse ? GridContext::CoarseIndex(fine, ctx, g) : fine;
    };

    // --- Biofilm constants ---
    bool seeding_step = do_biofilm && !seeded_ &&
        (wound_age == static_cast<uint64_t>(sp->biofilm_seed_delay));
    real_t bf_growth = sp->biofilm_growth_rate;
    real_t bf_K = sp->biofilm_carrying_capacity;
    real_t bf_seed = sp->biofilm_seed_amount;
    real_t bf_infl = sp->biofilm_inflammation_rate;

    // --- Other constants ---
    real_t baseline_rate = sp->diabetic_baseline_inflammation;
    real_t scar_rate = sp->scar_accumulation_rate;

    // MMP constants
    real_t mmp_col_degrade = do_mmp ? sp->mmp_collagen_degradation : 0;
    real_t mmp_fn_degrade = (do_mmp && do_fn) ?
        sp->mmp_fibronectin_degradation : 0;
    real_t mmp_el_degrade = do_elastin_degrade ?
        sp->elastin_mmp_degradation : 0;
    real_t mmp_fb_degrade = do_fibrin_degrade ?
        sp->hemostasis_mmp_degradation : 0;

    // Fibronectin serum source (plasma fibronectin seeps into open wound)
    real_t fn_serum = do_fn ? sp->fibronectin_serum_rate : 0;

    // Keratinocyte MMP-1 source (wound-edge migration front)
    real_t mmp_kerat = do_mmp ? sp->mmp_keratinocyte_rate : 0;

    // --- Compact loop: epidermal wound voxels only ---
    for (size_t idx : mask_.epi_wound) {

      // pH alkalinity for MMP and biofilm modulation
      real_t ph_val = ph_grid ?
          ph_grid->GetConcentration(sidx(ph_grid, idx)) : 0;

      // -- Biofilm dynamics (biofilm is structural) --
      if (do_biofilm) {
        real_t stratum_val = stratum_grid->GetConcentration(idx);
        if (stratum_val < 1.0) {
          size_t bf_idx = sidx(biofilm_grid, idx);
          real_t bf = biofilm_grid->GetConcentration(bf_idx);

          if (seeding_step) {
            biofilm_grid->ChangeConcentrationBy(bf_idx, bf_seed * coarse_w);
            bf += bf_seed * coarse_w;
          }

          if (bf > 1e-10) {
            // Alkaline pH promotes bacterial colonization
            real_t eff_growth = bf_growth * (1.0 + sp->ph_biofilm_boost * ph_val);
            real_t delta = eff_growth * bf * (1.0 - bf / bf_K);
            if (delta > 1e-10) {
              biofilm_grid->ChangeConcentrationBy(bf_idx, delta * coarse_w);
            }
            // Inflammation is fine grid: no correction
            real_t infl_delta = bf_infl * bf;
            if (infl_delta > 1e-10) {
              infl_grid->ChangeConcentrationBy(idx, infl_delta);
            }
          }
        }
      }

      // -- Baseline inflammation (diabetic) --
      // Gated by stratum: re-epithelialized tissue seals the AGE/RAGE source.
      // Open wound bed (stratum~0) gets full baseline; healed (stratum>=1) gets none.
      if (do_baseline) {
        real_t sv = stratum_grid->GetConcentration(idx);
        real_t gate = std::max(static_cast<real_t>(0),
                               static_cast<real_t>(1) - sv);
        if (gate > 1e-10) {
          infl_grid->ChangeConcentrationBy(idx, baseline_rate * gate);
        }
      }

      // -- MMP-mediated ECM degradation --
      // MMP (fine grid) degrades structural ECM fields at coarse resolution.
      // Alkaline pH amplifies MMP proteolytic activity (optimal pH ~7.0-8.0).
      if (do_mmp) {
        real_t mmp_val = mmp_grid->GetConcentration(idx);
        if (mmp_val > 1e-10) {
          real_t mmp_ph = 1.0 + sp->ph_mmp_boost * ph_val;
          real_t eff_mmp = mmp_val * mmp_ph;
          // Collagen degradation (structural)
          size_t col_si = sidx(col_grid, idx);
          real_t col_val = col_grid->GetConcentration(col_si);
          if (col_val > 1e-10) {
            real_t degrade = std::min(col_val, eff_mmp * mmp_col_degrade);
            if (degrade > 1e-10) {
              col_grid->ChangeConcentrationBy(col_si, -degrade * coarse_w);
            }
          }
          // Fibronectin degradation (structural)
          if (do_fn && fn_grid) {
            size_t fn_si = sidx(fn_grid, idx);
            real_t fn_val = fn_grid->GetConcentration(fn_si);
            if (fn_val > 1e-10) {
              real_t fn_loss = std::min(fn_val, eff_mmp * mmp_fn_degrade);
              if (fn_loss > 1e-10) {
                fn_grid->ChangeConcentrationBy(fn_si, -fn_loss * coarse_w);
              }
            }
          }
          // Elastin degradation (structural)
          if (do_elastin_degrade && elastin_grid) {
            size_t el_si = sidx(elastin_grid, idx);
            real_t el_val = elastin_grid->GetConcentration(el_si);
            if (el_val > 1e-10) {
              real_t el_loss = std::min(el_val, eff_mmp * mmp_el_degrade);
              if (el_loss > 1e-10) {
                elastin_grid->ChangeConcentrationBy(el_si, -el_loss * coarse_w);
              }
            }
          }
          // Fibrin degradation (structural)
          if (do_fibrin_degrade && fibrin_grid) {
            size_t fb_si = sidx(fibrin_grid, idx);
            real_t fb_val = fibrin_grid->GetConcentration(fb_si);
            if (fb_val > 1e-10) {
              real_t fb_loss = std::min(fb_val, eff_mmp * mmp_fb_degrade);
              if (fb_loss > 1e-10) {
                fibrin_grid->ChangeConcentrationBy(fb_si, -fb_loss * coarse_w);
              }
            }
          }
        }
      }

      // -- Fibronectin serum source (structural) --
      if (do_fn) {
        real_t sv = stratum_grid->GetConcentration(idx);
        if (sv < 1.0 && fn_serum > 1e-10) {
          size_t fn_si = sidx(fn_grid, idx);
          real_t gate = std::max(static_cast<real_t>(0),
                                 static_cast<real_t>(1) - sv);
          // Carrying capacity limits serum FN like fibroblast deposition
          real_t fn_cap = sp->fibronectin_carrying_capacity;
          if (fn_cap > 0) {
            real_t fn_cur = fn_grid->GetConcentration(fn_si);
            gate *= std::max(static_cast<real_t>(0), 1.0 - fn_cur / fn_cap);
          }
          if (gate > 1e-10) {
            fn_grid->ChangeConcentrationBy(fn_si, fn_serum * gate * coarse_w);
          }
        }
      }

      // -- Scar accumulation (structural, reads fine inflammation) --
      if (do_scar) {
        real_t inflammation = infl_grid->GetConcentration(idx);
        if (inflammation > 1e-10) {
          scar_grid->ChangeConcentrationBy(
              sidx(scar_grid, idx), inflammation * scar_rate * coarse_w);
        }
      }

      // -- Keratinocyte MMP-1 at wound migration front --
      // Keratinocytes produce MMP-1 to degrade basement membrane during
      // re-epithelialization (Pilcher et al. 1997). Maximal at migration
      // front where stratum transitions (0 < sv < 1).
      if (mmp_kerat > 0) {
        real_t sv = stratum_grid->GetConcentration(idx);
        real_t edge = sv * (1.0 - sv) * 4.0;  // peaks at sv=0.5
        if (edge > 1e-10) {
          mmp_grid->ChangeConcentrationBy(idx, mmp_kerat * edge);
        }
      }
    }

    if (seeding_step) {
      seeded_ = true;
      if (sp->debug_wound) {
        std::cout << "[wound] biofilm seeded step=" << step << std::endl;
      }
    }
    timer.Print("fused_post");
  }

 private:
  GridContext::WoundMaskData mask_;
  bool seeded_ = false;
};

}  // namespace skibidy
}  // namespace bdm

#endif  // FUSED_POST_H_
