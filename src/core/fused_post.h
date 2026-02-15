#ifndef FUSED_POST_H_
#define FUSED_POST_H_

#include "infra/util.h"
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
//        MMP-collagen degradation, decorin TGF-beta sequestration,
//        fibronectin serum source.
//
// Only iterates epidermal wound voxels. The ordering within each voxel
// preserves dependencies: biofilm writes inflammation -> baseline writes
// inflammation -> MMP degrades ECM -> fibronectin serum -> scar reads
// inflammation.
//
// Reads/writes: Biofilm, Inflammation, Scar, MMP, Collagen, TGF-beta,
//               Fibronectin, Elastin (8 fields total)
// Tight coupling is intentional: ~7x speedup from single-pass iteration.
// ---------------------------------------------------------------------------
struct FusedWoundPostOp : public StandaloneOperationImpl {
  BDM_OP_HEADER(FusedWoundPostOp);

  void operator()() override {
    auto* sim = Simulation::GetActive();
    auto* sp = sim->GetParam()->Get<SimParam>();
    PerfTimer timer(sp->debug_perf);
    auto* scheduler = sim->GetScheduler();
    uint64_t step = GetGlobalStep(sim);
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
    bool do_decorin = sp->fibroblast_enabled &&
                      sp->decorin_sequestration_rate > 0;
    bool do_tissue_tgfb = sp->fibroblast_enabled &&
                          sp->tgfb_tissue_clearance > 0;
    bool do_lymphatic = sp->lymphatic_clearance_rate > 0;
    bool do_perfusion_clear = sp->perfusion_clearance_rate > 0;
    bool do_prommp = sp->mmp_enabled;
    bool do_age = sp->glucose_enabled && sp->diabetic_mode;
    bool do_senescence = sp->senescence_enabled;
    bool do_neuropeptide = sp->neuropathy_enabled;

    if (!do_biofilm && !do_baseline && !do_scar && !do_mmp && !do_fn &&
        !do_elastin_degrade && !do_decorin && !do_tissue_tgfb &&
        !do_lymphatic && !do_perfusion_clear && !do_prommp && !do_age &&
        !do_senescence && !do_neuropeptide) return;

    auto* rm = sim->GetResourceManager();

    // --- Grid pointers ---
    DiffusionGrid* stratum_grid = rm->GetDiffusionGrid(fields::kStratumId);
    GridContext ctx(stratum_grid, sp);

    // Lazy-init wound mask
    if (mask_.mask.empty()) {
      real_t z_max = sp->volume_z_cornified + ctx.box_len;
      mask_ = GridContext::ComputeWoundMask(ctx, z_max, sp);
    }

    // pH grid for MMP and biofilm modulation
    DiffusionGrid* ph_grid = rm->GetDiffusionGrid(fields::kPHId);

    DiffusionGrid* biofilm_grid = nullptr;
    if (do_biofilm) {
      biofilm_grid = rm->GetDiffusionGrid(fields::kBiofilmId);
    }

    DiffusionGrid* infl_grid = nullptr;
    if (do_biofilm || do_baseline || do_scar || do_age || do_senescence ||
        do_neuropeptide) {
      if (sp->split_inflammation_enabled) {
        infl_grid = rm->GetDiffusionGrid(fields::kProInflammatoryId);
      } else {
        infl_grid = rm->GetDiffusionGrid(fields::kInflammationId);
      }
    }

    DiffusionGrid* scar_grid = nullptr;
    if (do_scar) {
      scar_grid = rm->GetDiffusionGrid(fields::kScarId);
    }

    DiffusionGrid* mmp_grid = nullptr;
    DiffusionGrid* col_grid = nullptr;
    if (do_mmp) {
      mmp_grid = rm->GetDiffusionGrid(fields::kMMPId);
    }
    bool do_scar_collagen = do_scar && sp->fibroblast_enabled;
    if (do_mmp || do_decorin || do_scar_collagen) {
      col_grid = rm->GetDiffusionGrid(fields::kCollagenId);
    }

    DiffusionGrid* tgfb_grid = nullptr;
    if (do_decorin || do_tissue_tgfb || do_lymphatic || do_perfusion_clear) {
      tgfb_grid = rm->GetDiffusionGrid(fields::kTGFBetaId);
    }

    DiffusionGrid* vasc_grid = nullptr;
    if (do_tissue_tgfb || do_lymphatic || do_neuropeptide) {
      vasc_grid = rm->GetDiffusionGrid(fields::kVascularId);
    }

    // Perfusion (vascular) grid for venous TGF-beta clearance.
    // The Vascular field represents blood vessel density / perfusion recovery.
    DiffusionGrid* perf_grid = nullptr;
    if (do_perfusion_clear) {
      perf_grid = rm->GetDiffusionGrid(fields::kVascularId);
    }

    // Precompute lymphatic maturity (sigmoid of wound age with onset delay).
    // Lymphangiogenesis follows angiogenesis with ~4 day lag; lymphatic
    // vessels drain interstitial fluid carrying soluble cytokines.
    // (Paavonen et al. 2000, doi:10.1002/1097-0215)
    real_t lymph_maturity = 0;
    if (do_lymphatic) {
      real_t wound_h = wound_age * sim->GetParam()->simulation_time_step;
      lymph_maturity = static_cast<real_t>(1) /
          (static_cast<real_t>(1) +
           std::exp(-sp->lymphatic_sigmoid_k * (wound_h - sp->lymphatic_onset_h)));
    }

    DiffusionGrid* fn_grid = nullptr;
    if (do_fn) {
      fn_grid = rm->GetDiffusionGrid(fields::kFibronectinId);
    }

    DiffusionGrid* elastin_grid = nullptr;
    if (do_elastin_degrade) {
      elastin_grid = rm->GetDiffusionGrid(fields::kElastinId);
    }

    bool do_fibrin_degrade = sp->hemostasis_enabled && do_mmp;
    DiffusionGrid* fibrin_grid = nullptr;
    if (do_fibrin_degrade) {
      fibrin_grid = rm->GetDiffusionGrid(fields::kFibrinId);
    }

    // TIMP grid (explicit MMP inhibitor field)
    DiffusionGrid* timp_grid = nullptr;
    if (do_mmp) {
      timp_grid = rm->GetDiffusionGrid(fields::kTIMPId);
    }

    // Pro-MMP zymogen grid
    DiffusionGrid* prommp_grid = nullptr;
    if (do_prommp) {
      prommp_grid = rm->GetDiffusionGrid(fields::kProMMPId);
    }

    // AGE accumulator grid
    DiffusionGrid* age_grid = nullptr;
    if (do_age) {
      age_grid = rm->GetDiffusionGrid(fields::kAGEId);
    }

    // Senescence accumulator grid
    DiffusionGrid* sen_grid = nullptr;
    if (do_senescence) {
      sen_grid = rm->GetDiffusionGrid(fields::kSenescenceId);
    }

    // Nerve density grid (for neuropeptide effects)
    DiffusionGrid* nerve_grid = nullptr;
    if (do_neuropeptide) {
      nerve_grid = rm->GetDiffusionGrid(fields::kNerveId);
    }

    // New physics grid pointers
    DiffusionGrid* temp_grid = nullptr;
    bool do_temp = sp->temperature_enabled;
    if (do_temp) {
      temp_grid = rm->GetDiffusionGrid(fields::kTemperatureId);
    }
    DiffusionGrid* glucose_grid = nullptr;
    bool do_glucose = sp->glucose_enabled;
    bool do_glucose_age = do_glucose && sp->diabetic_mode;
    bool do_glucose_biofilm = do_glucose && do_biofilm;
    if (do_glucose) {
      glucose_grid = rm->GetDiffusionGrid(fields::kGlucoseId);
    }
    DiffusionGrid* no_grid = nullptr;
    bool do_no_antimicrobial = sp->nitric_oxide_enabled && do_biofilm;
    if (do_no_antimicrobial) {
      no_grid = rm->GetDiffusionGrid(fields::kNitricOxideId);
    }

    // Multi-resolution: structural grids may use coarser resolution.
    // Precompute fine-to-coarse index mapping once (all structural grids
    // share the same coarse resolution, so one table suffices).
    bool coarse = sp->grid_resolution_structural > 0 &&
                  sp->grid_resolution_structural != sp->grid_resolution;
    real_t coarse_w = 1.0;
    if (coarse && coarse_map_.empty()) {
      real_t rf = static_cast<real_t>(sp->grid_resolution);
      real_t rc = static_cast<real_t>(sp->grid_resolution_structural);
      coarse_w = (rc * rc * rc) / (rf * rf * rf);
      // Build lookup table using any structural grid
      DiffusionGrid* ref = col_grid ? col_grid :
          (scar_grid ? scar_grid : (fn_grid ? fn_grid : nullptr));
      if (ref) {
        coarse_map_.resize(ctx.n);
        for (size_t i = 0; i < ctx.n; i++) {
          Real3 pos = {ctx.X(i), ctx.Y(i), ctx.Z(i)};
          coarse_map_[i] = ref->GetBoxIndex(pos);
        }
      }
    } else if (coarse) {
      real_t rf = static_cast<real_t>(sp->grid_resolution);
      real_t rc = static_cast<real_t>(sp->grid_resolution_structural);
      coarse_w = (rc * rc * rc) / (rf * rf * rf);
    }
    auto sidx = [&](DiffusionGrid*, size_t fine) -> size_t {
      return coarse ? coarse_map_[fine] : fine;
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

    // Precompute dermal z-slice offset for scar stratum stamping.
    // Instead of col_grid->GetValue({x,y,dermal_z}) (coordinate lookup +
    // interpolation per voxel per step), use index arithmetic:
    // dermal_fine_idx = (idx % res2) + dermal_z_offset
    size_t res2 = ctx.res * ctx.res;
    size_t dermal_z_offset = 0;
    if (do_scar_collagen) {
      int dz_i = static_cast<int>((sp->dermal_fibroblast_depth - ctx.lo) /
                                   ctx.box_len);
      dz_i = std::max(0, std::min(dz_i, static_cast<int>(ctx.res) - 1));
      dermal_z_offset = static_cast<size_t>(dz_i) * res2;
    }

    // Collagen is used by MMP, decorin, and scar -- hoist the coarse index.
    bool need_col = do_mmp || do_decorin || do_scar_collagen || do_lymphatic;

    // --- Compact loop: epidermal wound voxels only ---
    for (size_t idx : mask_.epi_wound) {

      // pH alkalinity for MMP and biofilm modulation
      real_t ph_val = ph_grid ?
          ph_grid->GetConcentration(sidx(ph_grid, idx)) : 0;

      // Hoisted collagen index (shared by MMP, decorin, scar)
      size_t col_si = need_col ? sidx(col_grid, idx) : 0;

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
            // Temperature Q10: warmer wound promotes bacterial growth
            if (do_temp && temp_grid) {
              real_t temp_val = temp_grid->GetConcentration(idx);
              real_t temp_c = temp_val * 37.0;  // denormalize to Celsius
              eff_growth *= std::pow(sp->temperature_q10_biofilm,
                                     (temp_c - 37.0) / 10.0);
            }
            // Glucose fuels bacterial metabolism
            if (do_glucose_biofilm && glucose_grid) {
              real_t gluc = glucose_grid->GetConcentration(idx);
              eff_growth *= (1.0 + sp->glucose_bacterial_consumption * gluc);
            }
            // NO antimicrobial: reactive nitrogen species suppress growth
            if (do_no_antimicrobial && no_grid) {
              real_t no_val = no_grid->GetConcentration(idx);
              eff_growth *= std::max(0.0, 1.0 - sp->no_antimicrobial_factor * no_val);
            }
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
          real_t eff_baseline = baseline_rate;
          // Glucose-driven AGE/RAGE inflammation: mechanistic diabetic pathway
          if (do_glucose_age && glucose_grid) {
            real_t gluc = glucose_grid->GetConcentration(idx);
            eff_baseline += sp->glucose_age_inflammation * gluc;
          }
          infl_grid->ChangeConcentrationBy(idx, eff_baseline * gate);
        }
      }

      // -- MMP-mediated ECM degradation --
      // MMP (fine grid) degrades structural ECM fields at coarse resolution.
      // Alkaline pH amplifies MMP proteolytic activity (optimal pH ~7.0-8.0).
      real_t total_degraded = 0;
      if (do_mmp) {
        real_t mmp_val = mmp_grid->GetConcentration(idx);
        if (mmp_val > 1e-10) {
          real_t mmp_ph = 1.0 + sp->ph_mmp_boost * ph_val;
          // Temperature Q10: enzyme kinetics scale with temperature
          real_t mmp_temp = 1.0;
          if (do_temp && temp_grid) {
            real_t temp_val = temp_grid->GetConcentration(idx);
            real_t temp_c = temp_val * 37.0;
            mmp_temp = std::pow(sp->temperature_q10_mmp,
                                (temp_c - 37.0) / 10.0);
          }
          real_t eff_mmp = mmp_val * mmp_ph * mmp_temp;
          // Collagen degradation (structural, uses hoisted col_si)
          // AGE cross-linking: glycated collagen is resistant to MMP cleavage
          // (Mott et al. 1997, doi:10.1038/ki.1997.455)
          real_t col_val = col_grid->GetConcentration(col_si);
          if (col_val > 1e-10) {
            real_t col_degrade_rate = mmp_col_degrade;
            if (do_age && age_grid) {
              real_t age_val = age_grid->GetConcentration(idx);
              col_degrade_rate *= std::max(static_cast<real_t>(0),
                  1.0 - sp->age_collagen_crosslink * age_val);
            }
            real_t degrade = std::min(col_val, eff_mmp * col_degrade_rate);
            if (degrade > 1e-10) {
              col_grid->ChangeConcentrationBy(col_si, -degrade * coarse_w);
              total_degraded += degrade;
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
                total_degraded += fn_loss;
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
                total_degraded += el_loss;
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
                total_degraded += fb_loss;
              }
            }
          }

          // Matrikine positive feedback: ECM degradation fragments
          // (endostatin, PGE2) stimulate MMP transcription via NF-kB.
          // Self-limiting: depends on available ECM substrate.
          // Visse & Nagase 2003 (doi:10.1161/01.res.0000070112.80711.3d)
          if (sp->matrikine_mmp_boost > 0 && total_degraded > 1e-10) {
            mmp_grid->ChangeConcentrationBy(
                idx, sp->matrikine_mmp_boost * total_degraded);
          }
        }
      }

      // -- Decorin-mediated TGF-beta sequestration --
      // Decorin, a small leucine-rich proteoglycan that co-deposits with
      // collagen fibrils, directly binds and neutralizes active TGF-beta1.
      // As collagen matures, decorin accumulates proportionally, creating
      // a negative feedback that resolves TGF-beta in late proliferative
      // and remodeling phases without affecting early activation timing.
      // (Yamaguchi et al. 1990, doi:10.1038/346281a0;
      //  Border & Ruoslahti 1992, doi:10.1038/357568a0)
      if (do_decorin) {
        real_t col_val = col_grid->GetConcentration(col_si);
        if (col_val > 1e-10) {
          real_t tgfb_val = tgfb_grid->GetConcentration(idx);
          if (tgfb_val > 1e-10) {
            real_t sink = sp->decorin_sequestration_rate * col_val * tgfb_val;
            tgfb_grid->ChangeConcentrationBy(idx, -std::min(sink, tgfb_val));
          }
        }
      }

      // -- Tissue receptor-mediated TGF-beta clearance --
      // Cells expressing TGF-beta receptors (TbRII/TbRI) internalize
      // active ligand via clathrin-mediated endocytosis. Clearance rate
      // scales with local tissue density, derived from epithelial
      // coverage (stratum) and vascular recovery. Open wound bed: low
      // density, low clearance. Healed tissue: high density, high
      // clearance.
      if (do_tissue_tgfb) {
        real_t tgfb_val = tgfb_grid->GetConcentration(idx);
        if (tgfb_val > 1e-10) {
          real_t sv = stratum_grid->GetConcentration(idx);
          real_t vasc = vasc_grid ? vasc_grid->GetConcentration(idx) : 0;
          real_t tissue_density = std::max(
              std::max(static_cast<real_t>(0), sv),
              std::max(static_cast<real_t>(0), vasc));
          tissue_density = std::min(static_cast<real_t>(1), tissue_density);
          if (tissue_density > 1e-10) {
            real_t sink = sp->tgfb_tissue_clearance * tissue_density * tgfb_val;
            tgfb_grid->ChangeConcentrationBy(idx, -std::min(sink, tgfb_val));
          }
        }
      }

      // -- Lymphatic drainage of TGF-beta --
      // Recovering lymphatic vessels drain interstitial fluid carrying
      // soluble TGF-beta. Lymphatic density scales with vascular recovery
      // gated by a sigmoid onset delay (lymphangiogenesis lags
      // angiogenesis by ~4 days). Provides a clearance pathway that is
      // absent in the early wound bed and ramps up during the
      // proliferative phase.
      // (Paavonen et al. 2000, doi:10.1016/s0002-9440(10)65021-3;
      //  Kataru et al. 2009, doi:10.1182/blood-2008-09-176776)
      // Lymphatic drainage of soluble TGF-beta. New lymphatic capillaries
      // form in granulation tissue (lymphangiogenesis) proportional to local
      // collagen density, with a sigmoid onset lag (~4 days post-wounding).
      // This couples lymphatic clearance to granulation tissue maturation
      // rather than pre-existing vasculature.
      // (Paavonen et al. 2000; Kataru et al. 2009)
      if (do_lymphatic && lymph_maturity > 1e-6) {
        real_t tgfb_val = tgfb_grid->GetConcentration(idx);
        if (tgfb_val > 1e-10) {
          real_t granulation = std::max(static_cast<real_t>(0),
              col_grid->GetConcentration(col_si));
          real_t sink = sp->lymphatic_clearance_rate * lymph_maturity *
                        granulation * tgfb_val;
          tgfb_grid->ChangeConcentrationBy(idx, -std::min(sink, tgfb_val));
        }
      }

      // -- Venous perfusion clearance of TGF-beta --
      // Blood flow through granulation tissue vasculature carries away
      // soluble TGF-beta1 via venous drainage. Clearance scales with local
      // perfusion (vascular density), providing a natural feedback: as
      // angiogenesis restores blood supply, growth factor clearance increases.
      // (Bramhall et al. 1999)
      if (do_perfusion_clear) {
        real_t tgfb_val = tgfb_grid->GetConcentration(idx);
        if (tgfb_val > 1e-10) {
          // Perfusion grid is structural; use sidx for potential multi-res
          real_t perf = perf_grid ?
              perf_grid->GetConcentration(sidx(perf_grid, idx)) : 0;
          if (perf > 1e-10) {
            real_t sink = sp->perfusion_clearance_rate * perf * tgfb_val;
            tgfb_grid->ChangeConcentrationBy(idx, -std::min(sink, tgfb_val));
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

      // -- Scar accumulation --
      // Collagen-driven (mechanistic): scar mirrors collagen, stratum
      // gets +5 offset for scar visualization. Fallback to inflammation
      // integral when fibroblasts disabled.
      if (do_scar) {
        if (do_scar_collagen) {
          // Mirror DERMAL collagen to epidermal scar field. Collagen is
          // deposited by myofibroblasts in the dermis; the epidermal scar
          // reflects the underlying dermal ECM (surface scar appearance).
          size_t dermal_fine = (idx % res2) + dermal_z_offset;
          size_t dermal_ci = sidx(col_grid, dermal_fine);
          real_t dermal_col = col_grid->GetConcentration(dermal_ci);
          size_t sc_si = sidx(scar_grid, idx);
          real_t sc_cur = scar_grid->GetConcentration(sc_si);
          real_t sc_delta = dermal_col - sc_cur;
          if (std::abs(sc_delta) > 1e-10) {
            scar_grid->ChangeConcentrationBy(sc_si, sc_delta);
          }
          // Stamp stratum +5 for scar visualization
          if (dermal_col > sp->scar_collagen_threshold) {
            real_t st_cur = stratum_grid->GetConcentration(idx);
            if (st_cur >= 1.0 && st_cur < 4.5) {
              real_t scar_val = std::round(st_cur) + 5.0;
              if (scar_val > st_cur) {
                stratum_grid->ChangeConcentrationBy(idx, scar_val - st_cur);
              }
            }
          }
        } else {
          // Inflammation integral fallback (fibroblasts disabled)
          real_t inflammation = infl_grid->GetConcentration(idx);
          if (inflammation > 1e-10) {
            scar_grid->ChangeConcentrationBy(
                sidx(scar_grid, idx), inflammation * scar_rate * coarse_w);
          }
        }
      }

      // -- Keratinocyte pro-MMP-1 at wound migration front --
      // Keratinocytes produce pro-MMP-1 to degrade basement membrane during
      // re-epithelialization (Pilcher et al. 1997). Deposited as zymogen;
      // activation in the pro-MMP cascade below.
      if (mmp_kerat > 0 && prommp_grid) {
        real_t sv = stratum_grid->GetConcentration(idx);
        real_t edge = sv * (1.0 - sv) * 4.0;  // peaks at sv=0.5
        if (edge > 1e-10) {
          prommp_grid->ChangeConcentrationBy(idx, mmp_kerat * edge);
        }
      }

      // -- Keratinocyte TIMP-1 at wound migration front --
      // Saarialho-Kere et al. 1995 (doi:10.1172/JCI117883)
      if (do_mmp && timp_grid && sp->timp_keratinocyte_rate > 0) {
        real_t sv = stratum_grid->GetConcentration(idx);
        real_t edge = sv * (1.0 - sv) * 4.0;
        if (edge > 1e-10) {
          real_t rate = sp->timp_keratinocyte_rate;
          if (sp->diabetic_mode) rate *= sp->diabetic_timp_production_factor;
          timp_grid->ChangeConcentrationBy(idx, rate * edge);
        }
      }

      // -- MMP-TIMP second-order neutralization --
      // Active MMP and TIMP form irreversible 1:1 stoichiometric complexes.
      // This replaces the implicit TIMP decay constant with explicit kinetics.
      // Brew et al. 2000 (doi:10.1016/S0167-4838(99)00252-5)
      if (do_mmp && timp_grid && sp->mmp_timp_inhibition_rate > 0) {
        real_t mmp_val = mmp_grid->GetConcentration(idx);
        real_t timp_val = timp_grid->GetConcentration(idx);
        if (mmp_val > 1e-10 && timp_val > 1e-10) {
          real_t neutralize = sp->mmp_timp_inhibition_rate * mmp_val * timp_val;
          real_t loss = std::min(neutralize, std::min(mmp_val, timp_val));
          if (loss > 1e-10) {
            mmp_grid->ChangeConcentrationBy(idx, -loss);
            timp_grid->ChangeConcentrationBy(idx, -loss);
          }
        }
      }

      // -- Pro-MMP zymogen activation cascade --
      // MMPs are secreted as inactive zymogens (pro-MMPs). Extracellular
      // activation occurs via: (1) basal plasmin-mediated cleavage of the
      // cysteine switch, and (2) autocatalytic MMP-3 (stromelysin) positive
      // feedback on other pro-MMPs.
      // Visse & Nagase 2003 (doi:10.1161/01.res.0000070112.80711.3d)
      if (do_prommp && prommp_grid) {
        real_t prommp_val = prommp_grid->GetConcentration(idx);
        if (prommp_val > 1e-10) {
          real_t mmp_val = mmp_grid ? mmp_grid->GetConcentration(idx) : 0;
          // Basal plasmin activation + autocatalytic (MMP-3 amplifies)
          real_t activation = sp->prommp_activation_rate * prommp_val +
                              sp->prommp_autocatalytic_rate * prommp_val * mmp_val;
          activation = std::min(activation, prommp_val);
          if (activation > 1e-10) {
            prommp_grid->ChangeConcentrationBy(idx, -activation);
            if (mmp_grid) {
              mmp_grid->ChangeConcentrationBy(idx, activation);
            }
          }
        }
      }

      // -- AGE formation from glucose (diabetic only) --
      // Non-enzymatic glycation: glucose reacts with lysine/arginine residues
      // on ECM proteins via Amadori rearrangement to form irreversible AGEs.
      // Rate proportional to glucose concentration and exposure time.
      // Brownlee 2005 (doi:10.2337/diabetes.54.6.1615)
      if (do_age && age_grid && glucose_grid) {
        real_t gluc = glucose_grid->GetConcentration(idx);
        if (gluc > 1e-10) {
          real_t age_formation = sp->glucose_age_rate * gluc;
          age_grid->ChangeConcentrationBy(idx, age_formation);
        }
        // RAGE-mediated NF-kB inflammation: AGEs bind RAGE receptor on
        // macrophages and endothelial cells, activating NF-kB and sustaining
        // M1 polarization. (Chavakis et al. 2004, doi:10.1016/j.micinf.2004.08.004)
        real_t age_val = age_grid->GetConcentration(idx);
        if (age_val > 1e-10 && infl_grid) {
          real_t sv = stratum_grid->GetConcentration(idx);
          real_t gate = std::max(static_cast<real_t>(0),
                                 static_cast<real_t>(1) - sv);
          if (gate > 1e-10) {
            infl_grid->ChangeConcentrationBy(idx,
                sp->age_rage_inflammation * age_val * gate);
          }
        }
      }

      // -- Senescence accumulation and SASP output --
      // Senescent cells accumulate from DNA damage (wound-induced), oxidative
      // stress (inflammation-driven), and glycative stress (AGE in diabetic).
      // SASP: senescent cells secrete pro-inflammatory cytokines, MMPs, and
      // TGF-beta, creating a paracrine feedback that impairs wound healing.
      // Demaria et al. 2014 (doi:10.1016/j.devcel.2014.11.012)
      // Coppe et al. 2008 (doi:10.1371/journal.pbio.0060301)
      if (do_senescence && sen_grid) {
        size_t sen_si = sidx(sen_grid, idx);

        // Accumulation: wound-induced DNA damage
        real_t sv = stratum_grid->GetConcentration(idx);
        real_t wound_gate = std::max(static_cast<real_t>(0),
                                     static_cast<real_t>(1) - sv);
        real_t accum = sp->senescence_wound_rate * wound_gate;

        // Accumulation: inflammation-driven (ROS, NF-kB)
        if (infl_grid) {
          real_t infl_val = infl_grid->GetConcentration(idx);
          if (infl_val > 1e-10) {
            accum += sp->senescence_infl_rate * infl_val;
          }
        }

        // Accumulation: AGE-driven (diabetic glycative stress)
        if (do_age && age_grid) {
          real_t age_val = age_grid->GetConcentration(idx);
          if (age_val > 1e-10) {
            accum += sp->senescence_age_rate * age_val;
          }
        }

        // Diabetic mode: accelerated senescence
        if (sp->diabetic_mode) {
          accum *= sp->diabetic_senescence_factor;
        }

        // Senolytic clearance (treatment pathway)
        real_t sen_val = sen_grid->GetConcentration(sen_si);
        if (sp->senolytic_clearance_rate > 0 && sen_val > 1e-10) {
          accum -= sp->senolytic_clearance_rate * sen_val;
        }

        if (std::abs(accum) > 1e-10) {
          real_t new_val = sen_val + accum;
          if (new_val < 0) accum = -sen_val;
          sen_grid->ChangeConcentrationBy(sen_si, accum * coarse_w);
        }

        // SASP output (only when senescent cells are present)
        sen_val = sen_grid->GetConcentration(sen_si);
        if (sen_val > 1e-10) {
          // Pro-inflammatory cytokines (IL-6, IL-8)
          if (infl_grid && sp->sasp_inflammation_rate > 0) {
            infl_grid->ChangeConcentrationBy(idx,
                sp->sasp_inflammation_rate * sen_val);
          }
          // MMP-3/MMP-9 (as pro-MMP zymogen)
          if (prommp_grid && sp->sasp_mmp_rate > 0) {
            prommp_grid->ChangeConcentrationBy(idx,
                sp->sasp_mmp_rate * sen_val);
          }
          // TGF-beta1 (fibrotic)
          if (tgfb_grid && sp->sasp_tgfb_rate > 0) {
            tgfb_grid->ChangeConcentrationBy(idx,
                sp->sasp_tgfb_rate * sen_val);
          }
        }
      }

      // -- Neuropeptide effects on wound healing --
      // Sensory nerve terminals release substance P and CGRP. Local nerve
      // density determines neuropeptide flux: denervated wound = reduced
      // neurogenic signaling = impaired healing.
      // Suvas 2017 (doi:10.4049/jimmunol.1601751)
      if (do_neuropeptide && nerve_grid) {
        size_t nerve_si = sidx(nerve_grid, idx);
        real_t nerve_val = nerve_grid->GetConcentration(nerve_si);
        if (nerve_val > 1e-10) {
          // Substance P: neurogenic inflammation (mast cell degranulation)
          if (infl_grid && sp->substance_p_inflammation > 0) {
            infl_grid->ChangeConcentrationBy(idx,
                sp->substance_p_inflammation * nerve_val);
          }
          // CGRP: vasodilation boosts local perfusion
          if (vasc_grid && sp->cgrp_vasodilation > 0) {
            real_t vasc_val = vasc_grid->GetConcentration(idx);
            real_t vasc_target = sp->perfusion_basal;
            if (vasc_val < vasc_target) {
              real_t boost = sp->cgrp_vasodilation * nerve_val;
              real_t headroom = vasc_target - vasc_val;
              real_t gain = std::min(headroom, boost);
              if (gain > 1e-10) {
                vasc_grid->ChangeConcentrationBy(idx, gain);
              }
            }
          }
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
  std::vector<size_t> coarse_map_;  // fine-to-coarse index lookup
  bool seeded_ = false;
};

}  // namespace skibidy
}  // namespace bdm

#endif  // FUSED_POST_H_
