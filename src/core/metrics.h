#ifndef METRICS_H_
#define METRICS_H_

#include <chrono>
#include <cmath>
#include <fstream>
#include <memory>
#include <string>

#include "tissue/keratinocyte.h"
#include "immune/immune_cell.h"
#include "fibroblast/fibroblast.h"
#include "tumor/tumor_cell.h"
#include "core/field_names.h"
#include "core/pde.h"
#include "infra/sim_param.h"
#include "infra/util.h"
#include "core/operation/operation.h"

namespace bdm {
namespace skibidy {

// ---------------------------------------------------------------------------
// MetricsExporter -- standalone operation that writes a CSV row of simulation
// metrics every metrics_interval steps.
//
// 53 columns -- see header string below for full list
// Note: blood, burn, and pressure modules modify existing fields (perfusion,
// inflammation, ROS) rather than creating new grids. Their effects show in
// existing columns (mean_perfusion_wound, mean_infl_wound, mean_ros_wound).
// ---------------------------------------------------------------------------
struct MetricsExporter : public StandaloneOperationImpl {
  BDM_OP_HEADER(MetricsExporter);

  void operator()() override {
    auto* sim = Simulation::GetActive();
    auto* scheduler = sim->GetScheduler();
    auto* sp = sim->GetParam()->Get<SimParam>();
    PerfTimer timer(sp->debug_perf);

    uint64_t step = GetGlobalStep(sim);

    // Write header on first call
    if (!ofs_) {
      std::string path = sim->GetOutputDir() + "/metrics.csv";
      ofs_ = std::make_shared<std::ofstream>(
          path, std::ios::out | std::ios::trunc);
      *ofs_ << "step,time_h,time_days,n_agents,n_stem,n_ta,n_g0,"
            << "n_basal,n_spinous,n_granular,n_cornified,"
            << "n_neutrophils,n_macrophages,"
            << "wound_closure_pct,mean_o2_wound,mean_ca_wound,mean_infl_wound,"
            << "scar_magnitude,mean_anti_infl_wound,"
            << "n_fibroblasts,n_activated_fibroblasts,n_myofibroblasts,mean_tgfb_wound,mean_collagen_wound,"
            << "mean_perfusion_wound,"
            << "n_tumor_cells,n_tumor_cycling,tumor_field_cells,"
            << "mean_biofilm_wound,mean_vegf_wound,"
            << "mean_mmp_wound,mean_fibronectin_wound,"
            << "mean_elastin_wound,mean_hyaluronan_wound,"
            << "mean_dermis_papillary,mean_dermis_reticular,mean_dermis_hypodermis,"
            << "mean_ph_wound,"
            << "mean_fibrin_wound,"
            << "mean_ecm_quality,mean_tissue_viability,"
            << "mean_temperature_wound,mean_glucose_wound,"
            << "mean_lactate_wound,mean_no_wound,mean_ros_wound,"
            << "mean_stiffness_wound,mean_lymphatic_wound,mean_edema_wound,mean_voltage_wound,"
            << "mean_tnf_alpha_wound,mean_il6_wound,mean_cartilage_wound,mean_synovial_wound,"
            << "mean_tcell_wound,mean_bone_wound,"
            << "mean_scab_wound"
            << std::endl;
    }

    if (sp->metrics_interval <= 0) return;
    if (step % static_cast<uint64_t>(sp->metrics_interval) != 0) return;

    // Throughput report (every metrics interval, ~1h sim time)
    if (sp->debug_perf) {
      auto now = std::chrono::high_resolution_clock::now();
      if (last_step_ > 0) {
        double secs = std::chrono::duration<double>(now - last_time_).count();
        uint64_t delta = step - last_step_;
        if (secs > 0) {
          std::cout << "[perf] step=" << step
                    << " rate=" << static_cast<int>(delta / secs)
                    << " steps/s" << std::endl;
        }
      }
      last_step_ = step;
      last_time_ = now;
    }

    auto* rm = sim->GetResourceManager();

    // --- Agent metrics (single pass) ---
    int n_agents = 0, n_stem = 0, n_ta = 0, n_g0 = 0;
    int n_basal = 0, n_spinous = 0, n_granular = 0, n_cornified = 0;
    int n_neutrophils = 0, n_macrophages = 0;
    int n_fibroblasts = 0, n_activated_fibroblasts = 0, n_myofibroblasts = 0;
    int n_tumor_cells = 0;
    int n_tumor_cycling = 0;

    rm->ForEachAgent([&](Agent* agent) {
      // Keratinocytes first: ~90% of agents, avoids 3 failed casts per cell.
      auto* cell = dynamic_cast<Keratinocyte*>(agent);
      if (cell) {
        n_agents++;
        if (cell->IsStem()) {
          n_stem++;
        } else if (cell->GetCyclePhase() == kG0) {
          n_g0++;
        } else {
          n_ta++;
        }
        switch (cell->GetStratum()) {
          case kBasal:     n_basal++; break;
          case kSpinous:   n_spinous++; break;
          case kGranular:  n_granular++; break;
          case kCornified: n_cornified++; break;
          default: break;
        }
        return;
      }

      // Immune cells
      auto* immune = dynamic_cast<ImmuneCell*>(agent);
      if (immune) {
        if (immune->GetImmuneCellType() == kNeutrophil) {
          n_neutrophils++;
        } else {
          n_macrophages++;
        }
        return;
      }

      // Fibroblasts
      auto* fibro = dynamic_cast<Fibroblast*>(agent);
      if (fibro) {
        n_fibroblasts++;
        auto fstate = fibro->GetFibroblastState();
        if (fstate == kFibroActivated || fstate == kMyofibroblast) {
          n_activated_fibroblasts++;
        }
        if (fstate == kMyofibroblast) {
          n_myofibroblasts++;
        }
        return;
      }

      // Tumor cells (skip agents pending deferred removal after handoff)
      auto* tumor = dynamic_cast<TumorCell*>(agent);
      if (tumor) {
        if (tumor->GetCyclePhase() == kG0 &&
            sp->tumor.handoff_delay > 0 &&
            tumor->GetG0Steps() > sp->tumor.handoff_delay) {
          return;
        }
        n_tumor_cells++;
        if (tumor->GetCyclePhase() != kG0) n_tumor_cycling++;
        return;
      }
    });

    // --- Wound field metrics ---
    real_t wound_closure = 0;
    real_t mean_o2 = 0;
    real_t mean_ca = 0;
    real_t mean_infl = 0;
    real_t scar_magnitude = 0;
    real_t mean_anti_infl = 0;
    real_t mean_tgfb = 0;
    real_t mean_collagen = 0;
    real_t mean_perfusion = 0;
    real_t mean_biofilm = 0;
    real_t mean_vegf = 0;
    real_t mean_mmp = 0;
    real_t mean_fn = 0;
    real_t mean_elastin = 0;
    real_t mean_ha = 0;
    real_t mean_dermis_pap = 0, mean_dermis_ret = 0, mean_dermis_hyp = 0;
    real_t mean_ph = 0;
    real_t mean_fibrin = 0;
    real_t mean_ecm_quality = 0;
    real_t mean_tissue_viability = 0;
    real_t mean_temperature = 0;
    real_t mean_glucose = 0;
    real_t mean_lactate = 0;
    real_t mean_no = 0;
    real_t mean_ros = 0;
    real_t mean_stiffness = 0;
    real_t mean_lymphatic = 0;
    real_t mean_edema = 0;
    real_t mean_voltage = 0;
    real_t mean_tnf_alpha = 0;
    real_t mean_il6 = 0;
    real_t mean_cartilage = 0;
    real_t mean_synovial = 0;
    real_t mean_tcell = 0;
    real_t mean_bone = 0;
    real_t mean_scab = 0;

    if (sp->wound.enabled) {
      auto* stratum_grid = rm->GetDiffusionGrid(fields::kStratumId);
      auto* o2_grid = rm->GetDiffusionGrid(fields::kOxygenId);
      auto* ca_grid = rm->GetDiffusionGrid(fields::kCalciumId);
      auto* scar_grid = rm->GetDiffusionGrid(fields::kScarId);
      auto* vasc_grid = rm->GetDiffusionGrid(fields::kVascularId);

      DiffusionGrid* tgfb_grid = nullptr;
      DiffusionGrid* col_grid = nullptr;
      if (sp->fibroblast.enabled) {
        tgfb_grid = rm->GetDiffusionGrid(fields::kTGFBetaId);
        col_grid = rm->GetDiffusionGrid(fields::kCollagenId);
      }

      DiffusionGrid* biofilm_grid = nullptr;
      if (sp->biofilm.enabled) {
        biofilm_grid = rm->GetDiffusionGrid(fields::kBiofilmId);
      }
      DiffusionGrid* vegf_grid = nullptr;
      if (sp->angiogenesis.enabled) {
        vegf_grid = rm->GetDiffusionGrid(fields::kVEGFId);
      }
      DiffusionGrid* mmp_grid = nullptr;
      if (sp->mmp.enabled) {
        mmp_grid = rm->GetDiffusionGrid(fields::kMMPId);
      }
      DiffusionGrid* fn_grid = nullptr;
      if (sp->fibronectin.enabled) {
        fn_grid = rm->GetDiffusionGrid(fields::kFibronectinId);
      }
      DiffusionGrid* elastin_grid = nullptr;
      if (sp->elastin.enabled) {
        elastin_grid = rm->GetDiffusionGrid(fields::kElastinId);
      }
      DiffusionGrid* ha_grid = nullptr;
      if (sp->hyaluronan.enabled) {
        ha_grid = rm->GetDiffusionGrid(fields::kHyaluronanId);
      }
      DiffusionGrid* dermis_grid = nullptr;
      if (sp->dermis.enabled) {
        dermis_grid = rm->GetDiffusionGrid(fields::kDermisId);
      }
      DiffusionGrid* ph_grid = rm->GetDiffusionGrid(fields::kPHId);
      DiffusionGrid* fibrin_grid = nullptr;
      if (sp->hemostasis.enabled) {
        fibrin_grid = rm->GetDiffusionGrid(fields::kFibrinId);
      }
      DiffusionGrid* temp_grid_m = nullptr;
      if (sp->temperature.enabled) {
        temp_grid_m = rm->GetDiffusionGrid(fields::kTemperatureId);
      }
      DiffusionGrid* glucose_grid_m = nullptr;
      if (sp->glucose_mod.enabled) {
        glucose_grid_m = rm->GetDiffusionGrid(fields::kGlucoseId);
      }
      DiffusionGrid* lactate_grid_m = nullptr;
      if (sp->lactate.enabled) {
        lactate_grid_m = rm->GetDiffusionGrid(fields::kLactateId);
      }
      DiffusionGrid* no_grid_m = nullptr;
      if (sp->nitric_oxide.enabled) {
        no_grid_m = rm->GetDiffusionGrid(fields::kNitricOxideId);
      }
      DiffusionGrid* ros_grid_m = nullptr;
      if (sp->ros.enabled) {
        ros_grid_m = rm->GetDiffusionGrid(fields::kROSId);
      }
      DiffusionGrid* stiff_grid_m = nullptr;
      if (sp->mechanotransduction.enabled) {
        stiff_grid_m = rm->GetDiffusionGrid(fields::kStiffnessId);
      }
      DiffusionGrid* lymph_grid_m = nullptr;
      DiffusionGrid* edema_grid_m = nullptr;
      if (sp->lymphatic.enabled) {
        lymph_grid_m = rm->GetDiffusionGrid(fields::kLymphaticId);
        edema_grid_m = rm->GetDiffusionGrid(fields::kEdemaId);
      }
      DiffusionGrid* volt_grid_m = nullptr;
      if (sp->bioelectric.enabled) {
        volt_grid_m = rm->GetDiffusionGrid(fields::kVoltageId);
      }
      DiffusionGrid* scab_grid_m = nullptr;
      if (sp->scab.enabled) {
        scab_grid_m = rm->GetDiffusionGrid(fields::kScabId);
      }
      DiffusionGrid* tnf_grid_m = nullptr;
      DiffusionGrid* il6_grid_m = nullptr;
      DiffusionGrid* cart_grid_m = nullptr;
      DiffusionGrid* syn_grid_m = nullptr;
      DiffusionGrid* tcell_grid_m = nullptr;
      DiffusionGrid* bone_grid_m = nullptr;
      if (sp->ra.enabled) {
        tnf_grid_m = rm->GetDiffusionGrid(fields::kTNFAlphaId);
        il6_grid_m = rm->GetDiffusionGrid(fields::kIL6Id);
        cart_grid_m = rm->GetDiffusionGrid(fields::kCartilageId);
        syn_grid_m = rm->GetDiffusionGrid(fields::kSynovialFluidId);
        tcell_grid_m = rm->GetDiffusionGrid(fields::kTCellDensityId);
        bone_grid_m = rm->GetDiffusionGrid(fields::kBoneId);
      }

      DiffusionGrid* infl_grid = nullptr;
      DiffusionGrid* anti_grid = nullptr;
      if (sp->inflammation.split_inflammation_enabled) {
        infl_grid = rm->GetDiffusionGrid(fields::kProInflammatoryId);
        anti_grid = rm->GetDiffusionGrid(fields::kAntiInflammatoryId);
      } else {
        infl_grid = rm->GetDiffusionGrid(fields::kInflammationId);
      }

      GridContext ctx(stratum_grid, sp);

      // Multi-resolution: structural grids may use coarser resolution
      bool coarse = sp->grid_resolution_structural > 0 &&
                    sp->grid_resolution_structural != sp->grid_resolution;
      auto sidx = [&](DiffusionGrid* g, size_t fine) -> size_t {
        return coarse ? GridContext::CoarseIndex(fine, ctx, g) : fine;
      };

      // Lazy-init wound mask (computed once, reused every metrics interval)
      if (mask_.mask.empty()) {
        real_t z_max = sp->volume_z_cornified + ctx.box_len;
        mask_ = GridContext::ComputeWoundMask(ctx, z_max, sp);
      }

      int total_voxels = 0, filled_voxels = 0, dermal_voxels = 0;
      real_t sum_o2 = 0, sum_ca = 0, sum_infl = 0, sum_scar = 0;
      real_t sum_anti = 0, sum_tgfb = 0, sum_col = 0, sum_perf = 0;
      real_t sum_biofilm = 0, sum_vegf = 0;
      real_t sum_mmp = 0, sum_fn = 0;
      real_t sum_elastin = 0, sum_ha = 0;
      real_t sum_ph = 0;
      real_t sum_fibrin = 0;
      real_t sum_ecm = 0;
      real_t sum_viability = 0;
      real_t sum_temp = 0, sum_gluc = 0, sum_lac = 0, sum_no = 0, sum_ros = 0;
      real_t sum_stiff = 0, sum_lymph = 0, sum_edema = 0, sum_volt = 0;
      real_t sum_tnf = 0, sum_il6 = 0, sum_cart = 0, sum_syn = 0;
      real_t sum_tcell = 0, sum_bone = 0;
      real_t sum_scab = 0;

      // Compact loops: iterate only wound voxels via precomputed index lists.
      for (size_t idx : mask_.epi_wound) {
        total_voxels++;
        // Wound void = -1, kBasal = 0, kSpinous = 1, ...
        // Any non-void voxel (stratum >= 0) counts as re-epithelialized.
        if (stratum_grid->GetConcentration(idx) > -0.5) {
          filled_voxels++;
        }
        sum_o2 += o2_grid->GetConcentration(idx);
        sum_ca += ca_grid->GetConcentration(idx);
        sum_infl += infl_grid->GetConcentration(idx);
        sum_scar += scar_grid->GetConcentration(sidx(scar_grid, idx));
        if (anti_grid) sum_anti += anti_grid->GetConcentration(idx);
        if (tgfb_grid) sum_tgfb += tgfb_grid->GetConcentration(idx);
        if (col_grid) sum_col += col_grid->GetConcentration(sidx(col_grid, idx));
        if (biofilm_grid) sum_biofilm += biofilm_grid->GetConcentration(sidx(biofilm_grid, idx));
        if (vegf_grid) sum_vegf += vegf_grid->GetConcentration(idx);
        if (mmp_grid) sum_mmp += mmp_grid->GetConcentration(idx);
        if (fn_grid) sum_fn += fn_grid->GetConcentration(sidx(fn_grid, idx));
        if (ph_grid) sum_ph += ph_grid->GetConcentration(sidx(ph_grid, idx));
        if (fibrin_grid) sum_fibrin += fibrin_grid->GetConcentration(sidx(fibrin_grid, idx));
        if (g_ecm_quality) sum_ecm += g_ecm_quality->GetConcentration(idx);
        if (g_tissue_viability) sum_viability += g_tissue_viability->GetConcentration(idx);
        if (temp_grid_m) sum_temp += temp_grid_m->GetConcentration(idx);
        if (glucose_grid_m) sum_gluc += glucose_grid_m->GetConcentration(idx);
        if (lactate_grid_m) sum_lac += lactate_grid_m->GetConcentration(idx);
        if (no_grid_m) sum_no += no_grid_m->GetConcentration(idx);
        if (ros_grid_m) sum_ros += ros_grid_m->GetConcentration(idx);
        if (volt_grid_m) sum_volt += volt_grid_m->GetConcentration(idx);
        if (scab_grid_m) sum_scab += scab_grid_m->GetConcentration(sidx(scab_grid_m, idx));
        if (tnf_grid_m) sum_tnf += tnf_grid_m->GetConcentration(idx);
        if (il6_grid_m) sum_il6 += il6_grid_m->GetConcentration(idx);
        if (tcell_grid_m) sum_tcell += tcell_grid_m->GetConcentration(idx);
      }
      for (size_t idx : mask_.dermal_wound) {
        dermal_voxels++;
        sum_perf += vasc_grid->GetConcentration(idx);
        if (tgfb_grid) sum_tgfb += tgfb_grid->GetConcentration(idx);
        if (col_grid) sum_col += col_grid->GetConcentration(sidx(col_grid, idx));
        if (mmp_grid) sum_mmp += mmp_grid->GetConcentration(idx);
        if (fn_grid) sum_fn += fn_grid->GetConcentration(sidx(fn_grid, idx));
        if (elastin_grid) sum_elastin += elastin_grid->GetConcentration(sidx(elastin_grid, idx));
        if (ha_grid) sum_ha += ha_grid->GetConcentration(sidx(ha_grid, idx));
        if (temp_grid_m) sum_temp += temp_grid_m->GetConcentration(idx);
        if (glucose_grid_m) sum_gluc += glucose_grid_m->GetConcentration(idx);
        if (lactate_grid_m) sum_lac += lactate_grid_m->GetConcentration(idx);
        if (no_grid_m) sum_no += no_grid_m->GetConcentration(idx);
        if (ros_grid_m) sum_ros += ros_grid_m->GetConcentration(idx);
        if (stiff_grid_m) sum_stiff += stiff_grid_m->GetConcentration(idx);
        if (lymph_grid_m) sum_lymph += lymph_grid_m->GetConcentration(sidx(lymph_grid_m, idx));
        if (edema_grid_m) sum_edema += edema_grid_m->GetConcentration(idx);
        if (tnf_grid_m) sum_tnf += tnf_grid_m->GetConcentration(idx);
        if (il6_grid_m) sum_il6 += il6_grid_m->GetConcentration(idx);
        if (cart_grid_m) sum_cart += cart_grid_m->GetConcentration(sidx(cart_grid_m, idx));
        if (syn_grid_m) sum_syn += syn_grid_m->GetConcentration(sidx(syn_grid_m, idx));
        if (tcell_grid_m) sum_tcell += tcell_grid_m->GetConcentration(idx);
        if (bone_grid_m) sum_bone += bone_grid_m->GetConcentration(sidx(bone_grid_m, idx));
      }

      if (total_voxels > 0) {
        wound_closure = 100.0 * static_cast<real_t>(filled_voxels) /
                        static_cast<real_t>(total_voxels);
        mean_o2 = sum_o2 / total_voxels;
        mean_ca = sum_ca / total_voxels;
        mean_infl = sum_infl / total_voxels;
        scar_magnitude = sum_scar / total_voxels;
        mean_anti_infl = sum_anti / total_voxels;
        mean_biofilm = sum_biofilm / total_voxels;
        mean_vegf = sum_vegf / total_voxels;
        if (ph_grid) mean_ph = sum_ph / total_voxels;
        if (fibrin_grid) mean_fibrin = sum_fibrin / total_voxels;
        if (scab_grid_m) mean_scab = sum_scab / total_voxels;
        if (g_ecm_quality) mean_ecm_quality = sum_ecm / total_voxels;
        if (g_tissue_viability) mean_tissue_viability = sum_viability / total_voxels;
      }
      // ECM fields span both layers (fibroblasts deposit in dermis)
      int all_wound = total_voxels + dermal_voxels;
      if (all_wound > 0) {
        mean_tgfb = sum_tgfb / all_wound;
        mean_collagen = sum_col / all_wound;
        mean_mmp = sum_mmp / all_wound;
        mean_fn = sum_fn / all_wound;
        mean_elastin = sum_elastin / all_wound;
        mean_ha = sum_ha / all_wound;
        if (temp_grid_m) mean_temperature = sum_temp / all_wound;
        if (glucose_grid_m) mean_glucose = sum_gluc / all_wound;
        if (lactate_grid_m) mean_lactate = sum_lac / all_wound;
        if (no_grid_m) mean_no = sum_no / all_wound;
        if (ros_grid_m) mean_ros = sum_ros / all_wound;
        if (stiff_grid_m) mean_stiffness = sum_stiff / all_wound;
        if (lymph_grid_m) mean_lymphatic = sum_lymph / all_wound;
        if (edema_grid_m) mean_edema = sum_edema / all_wound;
        if (tnf_grid_m) mean_tnf_alpha = sum_tnf / all_wound;
        if (il6_grid_m) mean_il6 = sum_il6 / all_wound;
        if (cart_grid_m) mean_cartilage = sum_cart / all_wound;
        if (syn_grid_m) mean_synovial = sum_syn / all_wound;
        if (tcell_grid_m) mean_tcell = sum_tcell / all_wound;
        if (bone_grid_m) mean_bone = sum_bone / all_wound;
      }
      if (total_voxels > 0) {
        if (volt_grid_m) mean_voltage = sum_volt / total_voxels;
      }
      if (dermal_voxels > 0) {
        mean_perfusion = sum_perf / dermal_voxels;
      }

      // Per-sub-layer dermis means (dermis is structural)
      if (dermis_grid) {
        real_t s_pap = 0, s_ret = 0, s_hyp = 0;
        for (size_t idx : mask_.papillary_wound)
          s_pap += dermis_grid->GetConcentration(sidx(dermis_grid, idx));
        for (size_t idx : mask_.reticular_wound)
          s_ret += dermis_grid->GetConcentration(sidx(dermis_grid, idx));
        for (size_t idx : mask_.hypodermis_wound)
          s_hyp += dermis_grid->GetConcentration(sidx(dermis_grid, idx));
        if (!mask_.papillary_wound.empty())
          mean_dermis_pap = s_pap / mask_.papillary_wound.size();
        if (!mask_.reticular_wound.empty())
          mean_dermis_ret = s_ret / mask_.reticular_wound.size();
        if (!mask_.hypodermis_wound.empty())
          mean_dermis_hyp = s_hyp / mask_.hypodermis_wound.size();
      }
    }

    // --- Tumor field cells (binary field: 1 filled voxel ≈ 1 handoff) ---
    int tumor_field_cells = 0;
    if (sp->tumor.enabled) {
      auto* tumor_grid = rm->GetDiffusionGrid(fields::kTumorId);
      if (tumor_grid) {
        size_t n_boxes = tumor_grid->GetNumBoxes();
        for (size_t i = 0; i < n_boxes; i++) {
          if (tumor_grid->GetConcentration(i) > 0.5) tumor_field_cells++;
        }
      }
    }

    // --- Write CSV row ---
    real_t dt = sim->GetParam()->simulation_time_step;
    *ofs_ << step << ","
          << step * dt << ","
          << step * dt / 24.0 << ","
          << n_agents << ","
          << n_stem << ","
          << n_ta << ","
          << n_g0 << ","
          << n_basal << ","
          << n_spinous << ","
          << n_granular << ","
          << n_cornified << ","
          << n_neutrophils << ","
          << n_macrophages << ","
          << wound_closure << ","
          << mean_o2 << ","
          << mean_ca << ","
          << mean_infl << ","
          << scar_magnitude << ","
          << mean_anti_infl << ","
          << n_fibroblasts << ","
          << n_activated_fibroblasts << ","
          << n_myofibroblasts << ","
          << mean_tgfb << ","
          << mean_collagen << ","
          << mean_perfusion << ","
          << n_tumor_cells << ","
          << n_tumor_cycling << ","
          << tumor_field_cells << ","
          << mean_biofilm << ","
          << mean_vegf << ","
          << mean_mmp << ","
          << mean_fn << ","
          << mean_elastin << ","
          << mean_ha << ","
          << mean_dermis_pap << ","
          << mean_dermis_ret << ","
          << mean_dermis_hyp << ","
          << mean_ph << ","
          << mean_fibrin << ","
          << mean_ecm_quality << ","
          << mean_tissue_viability << ","
          << mean_temperature << ","
          << mean_glucose << ","
          << mean_lactate << ","
          << mean_no << ","
          << mean_ros << ","
          << mean_stiffness << ","
          << mean_lymphatic << ","
          << mean_edema << ","
          << mean_voltage << ","
          << mean_tnf_alpha << ","
          << mean_il6 << ","
          << mean_cartilage << ","
          << mean_synovial << ","
          << mean_tcell << ","
          << mean_bone << ","
          << mean_scab
          << std::endl;
    timer.Print("metrics");
  }

  void Close() {
    if (ofs_) {
      ofs_->close();
      ofs_.reset();
    }
  }

 private:
  std::shared_ptr<std::ofstream> ofs_;
  GridContext::WoundMaskData mask_;
  uint64_t last_step_ = 0;
  std::chrono::high_resolution_clock::time_point last_time_;
};

}  // namespace skibidy
}  // namespace bdm

#endif  // METRICS_H_
