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
// 40 columns -- see header string below for full list
// ---------------------------------------------------------------------------
struct MetricsExporter : public StandaloneOperationImpl {
  BDM_OP_HEADER(MetricsExporter);

  void operator()() override {
    auto* sim = Simulation::GetActive();
    auto* scheduler = sim->GetScheduler();
    auto* sp = sim->GetParam()->Get<SimParam>();
    PerfTimer timer(sp->debug_perf);

    uint64_t step = scheduler->GetSimulatedSteps();

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
            << "mean_ecm_quality,mean_tissue_viability"
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
            sp->tumor_handoff_delay > 0 &&
            tumor->GetG0Steps() > sp->tumor_handoff_delay) {
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

    if (sp->wound_enabled) {
      auto* stratum_grid = rm->GetDiffusionGrid(fields::kStratum);
      auto* o2_grid = rm->GetDiffusionGrid(fields::kOxygen);
      auto* ca_grid = rm->GetDiffusionGrid(fields::kCalcium);
      auto* scar_grid = rm->GetDiffusionGrid(fields::kScar);
      auto* vasc_grid = rm->GetDiffusionGrid(fields::kVascular);

      DiffusionGrid* tgfb_grid = nullptr;
      DiffusionGrid* col_grid = nullptr;
      if (sp->fibroblast_enabled) {
        tgfb_grid = rm->GetDiffusionGrid(fields::kTGFBeta);
        col_grid = rm->GetDiffusionGrid(fields::kCollagen);
      }

      DiffusionGrid* biofilm_grid = nullptr;
      if (sp->biofilm_enabled) {
        biofilm_grid = rm->GetDiffusionGrid(fields::kBiofilm);
      }
      DiffusionGrid* vegf_grid = nullptr;
      if (sp->angiogenesis_enabled) {
        vegf_grid = rm->GetDiffusionGrid(fields::kVEGF);
      }
      DiffusionGrid* mmp_grid = nullptr;
      if (sp->mmp_enabled) {
        mmp_grid = rm->GetDiffusionGrid(fields::kMMP);
      }
      DiffusionGrid* fn_grid = nullptr;
      if (sp->fibronectin_enabled) {
        fn_grid = rm->GetDiffusionGrid(fields::kFibronectin);
      }
      DiffusionGrid* elastin_grid = nullptr;
      if (sp->elastin_enabled) {
        elastin_grid = rm->GetDiffusionGrid(fields::kElastin);
      }
      DiffusionGrid* ha_grid = nullptr;
      if (sp->hyaluronan_enabled) {
        ha_grid = rm->GetDiffusionGrid(fields::kHyaluronan);
      }
      DiffusionGrid* dermis_grid = nullptr;
      if (sp->dermis_enabled) {
        dermis_grid = rm->GetDiffusionGrid(fields::kDermis);
      }
      DiffusionGrid* ph_grid = rm->GetDiffusionGrid(fields::kPH);
      DiffusionGrid* fibrin_grid = nullptr;
      if (sp->hemostasis_enabled) {
        fibrin_grid = rm->GetDiffusionGrid(fields::kFibrin);
      }

      DiffusionGrid* infl_grid = nullptr;
      DiffusionGrid* anti_grid = nullptr;
      if (sp->split_inflammation_enabled) {
        infl_grid = rm->GetDiffusionGrid(fields::kProInflammatory);
        anti_grid = rm->GetDiffusionGrid(fields::kAntiInflammatory);
      } else {
        infl_grid = rm->GetDiffusionGrid(fields::kInflammation);
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

      // Compact loops: iterate only wound voxels via precomputed index lists.
      for (size_t idx : mask_.epi_wound) {
        total_voxels++;
        if (std::abs(stratum_grid->GetConcentration(idx)) > 0.5) {
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
    if (sp->tumor_enabled) {
      auto* tumor_grid = rm->GetDiffusionGrid(fields::kTumor);
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
          << mean_tissue_viability
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
