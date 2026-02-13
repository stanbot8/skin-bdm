#ifndef REGISTRATION_H_
#define REGISTRATION_H_

// Field and operation registration for skibidy.
//
// Extracted from skibidy.h so that adding a new module means editing this
// file (field registration + operation scheduling) rather than the main
// simulation orchestrator.

#include "core/composite_field.h"
#include "core/field_names.h"
#include "infra/sim_param.h"
#include "infra/volume.h"

// PDE headers
#include "perfusion/vascular.h"
#include "tissue/oxygen.h"
#include "tissue/stratum.h"
#include "tissue/calcium.h"
#include "tissue/kgf.h"
#include "tissue/water.h"
#include "inflammation/inflammation_pde.h"
#include "scar/scar_pde.h"
#include "fibroblast/tgfbeta_pde.h"
#include "fibroblast/collagen_pde.h"
#include "biofilm/biofilm_pde.h"
#include "angiogenesis/vegf_pde.h"
#include "mmp/mmp_pde.h"
#include "fibronectin/fibronectin_pde.h"
#include "elastin/elastin_pde.h"
#include "hyaluronan/hyaluronan_pde.h"
#include "dermis/dermis_pde.h"
#include "ph/ph_pde.h"
#include "hemostasis/hemostasis_pde.h"
#include "tumor/tumor_pde.h"

// Operation headers
#include "wound/wound_event.h"
#include "immune/immune_response.h"
#include "fibroblast/fibroblast_recruitment.h"
#include "tumor/tumor_initiation.h"
#include "diabetic/baseline_inflammation.h"
#include "biofilm/biofilm_op.h"
#include "angiogenesis/vegf_op.h"
#include "core/derived_field.h"
#include "core/derived_fields_op.h"
#include "core/fused_post.h"
#include "core/fused_source.h"
#include "core/hot_reload.h"
#include "core/metrics.h"

#include "core/operation/operation_registry.h"

namespace bdm {
namespace skibidy {

// Register all PDE channels into the CompositeField, initialize grids,
// and configure sub-cycling for performance.
inline void RegisterFields(Simulation* sim, const SimParam* sp,
                           CompositeField& fields) {
  // Core fields always needed (VascularPDE before OxygenPDE so O2 reads perfusion)
  fields.Add(std::make_unique<VascularPDE>());
  fields.Add(std::make_unique<OxygenPDE>());
  fields.Add(std::make_unique<StratumPDE>());
  // Wound-related fields only when wound is enabled
  if (sp->wound_enabled) {
    fields.Add(std::make_unique<CalciumPDE>());
    fields.Add(std::make_unique<KgfPDE>());
    fields.Add(std::make_unique<ScarPDE>());
    fields.Add(std::make_unique<WaterPDE>());
    if (sp->split_inflammation_enabled) {
      fields.Add(std::make_unique<ProInflammatoryPDE>(sp));
      fields.Add(std::make_unique<AntiInflammatoryPDE>(sp));
    } else {
      fields.Add(std::make_unique<InflammationPDE>(sp));
    }
    // Immune pressure: same dynamics as inflammation but excludes wound DAMPs.
    // Keratinocytes read this for suppression instead of the inflammation field.
    fields.Add(std::make_unique<ImmunePressurePDE>(sp));
  }
  // Fibroblast module: TGF-beta and Collagen fields
  if (sp->fibroblast_enabled) {
    fields.Add(std::make_unique<TGFBetaPDE>(sp));
    fields.Add(std::make_unique<CollagenPDE>(sp));
  }
  // Biofilm module: sessile biofilm density field
  if (sp->biofilm_enabled) {
    fields.Add(std::make_unique<BiofilmPDE>());
  }
  // Angiogenesis module: VEGF diffusion field
  if (sp->angiogenesis_enabled) {
    fields.Add(std::make_unique<VEGFPDE>(sp));
  }
  // MMP module: matrix metalloproteinase field
  if (sp->mmp_enabled) {
    fields.Add(std::make_unique<MMPPDE>(sp));
  }
  // Fibronectin module: provisional wound matrix
  if (sp->fibronectin_enabled) {
    fields.Add(std::make_unique<FibronectinPDE>(sp));
  }
  // Elastin module: elastic fiber network
  if (sp->elastin_enabled) {
    fields.Add(std::make_unique<ElastinPDE>());
  }
  // Hyaluronan module: HA / GAG ground substance
  if (sp->hyaluronan_enabled) {
    fields.Add(std::make_unique<HyaluronanPDE>());
  }
  // Dermis module: dermal tissue integrity
  if (sp->dermis_enabled) {
    fields.Add(std::make_unique<DermisPDE>());
  }
  // pH module: wound alkalinity field
  fields.Add(std::make_unique<PHPDE>(sp));
  // Hemostasis module: fibrin clot provisional matrix
  if (sp->hemostasis_enabled) {
    fields.Add(std::make_unique<FibrinPDE>(sp));
  }
  // Tumor module: binary tumor tissue field
  if (sp->tumor_enabled) {
    fields.Add(std::make_unique<TumorPDE>());
  }
  fields.InitAll(sim);

  // --- Performance: skip FTCS solver for decay-only fields (D=0, decay>0) ---
  auto* rm = sim->GetResourceManager();
  if (sp->fibronectin_enabled && sp->fibronectin_decay > 0) {
    rm->GetDiffusionGrid(fields::kFibronectin)->SetTimeStep(1e30);
  }
  if (sp->elastin_enabled && sp->elastin_decay > 0) {
    rm->GetDiffusionGrid(fields::kElastin)->SetTimeStep(1e30);
  }
  if (sp->hemostasis_enabled) {
    rm->GetDiffusionGrid(fields::kFibrin)->SetTimeStep(1e30);
  }

  // --- PDE sub-cycling: solve slow fields less frequently ---
  real_t dt = sim->GetParam()->simulation_time_step;
  if (sp->subcycle_slow > 1) {
    real_t slow_dt = sp->subcycle_slow * dt;
    if (sp->wound_enabled) {
      rm->GetDiffusionGrid(fields::kWater)->SetTimeStep(slow_dt);
    }
    if (sp->hyaluronan_enabled) {
      rm->GetDiffusionGrid(fields::kHyaluronan)->SetTimeStep(slow_dt);
    }
    if (sp->perfusion_diffusion > 0) {
      rm->GetDiffusionGrid(fields::kVascular)->SetTimeStep(slow_dt);
    }
    if (sp->dermis_enabled && sp->dermis_diffusion > 0) {
      rm->GetDiffusionGrid(fields::kDermis)->SetTimeStep(slow_dt);
    }
  }
  if (sp->subcycle_medium > 1) {
    real_t med_dt = sp->subcycle_medium * dt;
    if (sp->wound_enabled) {
      if (sp->split_inflammation_enabled) {
        rm->GetDiffusionGrid(fields::kProInflammatory)->SetTimeStep(med_dt);
        rm->GetDiffusionGrid(fields::kAntiInflammatory)->SetTimeStep(med_dt);
      } else {
        rm->GetDiffusionGrid(fields::kInflammation)->SetTimeStep(med_dt);
      }
      rm->GetDiffusionGrid(fields::kImmunePressure)->SetTimeStep(med_dt);
    }
    if (sp->fibroblast_enabled) {
      rm->GetDiffusionGrid(fields::kTGFBeta)->SetTimeStep(med_dt);
    }
    if (sp->mmp_enabled) {
      rm->GetDiffusionGrid(fields::kMMP)->SetTimeStep(med_dt);
    }
    if (sp->angiogenesis_enabled) {
      rm->GetDiffusionGrid(fields::kVEGF)->SetTimeStep(med_dt);
    }
  }
}

// Register all scheduled operations (wound events, immune response,
// fused source/post ops, metrics, hot-reload).
inline void RegisterOperations(Simulation* sim, const SimParam* sp,
                               CompositeField* fields,
                               DerivedField* ecm_quality = nullptr,
                               DerivedField* tissue_viability = nullptr,
                               DerivedField* wound_microenv = nullptr) {
  auto* scheduler = sim->GetScheduler();

  // --- Wound event (fires once at configured step) ---
  OperationRegistry::GetInstance()->AddOperationImpl(
      "WoundEvent", OpComputeTarget::kCpu, new WoundEvent(fields));
  auto* wound_op = NewOperation("WoundEvent");
  scheduler->ScheduleOp(wound_op, OpType::kPreSchedule);

  // --- Wound resolution (dissolves agents when healed) ---
  OperationRegistry::GetInstance()->AddOperationImpl(
      "WoundResolution", OpComputeTarget::kCpu, new WoundResolution());
  auto* resolve_op = NewOperation("WoundResolution");
  scheduler->ScheduleOp(resolve_op, OpType::kPostSchedule);

  // --- Immune response (delayed immune cell spawning) ---
  OperationRegistry::GetInstance()->AddOperationImpl(
      "ImmuneResponse", OpComputeTarget::kCpu, new ImmuneResponse());
  auto* immune_op = NewOperation("ImmuneResponse");
  scheduler->ScheduleOp(immune_op, OpType::kPreSchedule);

  // --- Fibroblast recruitment (delayed spawning at wound margin) ---
  if (sp->fibroblast_enabled) {
    OperationRegistry::GetInstance()->AddOperationImpl(
        "FibroblastRecruitment", OpComputeTarget::kCpu,
        new FibroblastRecruitment());
    auto* fibro_op = NewOperation("FibroblastRecruitment");
    scheduler->ScheduleOp(fibro_op, OpType::kPreSchedule);
  }

  // --- Tumor initiation (seeds tumor cluster at configured step) ---
  if (sp->tumor_enabled) {
    OperationRegistry::GetInstance()->AddOperationImpl(
        "TumorInitiation", OpComputeTarget::kCpu, new TumorInitiation());
    auto* tumor_op = NewOperation("TumorInitiation");
    scheduler->ScheduleOp(tumor_op, OpType::kPreSchedule);
  }

  // --- PDE source terms: fused single-pass when wound_enabled ---
  if (sp->wound_enabled) {
    OperationRegistry::GetInstance()->AddOperationImpl(
        "FusedWoundSourceOp", OpComputeTarget::kCpu,
        new FusedWoundSourceOp());
    auto* fused_src = NewOperation("FusedWoundSourceOp");
    scheduler->ScheduleOp(fused_src, OpType::kPreSchedule);

    bool need_post = sp->biofilm_enabled || sp->scar_proportional_enabled ||
                     sp->diabetic_mode || sp->mmp_enabled ||
                     sp->fibronectin_enabled || sp->hemostasis_enabled;
    if (need_post) {
      OperationRegistry::GetInstance()->AddOperationImpl(
          "FusedWoundPostOp", OpComputeTarget::kCpu, new FusedWoundPostOp());
      auto* fused_post = NewOperation("FusedWoundPostOp");
      scheduler->ScheduleOp(fused_post, OpType::kPostSchedule);
    }
  } else {
    OperationRegistry::GetInstance()->AddOperationImpl(
        "CompositeFieldOp", OpComputeTarget::kCpu,
        new CompositeFieldOp(fields));
    auto* field_op = NewOperation("CompositeFieldOp");
    scheduler->ScheduleOp(field_op, OpType::kPreSchedule);
  }

  // --- Derived composite fields (computed once per step) ---
  if (ecm_quality && tissue_viability && wound_microenv) {
    OperationRegistry::GetInstance()->AddOperationImpl(
        "DerivedFieldsOp", OpComputeTarget::kCpu,
        new DerivedFieldsOp(ecm_quality, tissue_viability, wound_microenv));
    auto* derived_op = NewOperation("DerivedFieldsOp");
    scheduler->ScheduleOp(derived_op, OpType::kPostSchedule);
  }

  // --- Metrics CSV export ---
  OperationRegistry::GetInstance()->AddOperationImpl(
      "MetricsExporter", OpComputeTarget::kCpu, new MetricsExporter());
  auto* metrics_op = NewOperation("MetricsExporter");
  scheduler->ScheduleOp(metrics_op, OpType::kPostSchedule);

  // --- Hot-reload: watch bdm.toml for runtime parameter changes ---
  if (sp->hot_reload) {
    OperationRegistry::GetInstance()->AddOperationImpl(
        "HotReloadOp", OpComputeTarget::kCpu, new HotReloadOp());
    auto* reload_op = NewOperation("HotReloadOp");
    scheduler->ScheduleOp(reload_op, OpType::kPreSchedule);
    std::cout << "[hot-reload] enabled -- edit bdm.toml while sim runs"
              << std::endl;
  }
}

}  // namespace skibidy
}  // namespace bdm

#endif  // REGISTRATION_H_
