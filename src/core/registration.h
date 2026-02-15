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
#include "mmp/prommp_pde.h"
#include "mmp/timp_pde.h"
#include "fibronectin/fibronectin_pde.h"
#include "elastin/elastin_pde.h"
#include "hyaluronan/hyaluronan_pde.h"
#include "dermis/dermis_pde.h"
#include "ph/ph_pde.h"
#include "hemostasis/hemostasis_pde.h"
#include "tumor/tumor_pde.h"
#include "temperature/temperature_pde.h"
#include "glucose/glucose_pde.h"
#include "glucose/age_pde.h"
#include "lactate/lactate_pde.h"
#include "nitric_oxide/nitric_oxide_pde.h"
#include "tissue/basal_density_pde.h"
#include "senescence/senescence_pde.h"
#include "neuropathy/nerve_pde.h"

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
#include "core/voxel_env.h"

#include "core/batch_diffusion.h"
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
  // MMP module: matrix metalloproteinase field + TIMP inhibitor field + pro-MMP zymogen
  if (sp->mmp_enabled) {
    fields.Add(std::make_unique<MMPPDE>(sp));
    fields.Add(std::make_unique<ProMMPPDE>(sp));
    fields.Add(std::make_unique<TIMPPDE>(sp));
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
  // Temperature module: wound thermal regulation
  if (sp->temperature_enabled) {
    fields.Add(std::make_unique<TemperaturePDE>(sp));
  }
  // Glucose module: metabolic substrate
  if (sp->glucose_enabled) {
    fields.Add(std::make_unique<GlucosePDE>(sp));
    // AGE field: accumulates in diabetic tissue from glucose glycation
    if (sp->diabetic_mode) {
      fields.Add(std::make_unique<AGEPDE>(sp));
    }
  }
  // Lactate module: hypoxia metabolite
  if (sp->lactate_enabled) {
    fields.Add(std::make_unique<LactatePDE>(sp));
  }
  // Nitric oxide module: immune antimicrobial and vasodilator
  if (sp->nitric_oxide_enabled) {
    fields.Add(std::make_unique<NitricOxidePDE>(sp));
  }
  // Basal density continuum: homeostatic keratinocyte population field.
  // D=0, decay=0 -- evolves via BasalDensityOp logistic PDE (not FTCS).
  if (sp->basal_density_enabled && sp->wound_enabled) {
    fields.Add(std::make_unique<BasalDensityPDE>());
  }
  // Senescence module: senescent cell density (non-diffusing accumulator)
  if (sp->senescence_enabled) {
    fields.Add(std::make_unique<SenescencePDE>(sp));
  }
  // Neuropathy module: nerve fiber density (slow diffusion = neurite extension)
  if (sp->neuropathy_enabled) {
    fields.Add(std::make_unique<NervePDE>(sp));
  }
  fields.InitAll(sim);

  // --- Performance: skip FTCS solver for decay-only fields (D=0, decay>0) ---
  auto* rm = sim->GetResourceManager();
  if (sp->fibronectin_enabled && sp->fibronectin_decay > 0) {
    rm->GetDiffusionGrid(fields::kFibronectinId)->SetTimeStep(1e30);
  }
  if (sp->elastin_enabled && sp->elastin_decay > 0) {
    rm->GetDiffusionGrid(fields::kElastinId)->SetTimeStep(1e30);
  }
  if (sp->hemostasis_enabled) {
    rm->GetDiffusionGrid(fields::kFibrinId)->SetTimeStep(1e30);
  }
  // BasalDensity: D=0, decay=0 -- evolved by BasalDensityOp (logistic, not FTCS).
  // Disable BDM's FTCS step entirely so the grid is only touched by the op.
  if (sp->basal_density_enabled && sp->wound_enabled) {
    rm->GetDiffusionGrid(fields::kBasalDensityId)->SetTimeStep(1e30);
  }
  // AGE: D=0, decay~0 -- structural accumulator, skip FTCS.
  if (sp->glucose_enabled && sp->diabetic_mode) {
    rm->GetDiffusionGrid(fields::kAGEId)->SetTimeStep(1e30);
  }
  // Senescence: D=0, non-diffusing -- evolved in fused_post (accumulation + SASP).
  if (sp->senescence_enabled) {
    rm->GetDiffusionGrid(fields::kSenescenceId)->SetTimeStep(1e30);
  }

  // --- PDE sub-cycling: solve slow fields less frequently ---
  real_t dt = sim->GetParam()->simulation_time_step;
  if (sp->subcycle_slow > 1) {
    real_t slow_dt = sp->subcycle_slow * dt;
    if (sp->wound_enabled) {
      rm->GetDiffusionGrid(fields::kWaterId)->SetTimeStep(slow_dt);
    }
    if (sp->hyaluronan_enabled) {
      rm->GetDiffusionGrid(fields::kHyaluronanId)->SetTimeStep(slow_dt);
    }
    if (sp->perfusion_diffusion > 0) {
      rm->GetDiffusionGrid(fields::kVascularId)->SetTimeStep(slow_dt);
    }
    if (sp->dermis_enabled && sp->dermis_diffusion > 0) {
      rm->GetDiffusionGrid(fields::kDermisId)->SetTimeStep(slow_dt);
    }
    if (sp->temperature_enabled) {
      rm->GetDiffusionGrid(fields::kTemperatureId)->SetTimeStep(slow_dt);
    }
    if (sp->glucose_enabled) {
      rm->GetDiffusionGrid(fields::kGlucoseId)->SetTimeStep(slow_dt);
    }
    if (sp->neuropathy_enabled && sp->neuropathy_diffusion > 0) {
      rm->GetDiffusionGrid(fields::kNerveId)->SetTimeStep(slow_dt);
    }
  }
  if (sp->subcycle_medium > 1) {
    real_t med_dt = sp->subcycle_medium * dt;
    if (sp->wound_enabled) {
      if (sp->split_inflammation_enabled) {
        rm->GetDiffusionGrid(fields::kProInflammatoryId)->SetTimeStep(med_dt);
        rm->GetDiffusionGrid(fields::kAntiInflammatoryId)->SetTimeStep(med_dt);
      } else {
        rm->GetDiffusionGrid(fields::kInflammationId)->SetTimeStep(med_dt);
      }
      rm->GetDiffusionGrid(fields::kImmunePressureId)->SetTimeStep(med_dt);
    }
    if (sp->fibroblast_enabled) {
      rm->GetDiffusionGrid(fields::kTGFBetaId)->SetTimeStep(med_dt);
    }
    if (sp->mmp_enabled) {
      rm->GetDiffusionGrid(fields::kMMPId)->SetTimeStep(med_dt);
      rm->GetDiffusionGrid(fields::kProMMPId)->SetTimeStep(med_dt);
      rm->GetDiffusionGrid(fields::kTIMPId)->SetTimeStep(med_dt);
    }
    if (sp->angiogenesis_enabled) {
      rm->GetDiffusionGrid(fields::kVEGFId)->SetTimeStep(med_dt);
    }
    if (sp->lactate_enabled) {
      rm->GetDiffusionGrid(fields::kLactateId)->SetTimeStep(med_dt);
    }
    if (sp->nitric_oxide_enabled) {
      rm->GetDiffusionGrid(fields::kNitricOxideId)->SetTimeStep(med_dt);
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
                     sp->fibronectin_enabled || sp->hemostasis_enabled ||
                     sp->senescence_enabled || sp->neuropathy_enabled;
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

  // --- Fused per-voxel environment snapshot (before behavior dispatch) ---
  OperationRegistry::GetInstance()->AddOperationImpl(
      "VoxelEnvFillOp", OpComputeTarget::kCpu, new VoxelEnvFillOp());
  auto* venv_op = NewOperation("VoxelEnvFillOp");
  scheduler->ScheduleOp(venv_op, OpType::kPreSchedule);

  // --- Basal density continuum evolution ---
  if (sp->basal_density_enabled && sp->wound_enabled) {
    OperationRegistry::GetInstance()->AddOperationImpl(
        "BasalDensityOp", OpComputeTarget::kCpu, new BasalDensityOp());
    auto* density_op = NewOperation("BasalDensityOp");
    scheduler->ScheduleOp(density_op, OpType::kPostSchedule);
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

  // --- Batch diffusion: replace BDM's ContinuumOp with lean scheduler ---
  // Must be registered AFTER sub-cycling SetTimeStep calls in RegisterFields
  // so that CategorizeGrids() can read the configured time steps.
  OperationRegistry::GetInstance()->AddOperationImpl(
      "BatchDiffusionOp", OpComputeTarget::kCpu, new BatchDiffusionOp());
  auto* batch_diff = NewOperation("BatchDiffusionOp");
  scheduler->ScheduleOp(batch_diff, OpType::kSchedule);

  // Unschedule BDM's default ContinuumOp ("continuum" is not protected).
  auto continuum_ops = scheduler->GetOps("continuum");
  if (!continuum_ops.empty()) {
    scheduler->UnscheduleOp(continuum_ops[0]);
  }
}

}  // namespace skibidy
}  // namespace bdm

#endif  // REGISTRATION_H_
