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
#include "ros/ros_pde.h"
#include "mechanotransduction/stiffness_pde.h"
#include "lymphatic/lymphatic_pde.h"
#include "lymphatic/edema_pde.h"
#include "bioelectric/voltage_pde.h"
#include "rheumatoid/tnf_alpha_pde.h"
#include "rheumatoid/cartilage_pde.h"
#include "rheumatoid/il6_pde.h"
#include "rheumatoid/synovial_fluid_pde.h"
#include "rheumatoid/tcell_density_pde.h"
#include "rheumatoid/bone_pde.h"
#include "scab/scab_pde.h"
#include "photon/fluence_pde.h"
#include "photon/opsin_pde.h"

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
  if (sp->wound.enabled) {
    fields.Add(std::make_unique<CalciumPDE>());
    fields.Add(std::make_unique<KgfPDE>());
    fields.Add(std::make_unique<ScarPDE>());
    fields.Add(std::make_unique<WaterPDE>());
    if (sp->inflammation.split_inflammation_enabled) {
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
  if (sp->fibroblast.enabled) {
    fields.Add(std::make_unique<TGFBetaPDE>(sp));
    fields.Add(std::make_unique<CollagenPDE>(sp));
  }
  // Biofilm module: sessile biofilm density field
  if (sp->biofilm.enabled) {
    fields.Add(std::make_unique<BiofilmPDE>());
  }
  // Angiogenesis module: VEGF diffusion field
  if (sp->angiogenesis.enabled) {
    fields.Add(std::make_unique<VEGFPDE>(sp));
  }
  // MMP module: matrix metalloproteinase field + TIMP inhibitor field + pro-MMP zymogen
  if (sp->mmp.enabled) {
    fields.Add(std::make_unique<MMPPDE>(sp));
    fields.Add(std::make_unique<ProMMPPDE>(sp));
    fields.Add(std::make_unique<TIMPPDE>(sp));
  }
  // Fibronectin module: provisional wound matrix
  if (sp->fibronectin.enabled) {
    fields.Add(std::make_unique<FibronectinPDE>(sp));
  }
  // Elastin module: elastic fiber network
  if (sp->elastin.enabled) {
    fields.Add(std::make_unique<ElastinPDE>());
  }
  // Hyaluronan module: HA / GAG ground substance
  if (sp->hyaluronan.enabled) {
    fields.Add(std::make_unique<HyaluronanPDE>());
  }
  // Dermis module: dermal tissue integrity
  if (sp->dermis.enabled) {
    fields.Add(std::make_unique<DermisPDE>());
  }
  // pH module: wound alkalinity field
  fields.Add(std::make_unique<PHPDE>(sp));
  // Hemostasis module: fibrin clot provisional matrix
  if (sp->hemostasis.enabled) {
    fields.Add(std::make_unique<FibrinPDE>(sp));
  }
  // Tumor module: binary tumor tissue field
  if (sp->tumor.enabled) {
    fields.Add(std::make_unique<TumorPDE>());
  }
  // Temperature module: wound thermal regulation
  if (sp->temperature.enabled) {
    fields.Add(std::make_unique<TemperaturePDE>(sp));
  }
  // Glucose module: metabolic substrate
  if (sp->glucose_mod.enabled) {
    fields.Add(std::make_unique<GlucosePDE>(sp));
    // AGE field: accumulates in diabetic tissue from glucose glycation
    if (sp->diabetic.mode) {
      fields.Add(std::make_unique<AGEPDE>(sp));
    }
  }
  // Lactate module: hypoxia metabolite
  if (sp->lactate.enabled) {
    fields.Add(std::make_unique<LactatePDE>(sp));
  }
  // Nitric oxide module: immune antimicrobial and vasodilator
  if (sp->nitric_oxide.enabled) {
    fields.Add(std::make_unique<NitricOxidePDE>(sp));
  }
  // Basal density continuum: homeostatic keratinocyte population field.
  // D=0, decay=0 -- evolves via BasalDensityOp logistic PDE (not FTCS).
  if (sp->basal_density_enabled && sp->wound.enabled) {
    fields.Add(std::make_unique<BasalDensityPDE>());
  }
  // Senescence module: senescent cell density (non-diffusing accumulator)
  if (sp->senescence.enabled) {
    fields.Add(std::make_unique<SenescencePDE>(sp));
  }
  // Neuropathy module: nerve fiber density (slow diffusion = neurite extension)
  if (sp->neuropathy.enabled) {
    fields.Add(std::make_unique<NervePDE>(sp));
  }
  // ROS module: reactive oxygen species (oxidative stress mediator)
  if (sp->ros.enabled) {
    fields.Add(std::make_unique<ROSPDE>(sp));
  }
  // Mechanotransduction module: tissue stiffness (derived from ECM)
  if (sp->mechanotransduction.enabled) {
    fields.Add(std::make_unique<StiffnessPDE>(sp));
  }
  // Lymphatic module: vessel density + interstitial edema
  if (sp->lymphatic.enabled) {
    fields.Add(std::make_unique<LymphaticPDE>(sp));
    fields.Add(std::make_unique<EdemaPDE>(sp));
  }
  // Bioelectric module: transepithelial potential (galvanotaxis)
  if (sp->bioelectric.enabled) {
    fields.Add(std::make_unique<VoltagePDE>(sp));
  }
  // Scab module: protective wound crust (dried fibrin/platelet/exudate)
  if (sp->scab.enabled) {
    fields.Add(std::make_unique<ScabPDE>(sp));
  }
  // Photon transport module: fluence rate + opsin activation fields
  if (sp->photon.enabled) {
    fields.Add(std::make_unique<FluencePDE>(sp));
    fields.Add(std::make_unique<OpsinPDE>(sp));
  }
  // Rheumatoid arthritis module: TNF-alpha + IL-6 cytokines + cartilage + synovial fluid
  if (sp->ra.enabled) {
    fields.Add(std::make_unique<TNFAlphaPDE>(sp));
    fields.Add(std::make_unique<IL6PDE>(sp));
    fields.Add(std::make_unique<CartilagePDE>(sp));
    fields.Add(std::make_unique<SynovialFluidPDE>(sp));
    fields.Add(std::make_unique<TCellDensityPDE>(sp));
    fields.Add(std::make_unique<BonePDE>(sp));
  }
  fields.InitAll(sim);

  // --- Performance: skip FTCS solver for decay-only fields (D=0, decay>0) ---
  auto* rm = sim->GetResourceManager();
  if (sp->fibronectin.enabled && sp->fibronectin.decay > 0) {
    rm->GetDiffusionGrid(fields::kFibronectinId)->SetTimeStep(1e30);
  }
  if (sp->elastin.enabled && sp->elastin.decay > 0) {
    rm->GetDiffusionGrid(fields::kElastinId)->SetTimeStep(1e30);
  }
  if (sp->hemostasis.enabled) {
    rm->GetDiffusionGrid(fields::kFibrinId)->SetTimeStep(1e30);
  }
  if (sp->scab.enabled) {
    rm->GetDiffusionGrid(fields::kScabId)->SetTimeStep(1e30);
  }
  // BasalDensity: D=0, decay=0 -- evolved by BasalDensityOp (logistic, not FTCS).
  // Disable BDM's FTCS step entirely so the grid is only touched by the op.
  if (sp->basal_density_enabled && sp->wound.enabled) {
    rm->GetDiffusionGrid(fields::kBasalDensityId)->SetTimeStep(1e30);
  }
  // AGE: D=0, decay~0 -- structural accumulator, skip FTCS.
  if (sp->glucose_mod.enabled && sp->diabetic.mode) {
    rm->GetDiffusionGrid(fields::kAGEId)->SetTimeStep(1e30);
  }
  // Senescence: D=0, non-diffusing -- evolved in fused_post (accumulation + SASP).
  if (sp->senescence.enabled) {
    rm->GetDiffusionGrid(fields::kSenescenceId)->SetTimeStep(1e30);
  }
  // Stiffness: D=0, decay=0 -- derived structural property, computed in fused_source.
  if (sp->mechanotransduction.enabled) {
    rm->GetDiffusionGrid(fields::kStiffnessId)->SetTimeStep(1e30);
  }
  // Edema: D=0, decay=0 -- evolved mechanistically in fused_source (leak vs drain).
  if (sp->lymphatic.enabled) {
    rm->GetDiffusionGrid(fields::kEdemaId)->SetTimeStep(1e30);
  }
  // Opsin: D=0 (membrane-bound), evolved by PhotonSourceHook kinetics.
  if (sp->photon.enabled) {
    rm->GetDiffusionGrid(fields::kOpsinId)->SetTimeStep(1e30);
  }
  // Cartilage + Synovial fluid: D=0, decay=0 -- structural fields, evolved mechanistically.
  if (sp->ra.enabled) {
    rm->GetDiffusionGrid(fields::kCartilageId)->SetTimeStep(1e30);
    rm->GetDiffusionGrid(fields::kSynovialFluidId)->SetTimeStep(1e30);
    rm->GetDiffusionGrid(fields::kBoneId)->SetTimeStep(1e30);
  }

  // --- PDE sub-cycling: solve slow fields less frequently ---
  real_t dt = sim->GetParam()->simulation_time_step;
  if (sp->subcycle_slow > 1) {
    real_t slow_dt = sp->subcycle_slow * dt;
    if (sp->wound.enabled) {
      rm->GetDiffusionGrid(fields::kWaterId)->SetTimeStep(slow_dt);
    }
    if (sp->hyaluronan.enabled) {
      rm->GetDiffusionGrid(fields::kHyaluronanId)->SetTimeStep(slow_dt);
    }
    if (sp->perfusion.diffusion > 0) {
      rm->GetDiffusionGrid(fields::kVascularId)->SetTimeStep(slow_dt);
    }
    if (sp->dermis.enabled && sp->dermis.diffusion > 0) {
      rm->GetDiffusionGrid(fields::kDermisId)->SetTimeStep(slow_dt);
    }
    if (sp->temperature.enabled) {
      rm->GetDiffusionGrid(fields::kTemperatureId)->SetTimeStep(slow_dt);
    }
    if (sp->glucose_mod.enabled) {
      rm->GetDiffusionGrid(fields::kGlucoseId)->SetTimeStep(slow_dt);
    }
    if (sp->neuropathy.enabled && sp->neuropathy.diffusion > 0) {
      rm->GetDiffusionGrid(fields::kNerveId)->SetTimeStep(slow_dt);
    }
    if (sp->lymphatic.enabled && sp->lymphatic.diffusion > 0) {
      rm->GetDiffusionGrid(fields::kLymphaticId)->SetTimeStep(slow_dt);
    }
  }
  if (sp->subcycle_medium > 1) {
    real_t med_dt = sp->subcycle_medium * dt;
    if (sp->wound.enabled) {
      if (sp->inflammation.split_inflammation_enabled) {
        rm->GetDiffusionGrid(fields::kProInflammatoryId)->SetTimeStep(med_dt);
        rm->GetDiffusionGrid(fields::kAntiInflammatoryId)->SetTimeStep(med_dt);
      } else {
        rm->GetDiffusionGrid(fields::kInflammationId)->SetTimeStep(med_dt);
      }
      rm->GetDiffusionGrid(fields::kImmunePressureId)->SetTimeStep(med_dt);
    }
    if (sp->fibroblast.enabled) {
      rm->GetDiffusionGrid(fields::kTGFBetaId)->SetTimeStep(med_dt);
    }
    if (sp->mmp.enabled) {
      rm->GetDiffusionGrid(fields::kMMPId)->SetTimeStep(med_dt);
      rm->GetDiffusionGrid(fields::kProMMPId)->SetTimeStep(med_dt);
      rm->GetDiffusionGrid(fields::kTIMPId)->SetTimeStep(med_dt);
    }
    if (sp->angiogenesis.enabled) {
      rm->GetDiffusionGrid(fields::kVEGFId)->SetTimeStep(med_dt);
    }
    if (sp->lactate.enabled) {
      rm->GetDiffusionGrid(fields::kLactateId)->SetTimeStep(med_dt);
    }
    if (sp->nitric_oxide.enabled) {
      rm->GetDiffusionGrid(fields::kNitricOxideId)->SetTimeStep(med_dt);
    }
    if (sp->ros.enabled) {
      if (auto* g = rm->GetDiffusionGrid(fields::kROSId)) g->SetTimeStep(med_dt);
    }
    if (sp->ra.enabled) {
      if (auto* g = rm->GetDiffusionGrid(fields::kTNFAlphaId)) g->SetTimeStep(med_dt);
      if (auto* g = rm->GetDiffusionGrid(fields::kIL6Id)) g->SetTimeStep(med_dt);
      if (auto* g = rm->GetDiffusionGrid(fields::kTCellDensityId)) g->SetTimeStep(med_dt);
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
  if (sp->fibroblast.enabled) {
    OperationRegistry::GetInstance()->AddOperationImpl(
        "FibroblastRecruitment", OpComputeTarget::kCpu,
        new FibroblastRecruitment());
    auto* fibro_op = NewOperation("FibroblastRecruitment");
    scheduler->ScheduleOp(fibro_op, OpType::kPreSchedule);
  }

  // --- Tumor initiation (seeds tumor cluster at configured step) ---
  if (sp->tumor.enabled) {
    OperationRegistry::GetInstance()->AddOperationImpl(
        "TumorInitiation", OpComputeTarget::kCpu, new TumorInitiation());
    auto* tumor_op = NewOperation("TumorInitiation");
    scheduler->ScheduleOp(tumor_op, OpType::kPreSchedule);
  }

  // --- PDE source terms: fused single-pass when wound_enabled ---
  if (sp->wound.enabled) {
    OperationRegistry::GetInstance()->AddOperationImpl(
        "FusedWoundSourceOp", OpComputeTarget::kCpu,
        new FusedWoundSourceOp());
    auto* fused_src = NewOperation("FusedWoundSourceOp");
    scheduler->ScheduleOp(fused_src, OpType::kPreSchedule);

    bool need_post = sp->biofilm.enabled || sp->scar.proportional_enabled ||
                     sp->diabetic.mode || sp->mmp.enabled ||
                     sp->fibronectin.enabled || sp->hemostasis.enabled ||
                     sp->senescence.enabled || sp->neuropathy.enabled ||
                     sp->ros.enabled || sp->mechanotransduction.enabled ||
                     sp->lymphatic.enabled || sp->ra.enabled;
    if (need_post) {
      OperationRegistry::GetInstance()->AddOperationImpl(
          "FusedWoundPostOp", OpComputeTarget::kCpu, new FusedWoundPostOp());
      auto* fused_post = NewOperation("FusedWoundPostOp");
      scheduler->ScheduleOp(fused_post);
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
  if (sp->basal_density_enabled && sp->wound.enabled) {
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
