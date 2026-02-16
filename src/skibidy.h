#ifndef SKIBIDY_H_
#define SKIBIDY_H_

#include <cmath>
#include <cstdlib>
#include <memory>

#include "tissue/basal_division.h"
#include "core/cell_cell_force.h"
#include "tissue/differentiation.h"
#include "tissue/keratinocyte.h"
#include "tissue/shedding.h"
#include "infra/sim_param.h"
#include "infra/util.h"
#include "core/registration.h"
#include "core/metrics.h"
#include "tumor/tumor_cell.h"

#include "core/environment/uniform_grid_environment.h"
#include "core/operation/mechanical_forces_op.h"

namespace bdm {
namespace skibidy {

inline int Simulate(int argc, const char** argv) {
  Param::RegisterParamGroup(new SimParam());

  auto set_param = [](Param* param) {
    param->bound_space = Param::BoundSpaceMode::kClosed;
    // Fixed domain locks DiffusionGrid size independent of agent presence.
    // Padded beyond tissue [0,30] for boundary headroom.
    param->min_bound = -10;
    param->max_bound = 40;
    param->simulation_time_step = 0.1;  // 1 step = 0.1 hours (6 minutes)
    // Gradients are computed on-the-fly where needed (chemotaxis only uses
    // 2 of ~16 fields). Precomputing all gradients every step is pure waste.
    param->calculate_gradients = false;
    // Per-axis bounds for position clamping (match grid range)
    auto* sp = const_cast<SimParam*>(param->Get<SimParam>());
    sp->bounds_min = {0, 0, -10};     // x/y: tissue range, z: dermis
    sp->bounds_max = {30, 30, 40};    // x/y: tissue range, z: above cornified
  };

  Simulation simulation(argc, argv, set_param);

  // Load skin parameters from bdm.toml using bundled toml++.
  // Bypasses BDM's config chain for v1.04 compatibility.
  {
    auto toml_config = toml::parse_file("bdm.toml");
    auto* sp_load = const_cast<SimParam*>(
        simulation.GetParam()->Get<SimParam>());
    sp_load->LoadConfig(toml_config);
    sp_load->ValidateConfig();
  }

  // Lock grid dimensions to param bounds (min_bound/max_bound).
  // Without this, BioDynaMo auto-sizes the grid to the agent bounding box,
  // which collapses volumes when no agents exist.
  auto* env = dynamic_cast<UniformGridEnvironment*>(simulation.GetEnvironment());
  auto* sp_early = simulation.GetParam()->Get<SimParam>();
  // Box must exceed max possible neighbor search radius (largest diameter +
  // growth headroom). Factor of 2 keeps it safe against any transient growth.
  env->SetBoxLength(static_cast<int32_t>(sp_early->division_diameter * 2));
  env->SetDetermineSimSize(false);

  auto* sp = simulation.GetParam()->Get<SimParam>();

  // --- Auto-register conditional visualization ---
  {
    auto* param = const_cast<Param*>(simulation.GetParam());
    // Keratinocyte exported as uniform spheres (no per-stratum coloring)
    param->visualize_agents["Keratinocyte"] = {};
    // ImmuneCell only exists in wound scenarios
    if (sp->wound_enabled) {
      param->visualize_agents["ImmuneCell"] = {"immune_type_"};
    }
    // TumorCell + Tumor field only when tumor module is enabled
    if (sp->tumor_enabled) {
      param->visualize_agents["TumorCell"] = {"tumor_type_"};
      Param::VisualizeDiffusion vd;
      vd.name = "Tumor";
      param->visualize_diffusion.push_back(vd);
    }
  }

  // --- Volume manager ---
  auto* vm = VolumeManager::Get();
  vm->Init(sp->volume_z_spinous, sp->volume_z_granular, sp->volume_z_cornified,
           sp->dermal_z_papillary, sp->dermal_z_reticular);
  // Enable agents in all epidermal layers for wound healing visualization.
  // Init() defaults suprabasal to continuum-only (strict UWYN); the
  // simulation overrides this so cells traverse all layers and dissolve
  // individually via per-cell UWYN handoff in Differentiation.
  vm->SetAgentsEnabled(kSpinous, true);
  vm->SetAgentsEnabled(kGranular, true);
  vm->SetAgentsEnabled(kCornified, true);

  // --- Register PDE fields and configure sub-cycling ---
  CompositeField fields;
  RegisterFields(&simulation, sp, fields);

  // --- Derived composite fields (lightweight arrays, not DiffusionGrids) ---
  auto* ref_grid = simulation.GetResourceManager()
      ->GetDiffusionGrid(fields::kOxygen);
  size_t ref_res = ref_grid->GetResolution();
  auto ref_dims = ref_grid->GetDimensions();
  real_t ref_lo = static_cast<real_t>(ref_dims[0]);
  real_t ref_box_len = ref_grid->GetBoxLength();

  auto ecm_quality = std::make_unique<DerivedField>(
      fields::kECMQuality, ref_res, ref_lo, ref_box_len);
  auto tissue_viability = std::make_unique<DerivedField>(
      fields::kTissueViability, ref_res, ref_lo, ref_box_len);
  auto wound_microenv = std::make_unique<DerivedField>(
      fields::kWoundMicroenv, ref_res, ref_lo, ref_box_len);

  // Set global accessors for agent reads (util.h)
  g_ecm_quality = ecm_quality.get();
  g_tissue_viability = tissue_viability.get();
  g_wound_microenv = wound_microenv.get();

  // --- Custom cell-cell force with adhesion ---
  auto* scheduler = simulation.GetScheduler();
  auto* custom_force = new CellCellForce();
  auto* mech_op = scheduler->GetOps("mechanical forces")[0];
  auto* force_impl = mech_op->GetImplementation<MechanicalForcesOp>();
  force_impl->SetInteractionForce(custom_force);

  // --- Register all scheduled operations ---
  RegisterOperations(&simulation, sp, &fields,
                     ecm_quality.get(), tissue_viability.get(),
                     wound_microenv.get());

  simulation.GetScheduler()->Simulate(sp->num_steps);

  // --- Final tumor cleanup: convert remaining agents to field ---
  if (sp->tumor_enabled) {
    auto* rm_cleanup = simulation.GetResourceManager();
    auto* param = simulation.GetParam();
    real_t margin = 1e-3;
    auto* tumor_grid = rm_cleanup->GetDiffusionGrid(fields::kTumor);
    auto* stratum_grid = rm_cleanup->GetDiffusionGrid(fields::kStratum);
    int converted = 0;

    rm_cleanup->ForEachAgent([&](Agent* agent) {
      auto* tumor = dynamic_cast<TumorCell*>(agent);
      if (!tumor) return;

      Real3 pos = tumor->GetPosition();
      for (int i = 0; i < 3; i++) {
        pos[i] = std::max(param->min_bound + margin,
                          std::min(pos[i], param->max_bound - margin));
      }

      size_t tidx = tumor_grid->GetBoxIndex(pos);
      real_t t_cur = tumor_grid->GetConcentration(tidx);
      if (t_cur < 0.5) {
        tumor_grid->ChangeConcentrationBy(tidx, 1.0 - t_cur);
      }

      if (sp->tumor_stratum_value > 0) {
        size_t sidx = stratum_grid->GetBoxIndex(pos);
        real_t s_cur = stratum_grid->GetConcentration(sidx);
        if (sp->tumor_stratum_value > s_cur) {
          stratum_grid->ChangeConcentrationBy(
              sidx, sp->tumor_stratum_value - s_cur);
        }
      }

      tumor->RemoveFromSimulation();
      converted++;
    });

    simulation.GetScheduler()->Simulate(1);
    std::cout << "Final tumor cleanup: " << converted
              << " agents converted to field" << std::endl;
  }

  // --- Summary ---
  auto* rm = simulation.GetResourceManager();
  int n_agents = 0;
  rm->ForEachAgent([&](Agent*) { n_agents++; });
  std::cout << "\n=== Simulation Summary (UWYN: continuum + event-driven agents) ==="
            << std::endl;
  std::cout << "Active agents: " << n_agents << std::endl;
  std::cout << "Continuum fields: " << fields.Size() << " PDE channels"
            << std::endl;
  std::cout << "Skin simulation completed successfully!" << std::endl;

  // Flush and close metrics CSV before opening
  auto* metrics_op = scheduler->GetOps("MetricsExporter")[0];
  auto* metrics_impl = metrics_op->GetImplementation<MetricsExporter>();
  if (metrics_impl) metrics_impl->Close();

  // Auto-open metrics CSV in default app (configurable in bdm.toml)
  if (sp->metrics_autoopen && sp->metrics_interval > 0) {
    std::string csv = simulation.GetOutputDir() + "/metrics.csv";
    std::string cmd = "xdg-open \"" + csv + "\" &";
    std::system(cmd.c_str());
  }

  return 0;
}

}  // namespace skibidy
}  // namespace bdm

#endif  // SKIBIDY_H_
