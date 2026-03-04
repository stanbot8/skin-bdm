#include "test_helpers.h"
#include "wound/wound_event.h"
#include "core/fused_source.h"
#include "core/metrics.h"

namespace bdm {
namespace skibidy {

TEST(WoundTest, MarginCellCount) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  sp->wound.trigger_step = 0;  // fire immediately

  // Set up fields (WoundEvent needs them for ApplyWoundAll + InitFromFields)
  auto* vm = VolumeManager::Get();
  vm->Init(sp->volume_z_spinous, sp->volume_z_granular,
           sp->volume_z_cornified);

  CompositeField fields;
  fields.Add(std::make_unique<CalciumPDE>());
  fields.Add(std::make_unique<KgfPDE>());
  fields.Add(std::make_unique<OxygenPDE>());
  fields.Add(std::make_unique<StratumPDE>());
  fields.Add(std::make_unique<ScarPDE>());
  fields.Add(std::make_unique<WaterPDE>());
  if (sp->inflammation.split_inflammation_enabled) {
    fields.Add(std::make_unique<ProInflammatoryPDE>(sp));
    fields.Add(std::make_unique<AntiInflammatoryPDE>(sp));
  } else {
    fields.Add(std::make_unique<InflammationPDE>(sp));
  }
  if (sp->fibronectin.enabled) fields.Add(std::make_unique<FibronectinPDE>(sp));
  if (sp->mmp.enabled) fields.Add(std::make_unique<MMPPDE>(sp));
  if (sp->fibroblast.enabled) {
    fields.Add(std::make_unique<TGFBetaPDE>(sp));
    fields.Add(std::make_unique<CollagenPDE>(sp));
  }
  if (sp->ros.enabled) fields.Add(std::make_unique<ROSPDE>(sp));
  if (sp->mechanotransduction.enabled) fields.Add(std::make_unique<StiffnessPDE>(sp));
  if (sp->lymphatic.enabled) {
    fields.Add(std::make_unique<LymphaticPDE>(sp));
    fields.Add(std::make_unique<EdemaPDE>(sp));
  }
  if (sp->bioelectric.enabled) fields.Add(std::make_unique<VoltagePDE>(sp));
  fields.InitAll(sim);

  // Register and schedule WoundEvent
  OperationRegistry::GetInstance()->AddOperationImpl(
      "WoundEvent_test", OpComputeTarget::kCpu, new WoundEvent(&fields));
  auto* wound_op = NewOperation("WoundEvent_test");
  sim->GetScheduler()->ScheduleOp(wound_op, OpType::kPreSchedule);

  sim->GetScheduler()->Simulate(1);

  // Expected: round(2*pi*r / d), min 6
  // r=5, d=5 -> round(2*pi*5/5) = round(6.28) = 6
  int expected = static_cast<int>(
      std::round(2.0 * M_PI * sp->wound.radius / sp->division_diameter));
  if (expected < 6) expected = 6;

  EXPECT_EQ(sim->GetResourceManager()->GetNumAgents(),
            static_cast<size_t>(expected));
  delete sim;
}

TEST(WoundResolutionTest, SafetyTimeout) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  sp->wound.trigger_step = 0;
  sp->num_steps = 25;  // safety net fires at num_steps - 1 = step 24

  auto* vm = VolumeManager::Get();
  vm->Init(sp->volume_z_spinous, sp->volume_z_granular,
           sp->volume_z_cornified);

  CompositeField fields;
  fields.Add(std::make_unique<CalciumPDE>());
  fields.Add(std::make_unique<KgfPDE>());
  fields.Add(std::make_unique<OxygenPDE>());
  fields.Add(std::make_unique<StratumPDE>());
  fields.Add(std::make_unique<ScarPDE>());
  fields.Add(std::make_unique<WaterPDE>());
  if (sp->inflammation.split_inflammation_enabled) {
    fields.Add(std::make_unique<ProInflammatoryPDE>(sp));
    fields.Add(std::make_unique<AntiInflammatoryPDE>(sp));
  } else {
    fields.Add(std::make_unique<InflammationPDE>(sp));
  }
  if (sp->fibronectin.enabled) fields.Add(std::make_unique<FibronectinPDE>(sp));
  if (sp->mmp.enabled) fields.Add(std::make_unique<MMPPDE>(sp));
  if (sp->fibroblast.enabled) {
    fields.Add(std::make_unique<TGFBetaPDE>(sp));
    fields.Add(std::make_unique<CollagenPDE>(sp));
  }
  if (sp->ros.enabled) fields.Add(std::make_unique<ROSPDE>(sp));
  if (sp->mechanotransduction.enabled) fields.Add(std::make_unique<StiffnessPDE>(sp));
  if (sp->lymphatic.enabled) {
    fields.Add(std::make_unique<LymphaticPDE>(sp));
    fields.Add(std::make_unique<EdemaPDE>(sp));
  }
  if (sp->bioelectric.enabled) fields.Add(std::make_unique<VoltagePDE>(sp));
  fields.InitAll(sim);

  // Schedule wound at step 0 to spawn cells
  OperationRegistry::GetInstance()->AddOperationImpl(
      "WoundEvent_res", OpComputeTarget::kCpu, new WoundEvent(&fields));
  sim->GetScheduler()->ScheduleOp(
      NewOperation("WoundEvent_res"), OpType::kPreSchedule);

  // Schedule resolution
  OperationRegistry::GetInstance()->AddOperationImpl(
      "WoundRes_test", OpComputeTarget::kCpu, new WoundResolution());
  sim->GetScheduler()->ScheduleOp(
      NewOperation("WoundRes_test"), OpType::kPreSchedule);

  // Step 1: wound fires, cells spawned
  sim->GetScheduler()->Simulate(1);
  size_t before = sim->GetResourceManager()->GetNumAgents();
  EXPECT_GT(before, 0u);

  // Run past safety timeout (step >= num_steps - 1 = 24).
  // The safety net dissolves the vast majority of agents. A few
  // daughters from same-step division may survive since WoundResolution
  // marks resolved_ and stops re-entering.
  sim->GetScheduler()->Simulate(24);

  size_t after = sim->GetResourceManager()->GetNumAgents();
  EXPECT_LT(after, before) << "Safety net should dissolve most agents";

  delete sim;
}

TEST(WoundInflammationSourceTest, OpenWoundProducesInflammation) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  SanitizeForUnitTest(sp);
  sp->wound.inflammation_source_rate = 0.01;  // strong rate for test visibility
  sp->wound.trigger_step = 0;
  sp->wound.center_x = 25;
  sp->wound.center_y = 25;
  sp->wound.radius = 10;
  sp->dermis.enabled = false;

  auto* rm = sim->GetResourceManager();

  // Register required grids (FusedWoundSourceOp needs all 5 core + inflammation)
  int res = sp->grid_resolution;
  ModelInitializer::DefineSubstance(fields::kVascularId, fields::kVascular, 0, 0, res);
  ModelInitializer::DefineSubstance(fields::kOxygenId, fields::kOxygen, 0, 0, res);
  ModelInitializer::DefineSubstance(fields::kStratumId, fields::kStratum, 0, 0, res);
  ModelInitializer::DefineSubstance(fields::kCalciumId, fields::kCalcium, 0, 0, res);
  ModelInitializer::DefineSubstance(fields::kWaterId, fields::kWater, 0, 0, res);
  ModelInitializer::DefineSubstance(fields::kInflammationId, fields::kInflammation, 0, 0, res);

  sim->GetEnvironment()->Update();

  auto* vasc_grid = rm->GetDiffusionGrid(fields::kVascularId);
  auto* stratum_grid = rm->GetDiffusionGrid(fields::kStratumId);
  auto* infl_grid = rm->GetDiffusionGrid(fields::kInflammationId);
  vasc_grid->Initialize();
  rm->GetDiffusionGrid(fields::kOxygenId)->Initialize();
  stratum_grid->Initialize();
  rm->GetDiffusionGrid(fields::kCalciumId)->Initialize();
  rm->GetDiffusionGrid(fields::kWaterId)->Initialize();
  infl_grid->Initialize();

  // Skip diffusion solver
  vasc_grid->SetTimeStep(1e30);
  rm->GetDiffusionGrid(fields::kOxygenId)->SetTimeStep(1e30);
  stratum_grid->SetTimeStep(1e30);
  rm->GetDiffusionGrid(fields::kCalciumId)->SetTimeStep(1e30);
  rm->GetDiffusionGrid(fields::kWaterId)->SetTimeStep(1e30);
  infl_grid->SetTimeStep(1e30);

  size_t n = vasc_grid->GetNumBoxes();
  GridContext ctx(vasc_grid, sp);

  // Set vascular=1 everywhere and stratum=0 in wound (open wound)
  for (size_t i = 0; i < n; i++) {
    vasc_grid->ChangeConcentrationBy(i, 1.0);
  }

  // Advance past wound_trigger_step=0
  sim->GetScheduler()->Simulate(1);

  // Run fused source op
  FusedWoundSourceOp op;
  op();

  // Verify inflammation increased in epidermal wound voxels
  real_t total_infl = 0;
  int wound_voxels = 0;
  for (size_t i = 0; i < n; i++) {
    real_t c = infl_grid->GetConcentration(i);
    if (c > 1e-10) {
      total_infl += c;
      wound_voxels++;
    }
  }

  EXPECT_GT(wound_voxels, 0) << "No wound voxels received inflammation";
  EXPECT_GT(total_infl, 0) << "No inflammation produced by wound source";

  delete sim;
}

TEST(WoundInflammationSourceTest, ClosedWoundNoSource) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  SanitizeForUnitTest(sp);
  sp->wound.inflammation_source_rate = 0.01;
  sp->wound.trigger_step = 0;
  sp->wound.center_x = 25;
  sp->wound.center_y = 25;
  sp->wound.radius = 10;
  sp->dermis.enabled = false;

  auto* rm = sim->GetResourceManager();

  int res = sp->grid_resolution;
  ModelInitializer::DefineSubstance(fields::kVascularId, fields::kVascular, 0, 0, res);
  ModelInitializer::DefineSubstance(fields::kOxygenId, fields::kOxygen, 0, 0, res);
  ModelInitializer::DefineSubstance(fields::kStratumId, fields::kStratum, 0, 0, res);
  ModelInitializer::DefineSubstance(fields::kCalciumId, fields::kCalcium, 0, 0, res);
  ModelInitializer::DefineSubstance(fields::kWaterId, fields::kWater, 0, 0, res);
  ModelInitializer::DefineSubstance(fields::kInflammationId, fields::kInflammation, 0, 0, res);

  sim->GetEnvironment()->Update();

  auto* vasc_grid = rm->GetDiffusionGrid(fields::kVascularId);
  auto* stratum_grid = rm->GetDiffusionGrid(fields::kStratumId);
  auto* infl_grid = rm->GetDiffusionGrid(fields::kInflammationId);
  vasc_grid->Initialize();
  rm->GetDiffusionGrid(fields::kOxygenId)->Initialize();
  stratum_grid->Initialize();
  rm->GetDiffusionGrid(fields::kCalciumId)->Initialize();
  rm->GetDiffusionGrid(fields::kWaterId)->Initialize();
  infl_grid->Initialize();

  vasc_grid->SetTimeStep(1e30);
  rm->GetDiffusionGrid(fields::kOxygenId)->SetTimeStep(1e30);
  stratum_grid->SetTimeStep(1e30);
  rm->GetDiffusionGrid(fields::kCalciumId)->SetTimeStep(1e30);
  rm->GetDiffusionGrid(fields::kWaterId)->SetTimeStep(1e30);
  infl_grid->SetTimeStep(1e30);

  size_t n = vasc_grid->GetNumBoxes();

  // Set vascular=1, stratum=1 everywhere (fully healed)
  for (size_t i = 0; i < n; i++) {
    vasc_grid->ChangeConcentrationBy(i, 1.0);
    stratum_grid->ChangeConcentrationBy(i, 1.0);
  }

  sim->GetScheduler()->Simulate(1);

  FusedWoundSourceOp op;
  op();

  // With stratum=1 everywhere, wound source should produce nothing
  real_t total_infl = 0;
  for (size_t i = 0; i < n; i++) {
    total_infl += infl_grid->GetConcentration(i);
  }

  EXPECT_NEAR(total_infl, 0, 1e-8) << "Closed wound should not produce inflammation";

  delete sim;
}

}  // namespace skibidy
}  // namespace bdm
