#include "test_helpers.h"

#include "core/voxel_env.h"
#include "immune/neutrophil_behavior.h"
#include "immune/macrophage_behavior.h"
#include "immune/immune_response.h"
#include "diabetic/baseline_inflammation.h"
#include "fibroblast/fibroblast.h"
#include "fibroblast/fibroblast_behavior.h"
#include "fibroblast/fibroblast_recruitment.h"

namespace bdm {
namespace skibidy {

TEST(SplitInflammationTest, ProAntiFieldsRegistered) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  sp->inflammation.split_inflammation_enabled = true;

  ProInflammatoryPDE pro(sp);
  pro.Init(sim);
  AntiInflammatoryPDE anti(sp);
  anti.Init(sim);

  auto* rm = sim->GetResourceManager();
  EXPECT_NE(rm->GetDiffusionGrid(fields::kProInflammatoryId), nullptr);
  EXPECT_NE(rm->GetDiffusionGrid(fields::kAntiInflammatoryId), nullptr);

  delete sim;
}

TEST(SplitInflammationTest, NeutrophilWritesToProField) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  SanitizeForUnitTest(sp);
  sp->inflammation.split_inflammation_enabled = true;

  ProInflammatoryPDE pro(sp);
  pro.Init(sim);
  AntiInflammatoryPDE anti(sp);
  anti.Init(sim);
  sim->GetScheduler()->Simulate(1);

  auto* rm = sim->GetResourceManager();
  auto* pro_grid = rm->GetDiffusionGrid(fields::kProInflammatoryId);
  real_t initial = pro_grid->GetValue({15, 15, 1.5});

  auto* cell = new ImmuneCell({15, 15, 1.5});
  cell->SetDiameter(3);
  cell->SetImmuneCellType(kNeutrophil);
  cell->AddBehavior(new NeutrophilBehavior());
  rm->AddAgent(cell);

  sim->GetScheduler()->Simulate(3);

  real_t after = pro_grid->GetValue({15, 15, 1.5});
  EXPECT_GT(after, initial);

  delete sim;
}

TEST(SplitInflammationTest, M2WritesToAntiField) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  SanitizeForUnitTest(sp);
  sp->inflammation.split_inflammation_enabled = true;
  sp->immune.macrophage_lifespan = 100;

  ProInflammatoryPDE pro(sp);
  pro.Init(sim);
  AntiInflammatoryPDE anti(sp);
  anti.Init(sim);
  sim->GetScheduler()->Simulate(1);

  auto* rm = sim->GetResourceManager();
  auto* anti_grid = rm->GetDiffusionGrid(fields::kAntiInflammatoryId);
  real_t initial = anti_grid->GetValue({15, 15, 1.5});

  auto* cell = new ImmuneCell({15, 15, 1.5});
  cell->SetDiameter(3);
  cell->SetImmuneCellType(kMacrophage);
  cell->SetState(kM2Resolving);
  cell->AddBehavior(new MacrophageBehavior());
  rm->AddAgent(cell);

  sim->GetScheduler()->Simulate(3);

  real_t after = anti_grid->GetValue({15, 15, 1.5});
  EXPECT_GT(after, initial);

  delete sim;
}

// ===========================================================================
// Extension 5: Diabetic mode
// ===========================================================================
TEST(DiabeticTest, ExtendedM1Duration) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  SanitizeForDiabeticTest(sp);
  sp->diabetic.mode = true;
  sp->immune.macrophage_m1_duration = 3;
  sp->diabetic.m1_duration_factor = 3.0;
  sp->immune.macrophage_lifespan = 100;
  sp->m1_transition_threshold = 0;  // disable cytokine path, test timer only

  SetupInflammation(sim);
  sim->GetScheduler()->Simulate(1);

  auto* cell = new ImmuneCell({15, 15, 1.5});
  cell->SetDiameter(3);
  cell->SetImmuneCellType(kMacrophage);
  cell->SetState(kM1Active);
  cell->AddBehavior(new MacrophageBehavior());
  sim->GetResourceManager()->AddAgent(cell);

  // Normal M1 duration=3, diabetic factor=3 -> effective=9
  // After 5 steps: age=5, effective M1 duration=9, still M1
  sim->GetScheduler()->Simulate(5);
  ImmuneState state_at_5 = kM1Active;
  sim->GetResourceManager()->ForEachAgent([&](Agent* a) {
    auto* ic = dynamic_cast<ImmuneCell*>(a);
    if (ic) state_at_5 = ic->GetState();
  });
  EXPECT_EQ(state_at_5, kM1Active);

  // After 5 more (total 10): age=10 > 9, transitions to M2
  sim->GetScheduler()->Simulate(5);
  ImmuneState state_at_10 = kM1Active;
  sim->GetResourceManager()->ForEachAgent([&](Agent* a) {
    auto* ic = dynamic_cast<ImmuneCell*>(a);
    if (ic) state_at_10 = ic->GetState();
  });
  EXPECT_EQ(state_at_10, kM2Resolving);

  delete sim;
}

TEST(DiabeticTest, ReducedResolutionRate) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  SanitizeForDiabeticTest(sp);
  sp->diabetic.mode = true;
  sp->diabetic.resolution_factor = 0.3;
  sp->immune.macrophage_lifespan = 100;

  SetupInflammation(sim);
  sim->GetScheduler()->Simulate(1);

  // Seed inflammation at a voxel
  auto* rm = sim->GetResourceManager();
  auto* infl_grid = GetInflammationGrid(sim);
  Real3 pos = {15, 15, 1.5};
  Real3 qpos = pos;
  real_t margin = 1e-3;
  auto* param = sim->GetParam();
  for (int i = 0; i < 3; i++) {
    qpos[i] = std::max(param->min_bound + margin,
                       std::min(qpos[i], param->max_bound - margin));
  }
  size_t idx = infl_grid->GetBoxIndex(qpos);
  infl_grid->ChangeConcentrationBy(idx, 1.0);

  real_t before = infl_grid->GetConcentration(idx);

  // M2 macrophage with diabetic mode (not added to RM -- direct behavior call)
  auto* cell = new ImmuneCell({15, 15, 1.5});
  cell->SetDiameter(3);
  cell->SetImmuneCellType(kMacrophage);
  cell->SetState(kM2Resolving);

  // Run behavior directly to test exact consumption amount
  MacrophageBehavior beh;
  beh.Run(cell);

  real_t after = infl_grid->GetConcentration(idx);
  real_t consumed = before - after;

  // Effective rate = resolution_rate * diabetic_factor = 0.005 * 0.3 = 0.0015
  real_t expected = sp->immune.resolution_rate * sp->diabetic.resolution_factor;
  EXPECT_NEAR(consumed, expected, 1e-6);

  delete cell;
  delete sim;
}

TEST(DiabeticTest, KeratinocyteProliferation) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  SanitizeForDiabeticTest(sp);
  sp->diabetic.mode = true;
  sp->diabetic.prolif_factor = 0.0;  // completely suppress proliferation

  SetupAllFields(sim);
  sim->GetScheduler()->Simulate(1);
  g_voxel_env.Fill(sim);

  auto* cell = new Keratinocyte({15, 15, 2});
  cell->SetDiameter(4);
  cell->SetStratum(kBasal);
  cell->SetDivisionsLeft(5);
  cell->SetCyclePhase(kG1);
  cell->SetPhaseElapsed(100);  // very high -> p would be >1 normally

  // Direct behavior call: factor=0 means p*0=0, no G1->S transition
  BasalDivision div;
  for (int i = 0; i < 50; i++) {
    div.Run(cell);
    EXPECT_EQ(cell->GetCyclePhase(), kG1)
        << "Cell should never leave G1 with prolif_factor=0";
  }

  delete cell;
  delete sim;
}

TEST(DiabeticTest, KeratinocyteMigration) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  SanitizeForDiabeticTest(sp);
  sp->wound.enabled = true;
  sp->migration_enabled = true;
  sp->wound.center_x = 15;
  sp->wound.center_y = 15;
  sp->wound.radius = 10;
  SetupAllFields(sim);
  sim->GetScheduler()->Simulate(1);
  g_voxel_env.Fill(sim);

  // Normal cell (diabetic off): measure tractor force
  sp->diabetic.mode = false;
  auto* cell1 = new Keratinocyte({20, 15, 2});
  cell1->SetDiameter(4);
  cell1->SetStratum(kBasal);

  Migration mig;
  mig.Run(cell1);
  Real3 force_normal = cell1->GetTractorForce();
  real_t mag_normal = std::sqrt(
      force_normal[0] * force_normal[0] + force_normal[1] * force_normal[1]);

  // Diabetic cell at same position
  sp->diabetic.mode = true;
  sp->diabetic.migration_factor = 0.5;

  auto* cell2 = new Keratinocyte({20, 15, 2});
  cell2->SetDiameter(4);
  cell2->SetStratum(kBasal);

  mig.Run(cell2);
  Real3 force_diabetic = cell2->GetTractorForce();
  real_t mag_diabetic = std::sqrt(
      force_diabetic[0] * force_diabetic[0] +
      force_diabetic[1] * force_diabetic[1]);

  EXPECT_GT(mag_normal, 0) << "Normal cell should have nonzero migration";
  EXPECT_NEAR(mag_diabetic, mag_normal * 0.5, mag_normal * 0.01);

  delete cell1;
  delete cell2;
  delete sim;
}

TEST(DiabeticTest, FibroblastDysfunction) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  SanitizeForDiabeticTest(sp);
  sp->diabetic.mode = true;
  sp->fibroblast.enabled = true;
  sp->fibroblast.activation_delay = 10;
  sp->diabetic.fibroblast_activation_factor = 3.0;
  sp->fibroblast.activation_threshold = 999;  // prevent TGF-beta trigger
  sp->fibroblast.lifespan = 10000;

  SetupAllFields(sim);
  sim->GetScheduler()->Simulate(1);

  auto* cell = new Fibroblast({15, 15, -3});
  cell->SetDiameter(4);
  cell->SetFibroblastState(kFibroQuiescent);

  FibroblastBehavior beh;
  // Run 15 steps: base delay=10, effective=30 with factor=3
  // At step 15 the cell should still be quiescent
  for (int i = 0; i < 15; i++) {
    beh.Run(cell);
  }
  EXPECT_EQ(cell->GetFibroblastState(), kFibroQuiescent)
      << "Fibroblast should stay quiescent with effective delay=30";

  delete cell;
  delete sim;
}

TEST(DiabeticTest, ExcessiveNeutrophils) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  SanitizeForDiabeticTest(sp);
  sp->wound.enabled = true;
  sp->wound.trigger_step = 0;
  sp->immune.neutrophil_spawn_delay = 1;
  sp->immune.neutrophil_spawn_waves = 1;       // single wave for count test
  sp->immune.neutrophil_apoptosis_rate = 0;    // disable stochastic death
  sp->diabetic.mode = true;
  sp->diabetic.neutrophil_factor = 2.0;
  sp->diabetic.neutrophil_waves_factor = 1.0;   // isolate count test from wave scaling
  sp->diabetic.neutrophil_window_factor = 1.0;  // isolate count test from window scaling

  SetupAllFields(sim);

  // Register ImmuneResponse as a scheduled op so agents get committed
  auto* resp = new ImmuneResponse();
  OperationRegistry::GetInstance()->AddOperationImpl(
      "ImmuneResponse_diabetic", OpComputeTarget::kCpu, resp);
  sim->GetScheduler()->ScheduleOp(
      NewOperation("ImmuneResponse_diabetic"), OpType::kPreSchedule);

  sim->GetScheduler()->Simulate(3);

  int diabetic_count = 0;
  sim->GetResourceManager()->ForEachAgent([&](Agent* a) {
    auto* ic = dynamic_cast<ImmuneCell*>(a);
    if (ic && ic->GetImmuneCellType() == kNeutrophil) diabetic_count++;
  });

  // Base count: round(2*PI*5/3) = 10, diabetic factor 2.0 -> ceil(10*2) = 20
  int base_count = static_cast<int>(std::round(
      2.0 * M_PI * sp->wound.radius / 3.0));
  if (base_count < 6) base_count = 6;
  int expected = static_cast<int>(std::ceil(base_count * 2.0));

  EXPECT_EQ(diabetic_count, expected);

  delete sim;
}

TEST(DiabeticTest, ExtendedNeutrophilLifespan) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  SanitizeForDiabeticTest(sp);
  sp->diabetic.mode = true;
  sp->immune.neutrophil_lifespan = 10;
  sp->diabetic.neutrophil_lifespan_factor = 2.0;
  sp->immune.neutrophil_apoptosis_rate = 0;    // disable stochastic death
  sp->immune.neutrophil_min_survival = 9999;   // disable stochastic death

  SetupInflammation(sim);
  sim->GetScheduler()->Simulate(1);

  // Neutrophil with age=15: past base lifespan (10) but within diabetic (20)
  auto* cell = new ImmuneCell({15, 15, 1.5});
  cell->SetDiameter(3);
  cell->SetImmuneCellType(kNeutrophil);
  // Manually set age past base lifespan
  for (int i = 0; i < 15; i++) cell->IncrementAge();

  NeutrophilBehavior beh;
  beh.Run(cell);

  // Cell should still be alive (age=16 < effective lifespan=20)
  // RemoveFromSimulation() sets removal flag; check it's not flagged.
  // We verify by checking the cell is still accessible (not removed).
  // If removed, IncrementAge would have been skipped (return early).
  EXPECT_EQ(cell->GetAge(), 16);

  delete cell;
  delete sim;
}

}  // namespace skibidy
}  // namespace bdm
