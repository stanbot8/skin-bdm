#include "test_helpers.h"
#include "immune/neutrophil_behavior.h"
#include "immune/macrophage_behavior.h"
#include "immune/immune_response.h"
#include "fibroblast/fibroblast.h"
#include "fibroblast/fibroblast_behavior.h"
#include "diabetic/baseline_inflammation.h"

namespace bdm {
namespace skibidy {

TEST(ImmuneCellTest, TypeAndStateDefaults) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* cell = new ImmuneCell({15, 15, 1.5});
  cell->SetDiameter(3);

  EXPECT_EQ(cell->GetImmuneCellType(), kNeutrophil);
  EXPECT_EQ(cell->GetState(), kM1Active);
  EXPECT_EQ(cell->GetAge(), 0);

  cell->SetImmuneCellType(kMacrophage);
  EXPECT_EQ(cell->GetImmuneCellType(), kMacrophage);

  cell->SetState(kM2Resolving);
  EXPECT_EQ(cell->GetState(), kM2Resolving);

  cell->IncrementAge();
  cell->IncrementAge();
  EXPECT_EQ(cell->GetAge(), 2);

  delete cell;
  delete sim;
}

TEST(NeutrophilBehaviorTest, HardCeilingLifespan) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  sp->immune.neutrophil_lifespan = 5;
  sp->immune.neutrophil_apoptosis_rate = 0;  // disable stochastic death, test ceiling
  SanitizeForUnitTest(sp);

  SetupInflammation(sim);
  sim->GetScheduler()->Simulate(1);  // populate grid

  auto* cell = new ImmuneCell({15, 15, 1.5});
  cell->SetDiameter(3);
  cell->SetImmuneCellType(kNeutrophil);
  cell->AddBehavior(new NeutrophilBehavior());
  sim->GetResourceManager()->AddAgent(cell);

  // Run 5 steps: age reaches 5, lifespan check is > 5, still alive
  sim->GetScheduler()->Simulate(5);
  EXPECT_EQ(sim->GetResourceManager()->GetNumAgents(), 1u);

  // Step 6: age = 6 > 5, removed
  sim->GetScheduler()->Simulate(1);
  EXPECT_EQ(sim->GetResourceManager()->GetNumAgents(), 0u);

  delete sim;
}

TEST(MacrophageBehaviorTest, M1ToM2TimerFallback) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  sp->immune.macrophage_m1_duration = 3;
  sp->immune.macrophage_lifespan = 100;
  sp->m1_transition_threshold = 0;  // disable cytokine path, test timer only
  SanitizeForUnitTest(sp);

  SetupInflammation(sim);
  sim->GetScheduler()->Simulate(1);

  auto* cell = new ImmuneCell({15, 15, 1.5});
  cell->SetDiameter(3);
  cell->SetImmuneCellType(kMacrophage);
  cell->SetState(kM1Active);
  cell->AddBehavior(new MacrophageBehavior());
  sim->GetResourceManager()->AddAgent(cell);

  // After 3 steps: age = 3, check is > 3, still M1
  sim->GetScheduler()->Simulate(3);
  // Must iterate to read state
  ImmuneState state_at_3 = kM1Active;
  sim->GetResourceManager()->ForEachAgent([&](Agent* a) {
    auto* ic = dynamic_cast<ImmuneCell*>(a);
    if (ic) state_at_3 = ic->GetState();
  });
  EXPECT_EQ(state_at_3, kM1Active);

  // Step 4: age = 4 > 3, transitions to M2
  sim->GetScheduler()->Simulate(1);
  ImmuneState state_at_4 = kM1Active;
  sim->GetResourceManager()->ForEachAgent([&](Agent* a) {
    auto* ic = dynamic_cast<ImmuneCell*>(a);
    if (ic) state_at_4 = ic->GetState();
  });
  EXPECT_EQ(state_at_4, kM2Resolving);

  delete sim;
}

TEST(MacrophageBehaviorTest, CytokineTriggeredTransition) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  sp->immune.macrophage_m1_duration = 9999;         // disable timer
  sp->m1_transition_threshold = 0.5;
  sp->m1_transition_min_age = 3;             // allow cytokine trigger after 3 steps
  sp->immune.macrophage_lifespan = 100;
  sp->inflammation.decay = 0;                // no decay, we control concentration
  SanitizeForUnitTest(sp);

  SetupInflammation(sim);
  sim->GetScheduler()->Simulate(1);

  // Seed inflammation above threshold at cell position
  auto* infl_grid = GetInflammationGrid(sim);
  Real3 pos = {15, 15, 1.5};
  size_t idx = infl_grid->GetBoxIndex(pos);
  infl_grid->ChangeConcentrationBy(idx, 2.0);  // well above 0.5

  auto* cell = new ImmuneCell(pos);
  cell->SetDiameter(3);
  cell->SetImmuneCellType(kMacrophage);
  cell->SetState(kM1Active);
  cell->AddBehavior(new MacrophageBehavior());
  sim->GetResourceManager()->AddAgent(cell);

  // Run 5 steps with high inflammation -- stays M1
  sim->GetScheduler()->Simulate(5);
  ImmuneState state_high = kM1Active;
  sim->GetResourceManager()->ForEachAgent([&](Agent* a) {
    auto* ic = dynamic_cast<ImmuneCell*>(a);
    if (ic && ic->GetImmuneCellType() == kMacrophage)
      state_high = ic->GetState();
  });
  EXPECT_EQ(state_high, kM1Active);

  // Clear inflammation below threshold
  real_t current = infl_grid->GetConcentration(idx);
  infl_grid->ChangeConcentrationBy(idx, -current);  // set to 0

  // Run 1 more step -- now below threshold and past min_age -> M2
  sim->GetScheduler()->Simulate(1);
  ImmuneState state_low = kM1Active;
  sim->GetResourceManager()->ForEachAgent([&](Agent* a) {
    auto* ic = dynamic_cast<ImmuneCell*>(a);
    if (ic && ic->GetImmuneCellType() == kMacrophage)
      state_low = ic->GetState();
  });
  EXPECT_EQ(state_low, kM2Resolving);

  delete sim;
}

TEST(MacrophageBehaviorTest, MinAgeGuardPreventsEarlyTransition) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  sp->immune.macrophage_m1_duration = 9999;         // disable timer
  sp->m1_transition_threshold = 999.0;       // always above local inflammation
  sp->m1_transition_min_age = 10;            // guard: must be in M1 for 10 steps
  sp->immune.macrophage_lifespan = 100;
  SanitizeForUnitTest(sp);

  SetupInflammation(sim);
  sim->GetScheduler()->Simulate(1);

  auto* cell = new ImmuneCell({15, 15, 1.5});
  cell->SetDiameter(3);
  cell->SetImmuneCellType(kMacrophage);
  cell->SetState(kM1Active);
  cell->AddBehavior(new MacrophageBehavior());
  sim->GetResourceManager()->AddAgent(cell);

  // After 5 steps: state_age=5 < min_age=10, stays M1 despite threshold
  sim->GetScheduler()->Simulate(5);
  ImmuneState state_early = kM1Active;
  sim->GetResourceManager()->ForEachAgent([&](Agent* a) {
    auto* ic = dynamic_cast<ImmuneCell*>(a);
    if (ic) state_early = ic->GetState();
  });
  EXPECT_EQ(state_early, kM1Active);

  // After 6 more (total 11): state_age=11 > 10, now transitions
  sim->GetScheduler()->Simulate(6);
  ImmuneState state_after = kM1Active;
  sim->GetResourceManager()->ForEachAgent([&](Agent* a) {
    auto* ic = dynamic_cast<ImmuneCell*>(a);
    if (ic) state_after = ic->GetState();
  });
  EXPECT_EQ(state_after, kM2Resolving);

  delete sim;
}

TEST(NeutrophilBehaviorTest, CytokineProduction) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = sim->GetParam()->Get<SimParam>();

  SetupInflammation(sim);
  sim->GetScheduler()->Simulate(1);

  // Check initial inflammation at cell position
  auto* rm = sim->GetResourceManager();
  auto* infl_grid = GetInflammationGrid(sim);
  real_t initial = infl_grid->GetValue({15, 15, 1.5});

  auto* cell = new ImmuneCell({15, 15, 1.5});
  cell->SetDiameter(3);
  cell->SetImmuneCellType(kNeutrophil);
  cell->AddBehavior(new NeutrophilBehavior());
  rm->AddAgent(cell);

  // Run a few steps so the neutrophil produces cytokines
  sim->GetScheduler()->Simulate(5);

  real_t after = infl_grid->GetValue({15, 15, 1.5});
  EXPECT_GT(after, initial);

  delete sim;
}

TEST(ImmuneResponseTest, NeutrophilSpawnCount) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  SanitizeForUnitTest(sp);
  sp->wound.trigger_step = 0;
  sp->immune.neutrophil_spawn_delay = 5;
  sp->immune.neutrophil_spawn_waves = 1;     // single wave for count test
  sp->immune.macrophage_spawn_delay = 9999;

  SetupInflammation(sim);

  auto* ir = new ImmuneResponse();
  OperationRegistry::GetInstance()->AddOperationImpl(
      "ImmuneResponse_test_neut", OpComputeTarget::kCpu, ir);
  sim->GetScheduler()->ScheduleOp(
      NewOperation("ImmuneResponse_test_neut"), OpType::kPreSchedule);

  // Before spawn delay: no agents
  sim->GetScheduler()->Simulate(5);
  EXPECT_EQ(sim->GetResourceManager()->GetNumAgents(), 0u);

  // Step 6: wound_age = 6 >= 5, neutrophils spawn
  sim->GetScheduler()->Simulate(1);

  // Expected count: round(2*pi*5 / 3) = round(10.47) = 10
  int expected = static_cast<int>(std::round(2.0 * M_PI * sp->wound.radius / 3.0));
  if (expected < 6) expected = 6;
  EXPECT_EQ(sim->GetResourceManager()->GetNumAgents(),
            static_cast<size_t>(expected));

  // Verify all are neutrophils
  sim->GetResourceManager()->ForEachAgent([](Agent* a) {
    auto* ic = dynamic_cast<ImmuneCell*>(a);
    ASSERT_NE(ic, nullptr);
    EXPECT_EQ(ic->GetImmuneCellType(), kNeutrophil);
    EXPECT_EQ(ic->GetState(), kM1Active);
  });

  delete sim;
}

TEST(ImmuneResponseTest, MacrophageContinuousRecruitment) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  SanitizeForUnitTest(sp);
  sp->wound.trigger_step = 0;
  sp->immune.neutrophil_spawn_delay = 9999;     // disable neutrophils
  sp->immune.macrophage_spawn_delay = 5;
  sp->immune.macrophage_spawn_threshold = 0.0;  // any inflammation triggers recruitment
  sp->immune.macrophage_spawn_rate = 1.0;       // high rate for deterministic test
  sp->immune.macrophage_apoptosis_rate = 0;     // disable death during test
  sp->immune.macrophage_lifespan = 9999;

  SetupInflammation(sim);

  // Initialize grids (boxes not created until first step)
  sim->GetScheduler()->Simulate(1);

  // Seed inflammation at wound center so recruitment has signal
  auto* rm = sim->GetResourceManager();
  auto* infl_grid = GetInflammationGrid(sim);
  Real3 center = {sp->wound.center_x, sp->wound.center_y, 1.5};
  size_t idx = infl_grid->GetBoxIndex(center);
  infl_grid->ChangeConcentrationBy(idx, 2.0);

  auto* ir = new ImmuneResponse();
  OperationRegistry::GetInstance()->AddOperationImpl(
      "ImmuneResponse_test_mac", OpComputeTarget::kCpu, ir);
  sim->GetScheduler()->ScheduleOp(
      NewOperation("ImmuneResponse_test_mac"), OpType::kPreSchedule);

  // Before spawn delay: no agents (steps 1-5, wound_age 1-4 < delay=5)
  sim->GetScheduler()->Simulate(4);
  EXPECT_EQ(rm->GetNumAgents(), 0u);

  // Run past delay with inflammation present: macrophages should appear
  sim->GetScheduler()->Simulate(20);
  EXPECT_GT(rm->GetNumAgents(), 0u);

  // Verify all are macrophages
  rm->ForEachAgent([](Agent* a) {
    auto* ic = dynamic_cast<ImmuneCell*>(a);
    ASSERT_NE(ic, nullptr);
    EXPECT_EQ(ic->GetImmuneCellType(), kMacrophage);
  });

  delete sim;
}

TEST(ImmuneResponseTest, NoRecruitmentWithoutInflammation) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  SanitizeForUnitTest(sp);
  sp->wound.trigger_step = 0;
  sp->immune.neutrophil_spawn_delay = 9999;
  sp->immune.macrophage_spawn_delay = 1;
  sp->immune.macrophage_spawn_threshold = 0.1;
  sp->immune.macrophage_spawn_rate = 1.0;

  SetupInflammation(sim);
  // No inflammation seeded -- field is zero

  auto* ir = new ImmuneResponse();
  OperationRegistry::GetInstance()->AddOperationImpl(
      "ImmuneResponse_test_mac_noinfl", OpComputeTarget::kCpu, ir);
  sim->GetScheduler()->ScheduleOp(
      NewOperation("ImmuneResponse_test_mac_noinfl"), OpType::kPreSchedule);

  sim->GetScheduler()->Simulate(50);
  EXPECT_EQ(sim->GetResourceManager()->GetNumAgents(), 0u);

  delete sim;
}

TEST(ImmuneResponseTest, MultiWaveSpawnCount) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  SanitizeForUnitTest(sp);
  sp->wound.trigger_step = 0;
  sp->immune.neutrophil_spawn_delay = 1;
  sp->immune.neutrophil_spawn_waves = 3;
  sp->immune.neutrophil_spawn_window = 6;   // waves at steps 1, 4, 7
  sp->immune.neutrophil_apoptosis_rate = 0;  // disable stochastic death
  sp->immune.neutrophil_lifespan = 9999;    // don't die during test
  sp->immune.macrophage_spawn_delay = 9999;

  SetupInflammation(sim);

  auto* ir = new ImmuneResponse();
  OperationRegistry::GetInstance()->AddOperationImpl(
      "ImmuneResponse_test_wave", OpComputeTarget::kCpu, ir);
  sim->GetScheduler()->ScheduleOp(
      NewOperation("ImmuneResponse_test_wave"), OpType::kPreSchedule);

  int total = static_cast<int>(
      std::round(2.0 * M_PI * sp->wound.radius / 3.0));
  if (total < 6) total = 6;
  int per_wave = total / 3;

  // After step 2: first wave spawned (wound_age=2 >= 1)
  sim->GetScheduler()->Simulate(2);
  EXPECT_EQ(sim->GetResourceManager()->GetNumAgents(),
            static_cast<size_t>(per_wave));

  // After step 5: second wave spawned (wound_age=5 >= 4)
  sim->GetScheduler()->Simulate(3);
  EXPECT_EQ(sim->GetResourceManager()->GetNumAgents(),
            static_cast<size_t>(per_wave * 2));

  // After step 8: third wave (wound_age=8 >= 7), last wave gets remainder
  sim->GetScheduler()->Simulate(3);
  EXPECT_EQ(sim->GetResourceManager()->GetNumAgents(),
            static_cast<size_t>(total));

  delete sim;
}

TEST(NeutrophilBehaviorTest, StochasticApoptosis) {
  auto* sim = CreateTestSim(TEST_NAME);
  sim->GetRandom()->SetSeed(42);  // fixed seed: deterministic across suite runs
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  SanitizeForUnitTest(sp);
  sp->wound.enabled = false;              // disable migration (not under test)
  sp->immune.neutrophil_lifespan = 9999;         // disable hard ceiling
  sp->immune.neutrophil_min_survival = 0;        // apoptosis from step 1
  sp->immune.neutrophil_apoptosis_rate = 0.01;   // ~1% per step

  SetupInflammation(sim);
  sim->GetScheduler()->Simulate(1);

  // Spawn 100 neutrophils
  for (int i = 0; i < 100; i++) {
    auto* cell = new ImmuneCell({15.0 + i * 0.01, 15, 1.5});
    cell->SetDiameter(3);
    cell->SetImmuneCellType(kNeutrophil);
    cell->AddBehavior(new NeutrophilBehavior());
    sim->GetResourceManager()->AddAgent(cell);
  }

  sim->GetScheduler()->Simulate(1);  // commit agents (apoptosis may fire)
  size_t initial = sim->GetResourceManager()->GetNumAgents();
  EXPECT_GE(initial, 95u);   // most survive the commit step
  EXPECT_LE(initial, 100u);

  // Run 100 steps: with 1% rate, expect ~63% dead (1 - e^{-1})
  sim->GetScheduler()->Simulate(100);
  size_t survivors = sim->GetResourceManager()->GetNumAgents();

  // Survivors should be roughly 37 (e^{-1} ~ 0.37), allow wide margin
  EXPECT_GT(survivors, 10u);   // not all dead
  EXPECT_LT(survivors, 70u);   // substantial deaths occurred

  delete sim;
}

TEST(EfferocytosisTest, M1EngulfsDyingNeutrophil) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  SanitizeForUnitTest(sp);
  sp->efferocytosis_enabled = true;
  sp->efferocytosis_radius = 2.0;  // must fit within env box length
  sp->efferocytosis_age_fraction = 0.5;
  sp->immune.neutrophil_lifespan = 10;
  sp->immune.macrophage_m1_duration = 9999;  // disable timer
  sp->immune.macrophage_lifespan = 100;

  SetupInflammation(sim);
  sim->GetScheduler()->Simulate(1);

  // Dying neutrophil (age 6 >= 10 * 0.5 = 5)
  auto* neutrophil = new ImmuneCell({15, 15, 1.5});
  neutrophil->SetDiameter(3);
  neutrophil->SetImmuneCellType(kNeutrophil);
  for (int i = 0; i < 6; i++) neutrophil->IncrementAge();
  neutrophil->AddBehavior(new NeutrophilBehavior());
  sim->GetResourceManager()->AddAgent(neutrophil);

  // M1 macrophage nearby
  auto* macrophage = new ImmuneCell({16, 15, 1.5});
  macrophage->SetDiameter(3);
  macrophage->SetImmuneCellType(kMacrophage);
  macrophage->SetState(kM1Active);
  macrophage->AddBehavior(new MacrophageBehavior());
  sim->GetResourceManager()->AddAgent(macrophage);

  // Step 1: macrophage efferocytosis sets neutrophil age to 999999.
  // Step 2: neutrophil's own behavior sees age > lifespan, removes itself.
  sim->GetScheduler()->Simulate(2);

  // Macrophage should have transitioned to M2 and neutrophil removed
  int n_remaining = 0;
  ImmuneState mac_state = kM1Active;
  sim->GetResourceManager()->ForEachAgent([&](Agent* a) {
    n_remaining++;
    auto* ic = dynamic_cast<ImmuneCell*>(a);
    if (ic && ic->GetImmuneCellType() == kMacrophage) {
      mac_state = ic->GetState();
    }
  });
  // Neutrophil engulfed -> removed, only macrophage remains
  EXPECT_EQ(n_remaining, 1);
  EXPECT_EQ(mac_state, kM2Resolving);

  delete sim;
}

TEST(EfferocytosisTest, YoungNeutrophilNotEngulfed) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  SanitizeForUnitTest(sp);
  sp->efferocytosis_enabled = true;
  sp->efferocytosis_radius = 2.0;  // must fit within env box length
  sp->efferocytosis_age_fraction = 0.8;
  sp->immune.neutrophil_lifespan = 100;
  sp->immune.macrophage_m1_duration = 9999;
  sp->immune.macrophage_lifespan = 100;

  SetupInflammation(sim);
  sim->GetScheduler()->Simulate(1);

  // Young neutrophil (age 1, threshold = 100 * 0.8 = 80)
  auto* neutrophil = new ImmuneCell({15, 15, 1.5});
  neutrophil->SetDiameter(3);
  neutrophil->SetImmuneCellType(kNeutrophil);
  neutrophil->AddBehavior(new NeutrophilBehavior());
  sim->GetResourceManager()->AddAgent(neutrophil);

  // M1 macrophage nearby
  auto* macrophage = new ImmuneCell({16, 15, 1.5});
  macrophage->SetDiameter(3);
  macrophage->SetImmuneCellType(kMacrophage);
  macrophage->SetState(kM1Active);
  macrophage->AddBehavior(new MacrophageBehavior());
  sim->GetResourceManager()->AddAgent(macrophage);

  sim->GetScheduler()->Simulate(1);

  // Both should survive, macrophage stays M1
  EXPECT_EQ(sim->GetResourceManager()->GetNumAgents(), 2u);
  sim->GetResourceManager()->ForEachAgent([](Agent* a) {
    auto* ic = dynamic_cast<ImmuneCell*>(a);
    if (ic && ic->GetImmuneCellType() == kMacrophage) {
      EXPECT_EQ(ic->GetState(), kM1Active);
    }
  });

  delete sim;
}

TEST(MacrophageBehaviorTest, M2ProducesTGFBeta) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  SanitizeForUnitTest(sp);
  sp->fibroblast.enabled = true;
  sp->immune.macrophage_lifespan = 100;

  SetupInflammation(sim);
  TGFBetaPDE tgfb(sp);
  tgfb.Init(sim);
  CollagenPDE col(sp);
  col.Init(sim);
  sim->GetScheduler()->Simulate(1);

  auto* rm = sim->GetResourceManager();
  auto* tgfb_grid = rm->GetDiffusionGrid(fields::kTGFBetaId);
  real_t before = tgfb_grid->GetValue({15, 15, 1.5});

  auto* cell = new ImmuneCell({15, 15, 1.5});
  cell->SetDiameter(3);
  cell->SetImmuneCellType(kMacrophage);
  cell->SetState(kM2Resolving);
  cell->AddBehavior(new MacrophageBehavior());
  rm->AddAgent(cell);

  sim->GetScheduler()->Simulate(3);

  real_t after = tgfb_grid->GetValue({15, 15, 1.5});
  EXPECT_GT(after, before);

  delete sim;
}

}  // namespace skibidy
}  // namespace bdm
