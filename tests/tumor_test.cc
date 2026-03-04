#include "test_helpers.h"

#include "tumor/tumor_cell.h"
#include "tumor/tumor_behavior.h"
#include "tumor/tumor_initiation.h"
#include "tumor/tumor_pde.h"

namespace bdm {
namespace skibidy {

TEST(TumorCellTest, DefaultState) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* cell = new TumorCell({15, 15, 2});
  cell->SetDiameter(4);

  EXPECT_EQ(cell->GetCyclePhase(), kG1);
  EXPECT_DOUBLE_EQ(cell->GetPhaseElapsed(), 0.0);
  EXPECT_DOUBLE_EQ(cell->GetCellAge(), 0.0);
  EXPECT_EQ(cell->GetTumorType(), 0);

  delete cell;
  delete sim;
}

TEST(TumorCellTest, DaughterInheritance) {
  auto* sim = CreateTestSim(TEST_NAME);
  SetupAllFields(sim);
  sim->GetScheduler()->Simulate(1);

  auto* mother = new TumorCell({15, 15, 2});
  mother->SetDiameter(6);  // large enough to divide
  mother->SetTumorType(0);
  mother->SetCyclePhase(kS);
  mother->SetPhaseElapsed(5.0);
  mother->SetCellAge(10.0);
  mother->AddBehavior(new TumorBehavior());
  sim->GetResourceManager()->AddAgent(mother);

  // Force division via Divide() -- returns the daughter directly
  // (daughter is in the execution context, not yet committed to RM,
  // so ForEachAgent won't see it; use the return value instead)
  auto* daughter = dynamic_cast<TumorCell*>(mother->Divide({1, 0, 0}));

  ASSERT_NE(daughter, nullptr);
  EXPECT_EQ(daughter->GetCyclePhase(), kG1);
  EXPECT_DOUBLE_EQ(daughter->GetPhaseElapsed(), 0.0);
  EXPECT_DOUBLE_EQ(daughter->GetCellAge(), 0.0);
  EXPECT_EQ(daughter->GetTumorType(), 0);

  delete sim;
}

TEST(TumorBehaviorTest, CycleFasterThanNormal) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  sp->tumor.cycle_factor = 0.5;  // 2x faster

  SetupAllFields(sim);
  sim->GetScheduler()->Simulate(1);

  // Normal G1: transition probability = phase_elapsed / (g1_duration * cf)
  // With cf=0.5, effective G1 = 3.5h vs normal 7h -> faster transition
  auto* cell = new TumorCell({15, 15, 2});
  cell->SetDiameter(4);
  cell->SetCyclePhase(kG1);
  cell->SetPhaseElapsed(0);

  // Advance phase_elapsed by running behavior multiple times
  real_t dt = sim->GetParam()->simulation_time_step;
  TumorBehavior beh;
  int steps_to_transition = 0;
  for (int i = 0; i < 1000 && cell->GetCyclePhase() == kG1; i++) {
    beh.Run(cell);
    steps_to_transition++;
  }

  // With cf=0.5, effective g1=3.5h, and dt=0.1h, transition should happen
  // much sooner than normal g1=7h (70 steps). Allow stochastic range.
  EXPECT_LT(steps_to_transition, 100);

  delete cell;
  delete sim;
}

TEST(TumorBehaviorTest, DividesWithHighNeighborCount) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  sp->tumor.max_neighbors = 30;
  sp->tumor.ci_steepness = 0;      // hard threshold (deterministic)
  sp->tumor.apoptosis_rate = 0;    // no stochastic death
  sp->tumor.cycle_factor = 0.01;   // near-certain M-phase exit in one step

  SetupAllFields(sim);

  // Force environment box_length above the neighbor search radius.
  // TumorBehavior searches sqrt(2.25 * d^2) = 1.5*d = 9 for d=6.
  // BDM box_length = largest agent diameter; a sizing agent far from
  // the test area ensures box_length >= 20 without affecting counts.
  auto* sizer = new TumorCell({0.5, 0.5, 0.5});
  sizer->SetDiameter(20);
  sim->GetResourceManager()->AddAgent(sizer);
  sim->GetScheduler()->Simulate(1);

  // Tumor in M phase, ready to divide
  auto* cell = new TumorCell({15, 15, 2});
  cell->SetDiameter(6);  // above division_diameter
  cell->SetCyclePhase(kM);
  cell->SetPhaseElapsed(10);  // well past m_duration*cf
  cell->AddBehavior(new TumorBehavior());
  sim->GetResourceManager()->AddAgent(cell);

  // Add 20 neighbors (above normal max_neighbors=14, below tumor_max=30)
  for (int i = 0; i < 20; i++) {
    auto* neighbor = new TumorCell({15.0 + (i % 5) * 0.5,
                                     15.0 + (i / 5) * 0.5, 2});
    neighbor->SetDiameter(4);
    neighbor->AddBehavior(new TumorBehavior());
    sim->GetResourceManager()->AddAgent(neighbor);
  }

  size_t before = sim->GetResourceManager()->GetNumAgents();
  sim->GetScheduler()->Simulate(1);
  size_t after = sim->GetResourceManager()->GetNumAgents();

  // Should have divided (added at least 1 daughter)
  EXPECT_GT(after, before);

  delete sim;
}

TEST(TumorBehaviorTest, NoDivisionLimit) {
  auto* sim = CreateTestSim(TEST_NAME);
  SetupAllFields(sim);

  // Sizing agent: force box_length above tumor neighbor search radius
  auto* sizer = new TumorCell({0.5, 0.5, 0.5});
  sizer->SetDiameter(20);
  sim->GetResourceManager()->AddAgent(sizer);
  sim->GetScheduler()->Simulate(1);

  // TumorCell has no divisions_left_ concept -- it always divides
  auto* cell = new TumorCell({15, 15, 2});
  cell->SetDiameter(4);
  cell->AddBehavior(new TumorBehavior());
  sim->GetResourceManager()->AddAgent(cell);

  // Run 500 steps: tumor should proliferate
  sim->GetScheduler()->Simulate(500);

  EXPECT_GT(sim->GetResourceManager()->GetNumAgents(), 1u);

  delete sim;
}

TEST(TumorBehaviorTest, IgnoresStratum) {
  auto* sim = CreateTestSim(TEST_NAME);
  SetupAllFields(sim);
  sim->GetScheduler()->Simulate(1);

  // Place tumor cell high (spinous zone) -- should still cycle
  auto* cell = new TumorCell({15, 15, 10});
  cell->SetDiameter(4);
  cell->SetCyclePhase(kG1);
  cell->SetPhaseElapsed(0);

  TumorBehavior beh;
  beh.Run(cell);

  // Should still be in cell cycle (not G0)
  // Phase elapsed should have advanced (dt added)
  EXPECT_GT(cell->GetPhaseElapsed(), 0);
  // Could still be G1 or S (stochastic), but NOT G0
  EXPECT_NE(cell->GetCyclePhase(), kG0);

  delete cell;
  delete sim;
}

TEST(TumorBehaviorTest, HypoxiaSuppression) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  sp->oxygen_prolif_threshold = 0.3;

  // Initialize O2 field to zero (severe hypoxia)
  ModelInitializer::DefineSubstance(fields::kOxygenId, fields::kOxygen,
      sp->oxygen_diffusion, sp->oxygen_decay, 10);
  ModelInitializer::InitializeSubstance(fields::kOxygenId,
      [](real_t, real_t, real_t) { return 0.0; });  // zero O2
  sim->GetScheduler()->Simulate(1);

  auto* cell = new TumorCell({15, 15, 2});
  cell->SetDiameter(4);
  cell->SetCyclePhase(kG1);
  cell->SetPhaseElapsed(0);

  // Run many steps -- with zero O2, G1->S transition probability is 0
  TumorBehavior beh;
  for (int i = 0; i < 200; i++) {
    beh.Run(cell);
  }

  // Cell should still be in G1 (O2 suppresses transition completely)
  EXPECT_EQ(cell->GetCyclePhase(), kG1);

  delete cell;
  delete sim;
}

TEST(TumorInitiationTest, SpawnTimingAndCount) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  sp->tumor.enabled = true;
  sp->tumor.seed_time = 5;
  sp->tumor.seed_count = 8;

  SetupAllFields(sim);

  auto* ti = new TumorInitiation();
  OperationRegistry::GetInstance()->AddOperationImpl(
      "TumorInitiation_test", OpComputeTarget::kCpu, ti);
  sim->GetScheduler()->ScheduleOp(
      NewOperation("TumorInitiation_test"), OpType::kPreSchedule);

  // Before seed time: no agents
  sim->GetScheduler()->Simulate(5);
  EXPECT_EQ(sim->GetResourceManager()->GetNumAgents(), 0u);

  // Step 6: step = 5 >= seed_time=5, tumor seeds
  sim->GetScheduler()->Simulate(1);
  EXPECT_EQ(sim->GetResourceManager()->GetNumAgents(), 8u);

  // Verify all are TumorCell
  sim->GetResourceManager()->ForEachAgent([](Agent* a) {
    auto* tc = dynamic_cast<TumorCell*>(a);
    ASSERT_NE(tc, nullptr);
    EXPECT_EQ(tc->GetCyclePhase(), kG1);
  });

  delete sim;
}

}  // namespace skibidy
}  // namespace bdm
