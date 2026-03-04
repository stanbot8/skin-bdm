#include "test_helpers.h"
#include "fibroblast/fibroblast.h"
#include "fibroblast/fibroblast_behavior.h"
#include "fibroblast/fibroblast_recruitment.h"
#include "elastin/elastin_pde.h"
#include "hyaluronan/hyaluronan_pde.h"

namespace bdm {
namespace skibidy {

TEST(FibroblastTest, TypeAndStateDefaults) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* cell = new Fibroblast({15, 15, 2});
  cell->SetDiameter(4);

  EXPECT_EQ(cell->GetFibroblastState(), kFibroQuiescent);
  EXPECT_EQ(cell->GetAge(), 0);
  EXPECT_EQ(cell->GetStateAge(), 0);

  delete cell;
  delete sim;
}

TEST(FibroblastTest, StateTransitionResetsStateAge) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* cell = new Fibroblast({15, 15, 2});

  cell->IncrementStateAge();
  cell->IncrementStateAge();
  cell->IncrementStateAge();
  EXPECT_EQ(cell->GetStateAge(), 3);

  cell->SetFibroblastState(kFibroActivated);
  EXPECT_EQ(cell->GetStateAge(), 0);

  cell->IncrementStateAge();
  cell->SetFibroblastState(kMyofibroblast);
  EXPECT_EQ(cell->GetStateAge(), 0);

  delete cell;
  delete sim;
}

TEST(FibroblastBehaviorTest, ActivationByTGFBeta) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  SanitizeForUnitTest(sp);
  sp->fibroblast.enabled = true;
  sp->fibroblast.activation_threshold = 0.1;
  sp->fibroblast.activation_delay = 9999;  // disable timeout

  TGFBetaPDE tgfb(sp);
  tgfb.Init(sim);
  CollagenPDE col(sp);
  col.Init(sim);
  sim->GetScheduler()->Simulate(1);

  // Seed TGF-beta above threshold at cell position
  auto* rm = sim->GetResourceManager();
  auto* tgfb_grid = rm->GetDiffusionGrid(fields::kTGFBetaId);
  Real3 pos = {15, 15, 2};
  size_t idx = tgfb_grid->GetBoxIndex(pos);
  tgfb_grid->ChangeConcentrationBy(idx, 0.5);

  auto* cell = new Fibroblast({15, 15, 2});
  cell->SetDiameter(4);
  cell->SetFibroblastState(kFibroQuiescent);

  FibroblastBehavior beh;
  beh.Run(cell);

  EXPECT_EQ(cell->GetFibroblastState(), kFibroActivated);

  delete cell;
  delete sim;
}

TEST(FibroblastBehaviorTest, ActivationByTimeout) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  SanitizeForUnitTest(sp);
  sp->fibroblast.enabled = true;
  sp->fibroblast.activation_threshold = 999;  // very high, won't trigger
  sp->fibroblast.activation_delay = 3;

  TGFBetaPDE tgfb(sp);
  tgfb.Init(sim);
  CollagenPDE col(sp);
  col.Init(sim);
  sim->GetScheduler()->Simulate(1);

  auto* cell = new Fibroblast({15, 15, 2});
  cell->SetDiameter(4);
  cell->SetFibroblastState(kFibroQuiescent);
  cell->AddBehavior(new FibroblastBehavior());
  sim->GetResourceManager()->AddAgent(cell);

  // After 3 steps: state_age=3, check is > 3, still quiescent
  sim->GetScheduler()->Simulate(3);
  FibroblastState state_at_3 = kFibroQuiescent;
  sim->GetResourceManager()->ForEachAgent([&](Agent* a) {
    auto* f = dynamic_cast<Fibroblast*>(a);
    if (f) state_at_3 = f->GetFibroblastState();
  });
  EXPECT_EQ(state_at_3, kFibroQuiescent);

  // Step 4: state_age=4 > 3, activates
  sim->GetScheduler()->Simulate(1);
  FibroblastState state_at_4 = kFibroQuiescent;
  sim->GetResourceManager()->ForEachAgent([&](Agent* a) {
    auto* f = dynamic_cast<Fibroblast*>(a);
    if (f) state_at_4 = f->GetFibroblastState();
  });
  EXPECT_EQ(state_at_4, kFibroActivated);

  delete sim;
}

TEST(FibroblastBehaviorTest, MyofibroblastTransition) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  SanitizeForUnitTest(sp);
  sp->fibroblast.enabled = true;
  sp->fibroblast.myofibroblast_threshold = 0.2;
  sp->fibroblast.myofibroblast_delay = 2;
  sp->fibroblast.lifespan = 9999;

  TGFBetaPDE tgfb(sp);
  tgfb.Init(sim);
  CollagenPDE col(sp);
  col.Init(sim);
  sim->GetScheduler()->Simulate(1);

  // Seed TGF-beta above myofibroblast threshold
  auto* rm = sim->GetResourceManager();
  auto* tgfb_grid = rm->GetDiffusionGrid(fields::kTGFBetaId);
  size_t idx = tgfb_grid->GetBoxIndex(Real3{15, 15, 2});
  tgfb_grid->ChangeConcentrationBy(idx, 1.0);

  auto* cell = new Fibroblast({15, 15, 2});
  cell->SetDiameter(4);
  cell->SetFibroblastState(kFibroActivated);
  cell->AddBehavior(new FibroblastBehavior());
  rm->AddAgent(cell);

  // After 2 steps: state_age=2, check is > 2, still activated
  sim->GetScheduler()->Simulate(2);
  FibroblastState state = kFibroActivated;
  rm->ForEachAgent([&](Agent* a) {
    auto* f = dynamic_cast<Fibroblast*>(a);
    if (f) state = f->GetFibroblastState();
  });
  EXPECT_EQ(state, kFibroActivated);

  // Re-seed TGF-beta (may have decayed)
  tgfb_grid->ChangeConcentrationBy(idx, 1.0);

  // Step 3: state_age=3 > 2, transitions to myofibroblast
  sim->GetScheduler()->Simulate(1);
  rm->ForEachAgent([&](Agent* a) {
    auto* f = dynamic_cast<Fibroblast*>(a);
    if (f) state = f->GetFibroblastState();
  });
  EXPECT_EQ(state, kMyofibroblast);

  delete sim;
}

TEST(FibroblastBehaviorTest, CollagenDeposition) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  SanitizeForUnitTest(sp);
  sp->fibroblast.enabled = true;
  sp->fibroblast.collagen_deposition_rate = 0.1;
  sp->fibroblast.tgfb_rate = 0;  // don't write TGF-beta (isolate collagen test)
  sp->fibroblast.lifespan = 9999;
  sp->fibroblast.apoptosis_threshold = 0;  // prevent removal

  TGFBetaPDE tgfb(sp);
  tgfb.Init(sim);
  CollagenPDE col(sp);
  col.Init(sim);
  sim->GetScheduler()->Simulate(1);

  // Seed TGF-beta
  auto* rm = sim->GetResourceManager();
  auto* tgfb_grid = rm->GetDiffusionGrid(fields::kTGFBetaId);
  Real3 pos = {15, 15, 2};
  size_t idx = tgfb_grid->GetBoxIndex(pos);
  tgfb_grid->ChangeConcentrationBy(idx, 2.0);

  auto* cell = new Fibroblast(pos);
  cell->SetDiameter(4);
  cell->SetFibroblastState(kMyofibroblast);

  FibroblastBehavior beh;
  beh.Run(cell);

  auto* col_grid = rm->GetDiffusionGrid(fields::kCollagenId);
  size_t col_idx = col_grid->GetBoxIndex(pos);
  real_t col_val = col_grid->GetConcentration(col_idx);

  // collagen += rate = 0.1 (constant deposit, not TGF-beta-dependent)
  EXPECT_NEAR(col_val, 0.1, 0.01);

  delete cell;
  delete sim;
}

TEST(FibroblastBehaviorTest, TGFBetaProduction) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  SanitizeForUnitTest(sp);
  sp->fibroblast.enabled = true;
  sp->fibroblast.tgfb_rate = 0.01;
  sp->fibroblast.lifespan = 9999;
  sp->fibroblast.apoptosis_threshold = 0;

  TGFBetaPDE tgfb(sp);
  tgfb.Init(sim);
  CollagenPDE col(sp);
  col.Init(sim);
  sim->GetScheduler()->Simulate(1);

  auto* rm = sim->GetResourceManager();
  auto* tgfb_grid = rm->GetDiffusionGrid(fields::kTGFBetaId);
  Real3 pos = {15, 15, 2};
  real_t before = tgfb_grid->GetValue(pos);

  auto* cell = new Fibroblast(pos);
  cell->SetDiameter(4);
  cell->SetFibroblastState(kMyofibroblast);

  FibroblastBehavior beh;
  beh.Run(cell);

  real_t after = tgfb_grid->GetValue(pos);
  EXPECT_GT(after, before);

  delete cell;
  delete sim;
}

TEST(FibroblastBehaviorTest, ApoptosisOnLowTGFBeta) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  SanitizeForUnitTest(sp);
  sp->fibroblast.enabled = true;
  sp->fibroblast.apoptosis_threshold = 0.5;
  sp->fibroblast.min_lifespan = 3;
  sp->fibroblast.lifespan = 9999;

  TGFBetaPDE tgfb(sp);
  tgfb.Init(sim);
  CollagenPDE col(sp);
  col.Init(sim);
  sim->GetScheduler()->Simulate(1);

  // No TGF-beta seeded (all zero, below threshold 0.5)
  auto* cell = new Fibroblast({15, 15, 2});
  cell->SetDiameter(4);
  cell->SetFibroblastState(kMyofibroblast);
  cell->AddBehavior(new FibroblastBehavior());
  sim->GetResourceManager()->AddAgent(cell);

  // 3 steps: age=3, check is > 3 for min_lifespan, still alive
  sim->GetScheduler()->Simulate(3);
  EXPECT_EQ(sim->GetResourceManager()->GetNumAgents(), 1u);

  // Step 4: age=4 > 3, TGF-beta < 0.5, apoptosis
  sim->GetScheduler()->Simulate(1);
  EXPECT_EQ(sim->GetResourceManager()->GetNumAgents(), 0u);

  delete sim;
}

TEST(FibroblastBehaviorTest, LifespanHardLimit) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  SanitizeForUnitTest(sp);
  sp->fibroblast.enabled = true;
  sp->fibroblast.lifespan = 5;
  sp->fibroblast.apoptosis_threshold = 0;  // prevent apoptosis path

  TGFBetaPDE tgfb(sp);
  tgfb.Init(sim);
  CollagenPDE col(sp);
  col.Init(sim);
  sim->GetScheduler()->Simulate(1);

  auto* cell = new Fibroblast({15, 15, 2});
  cell->SetDiameter(4);
  cell->SetFibroblastState(kFibroQuiescent);
  cell->AddBehavior(new FibroblastBehavior());
  sim->GetResourceManager()->AddAgent(cell);

  // 5 steps: age=5, check is > 5, still alive
  sim->GetScheduler()->Simulate(5);
  EXPECT_EQ(sim->GetResourceManager()->GetNumAgents(), 1u);

  // Step 6: age=6 > 5, removed
  sim->GetScheduler()->Simulate(1);
  EXPECT_EQ(sim->GetResourceManager()->GetNumAgents(), 0u);

  delete sim;
}

TEST(FibroblastBehaviorTest, MigrationTowardCenter) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  SanitizeForUnitTest(sp);
  sp->fibroblast.enabled = true;
  sp->fibroblast.lifespan = 9999;
  sp->fibroblast.activation_threshold = 0;  // auto-activate immediately

  TGFBetaPDE tgfb(sp);
  tgfb.Init(sim);
  CollagenPDE col(sp);
  col.Init(sim);
  sim->GetScheduler()->Simulate(1);

  // Cell east of wound center (center = 25,25 in test domain)
  auto* cell = new Fibroblast({30, 25, 2});
  cell->SetDiameter(4);
  cell->SetFibroblastState(kFibroActivated);

  FibroblastBehavior beh;
  beh.Run(cell);

  auto tf = cell->GetTractorForce();
  EXPECT_LT(tf[0], 0);          // toward center (west)
  EXPECT_NEAR(tf[1], 0, 1e-6);  // no y offset

  delete cell;
  delete sim;
}

TEST(FibroblastRecruitmentTest, SpawnTimingAndCount) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  sp->wound.trigger_step = 0;
  sp->fibroblast.enabled = true;
  sp->fibroblast.spawn_delay = 5;
  sp->fibroblast.spawn_waves = 1;  // single wave for deterministic count

  TGFBetaPDE tgfb(sp);
  tgfb.Init(sim);
  CollagenPDE col(sp);
  col.Init(sim);

  auto* fr = new FibroblastRecruitment();
  OperationRegistry::GetInstance()->AddOperationImpl(
      "FibroblastRecruitment_test", OpComputeTarget::kCpu, fr);
  sim->GetScheduler()->ScheduleOp(
      NewOperation("FibroblastRecruitment_test"), OpType::kPreSchedule);

  // Before spawn delay: no agents
  sim->GetScheduler()->Simulate(5);
  EXPECT_EQ(sim->GetResourceManager()->GetNumAgents(), 0u);

  // Step 6: wound_age = 6 >= 5, fibroblasts spawn
  sim->GetScheduler()->Simulate(1);

  real_t outer_r = sp->wound.radius + sp->dermal_fibroblast_margin;
  int expected = static_cast<int>(
      std::ceil(2.0 * M_PI * outer_r / sp->fibroblast.diameter *
                sp->fibroblast.density_factor));
  if (expected < 6) expected = 6;
  EXPECT_EQ(sim->GetResourceManager()->GetNumAgents(),
            static_cast<size_t>(expected));

  // Verify all are fibroblasts in quiescent state
  sim->GetResourceManager()->ForEachAgent([](Agent* a) {
    auto* f = dynamic_cast<Fibroblast*>(a);
    ASSERT_NE(f, nullptr);
    EXPECT_EQ(f->GetFibroblastState(), kFibroQuiescent);
  });

  delete sim;
}

TEST(FibroblastBehaviorTest, TropoelastinProduction) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  SanitizeForUnitTest(sp);
  sp->fibroblast.enabled = true;
  sp->elastin.enabled = true;
  sp->elastin.production_rate = 0.01;
  sp->fibroblast.tgfb_rate = 0;
  sp->fibroblast.lifespan = 9999;
  sp->fibroblast.apoptosis_threshold = 0;

  TGFBetaPDE tgfb(sp);
  tgfb.Init(sim);
  CollagenPDE col(sp);
  col.Init(sim);
  ElastinPDE elastin;
  elastin.Init(sim);
  sim->GetScheduler()->Simulate(1);

  auto* rm = sim->GetResourceManager();
  Real3 pos = {15, 15, 2};
  auto* el_grid = rm->GetDiffusionGrid(fields::kElastinId);
  real_t before = el_grid->GetValue(pos);

  auto* cell = new Fibroblast(pos);
  cell->SetDiameter(4);
  cell->SetFibroblastState(kFibroActivated);

  FibroblastBehavior beh;
  beh.Run(cell);

  real_t after = el_grid->GetValue(pos);
  EXPECT_GT(after, before);

  delete cell;
  delete sim;
}

TEST(FibroblastBehaviorTest, HyaluronanSynthesis) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  SanitizeForUnitTest(sp);
  sp->fibroblast.enabled = true;
  sp->hyaluronan.enabled = true;
  sp->hyaluronan.production_rate = 0.01;
  sp->fibroblast.tgfb_rate = 0;
  sp->fibroblast.lifespan = 9999;
  sp->fibroblast.apoptosis_threshold = 0;

  TGFBetaPDE tgfb(sp);
  tgfb.Init(sim);
  CollagenPDE col(sp);
  col.Init(sim);
  HyaluronanPDE ha;
  ha.Init(sim);
  sim->GetScheduler()->Simulate(1);

  auto* rm = sim->GetResourceManager();
  Real3 pos = {15, 15, 2};
  auto* ha_grid = rm->GetDiffusionGrid(fields::kHyaluronanId);
  real_t before = ha_grid->GetValue(pos);

  auto* cell = new Fibroblast(pos);
  cell->SetDiameter(4);
  cell->SetFibroblastState(kFibroActivated);

  FibroblastBehavior beh;
  beh.Run(cell);

  real_t after = ha_grid->GetValue(pos);
  EXPECT_GT(after, before);

  delete cell;
  delete sim;
}

TEST(FibroblastBehaviorTest, QuiescentNoECMProduction) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  SanitizeForUnitTest(sp);
  sp->fibroblast.enabled = true;
  sp->elastin.enabled = true;
  sp->hyaluronan.enabled = true;
  sp->elastin.production_rate = 0.01;
  sp->hyaluronan.production_rate = 0.01;
  sp->fibroblast.tgfb_rate = 0;
  sp->fibroblast.lifespan = 9999;
  sp->fibroblast.apoptosis_threshold = 0;
  sp->fibroblast.activation_delay = 9999;
  sp->fibroblast.activation_threshold = 999;

  TGFBetaPDE tgfb(sp);
  tgfb.Init(sim);
  CollagenPDE col(sp);
  col.Init(sim);
  ElastinPDE elastin;
  elastin.Init(sim);
  HyaluronanPDE ha;
  ha.Init(sim);
  sim->GetScheduler()->Simulate(1);

  auto* rm = sim->GetResourceManager();
  Real3 pos = {15, 15, 2};
  auto* el_grid = rm->GetDiffusionGrid(fields::kElastinId);
  auto* ha_grid = rm->GetDiffusionGrid(fields::kHyaluronanId);
  real_t el_before = el_grid->GetValue(pos);
  real_t ha_before = ha_grid->GetValue(pos);

  auto* cell = new Fibroblast(pos);
  cell->SetDiameter(4);
  // Leave in quiescent state (default)

  FibroblastBehavior beh;
  beh.Run(cell);

  // Quiescent fibroblasts should not produce elastin or HA
  EXPECT_NEAR(el_grid->GetValue(pos), el_before, 1e-10);
  EXPECT_NEAR(ha_grid->GetValue(pos), ha_before, 1e-10);

  delete cell;
  delete sim;
}

}  // namespace skibidy
}  // namespace bdm
