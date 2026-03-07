// Rheumatoid arthritis module tests (study-scoped).
// Tests IL-6/TNF-alpha field registration, flare dynamics, pannus boost,
// M1 prolongation, cartilage degradation, and treatment clearance.

#include "test_helpers.h"
#include "rheumatoid/tnf_alpha_pde.h"
#include "rheumatoid/il6_pde.h"
#include "rheumatoid/cartilage_pde.h"
#include "rheumatoid/source_hook.h"
#include "rheumatoid/post_hook.h"
#include "core/grid_registry.h"
#include "core/voxel_snapshot.h"
#include "core/signal_board.h"
#include "fibroblast/fibroblast.h"
#include "fibroblast/fibroblast_behavior.h"
#include "fibroblast/fibroblast_recruitment.h"
#include "immune/macrophage_behavior.h"
#include "immune/immune_response.h"

namespace bdm {
namespace skibidy {

// ===========================================================================
// RA field registration
// ===========================================================================
TEST(RAFieldTest, TNFAlphaAndIL6Registered) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  sp->ra.enabled = true;

  TNFAlphaPDE tnf(sp);
  tnf.Init(sim);
  IL6PDE il6(sp);
  il6.Init(sim);
  CartilagePDE cart(sp);
  cart.Init(sim);

  auto* rm = sim->GetResourceManager();
  EXPECT_NE(rm->GetDiffusionGrid(fields::kTNFAlphaId), nullptr);
  EXPECT_NE(rm->GetDiffusionGrid(fields::kIL6Id), nullptr);
  EXPECT_NE(rm->GetDiffusionGrid(fields::kCartilageId), nullptr);

  delete sim;
}

TEST(RAFieldTest, CartilageInitializedInDermis) {
  Param::RegisterParamGroup(new SimParam());
  auto set_param = [](Param* param) {
    param->bound_space = Param::BoundSpaceMode::kClosed;
    param->min_bound = -10;
    param->max_bound = 50;
    param->simulation_time_step = 0.1;
    param->export_visualization = false;
  };
  auto* sim = new Simulation(TEST_NAME, set_param);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  sp->ra.enabled = true;
  sp->ra.cartilage_basal = 1.0;
  sp->grid_resolution_structural = 0;

  CartilagePDE cart(sp);
  cart.Init(sim);
  sim->GetScheduler()->Simulate(1);

  auto* grid = sim->GetResourceManager()->GetDiffusionGrid(fields::kCartilageId);
  GridContext ctx(grid, sp);

  int dermal_count = 0;
  real_t dermal_sum = 0;
  int epi_count = 0;
  real_t epi_sum = 0;
  for (size_t idx = 0; idx < ctx.n; idx++) {
    real_t z = ctx.Z(idx);
    real_t val = grid->GetConcentration(idx);
    if (z < 0) {
      dermal_count++;
      dermal_sum += val;
    } else {
      epi_count++;
      epi_sum += val;
    }
  }

  // Dermal voxels should have cartilage = basal
  ASSERT_GT(dermal_count, 0);
  EXPECT_NEAR(dermal_sum / dermal_count, 1.0, 0.01);

  // Epidermal voxels should have no cartilage
  if (epi_count > 0) {
    EXPECT_NEAR(epi_sum / epi_count, 0.0, 0.01);
  }

  delete sim;
}

// ===========================================================================
// Flare dynamics (sigmoid activation)
// ===========================================================================
TEST(RAFlareTest, SigmoidRampsUp) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  SanitizeForUnitTest(sp);
  sp->ra.enabled = true;
  sp->ra.flare_delay = 100;
  sp->ra.flare_steepness = 0.05;
  sp->wound.trigger_step = 10;

  SetupInflammation(sim);
  TNFAlphaPDE tnf(sp); tnf.Init(sim);
  IL6PDE il6(sp); il6.Init(sim);

  // Before flare onset: factor should be low
  sim->GetScheduler()->Simulate(50);  // step 50, wound_age = 40 << delay 100
  {
    GridRegistry reg;
    reg.Fill(sim);
    SignalBoard sig;
    RASourceHook hook;
    hook.Init(reg, sig);
    EXPECT_LT(hook.flare_factor, 0.1);  // well before delay
  }

  // After flare onset: factor should be high
  sim->GetScheduler()->Simulate(200);  // step 250, wound_age = 240 >> delay 100
  {
    GridRegistry reg;
    reg.Fill(sim);
    SignalBoard sig;
    RASourceHook hook;
    hook.Init(reg, sig);
    EXPECT_GT(hook.flare_factor, 0.9);  // well past delay
  }

  delete sim;
}

// ===========================================================================
// Pannus fibroblast boost
// ===========================================================================
TEST(RAPannusTest, FibroblastCountBoosted) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  sp->wound.trigger_step = 0;
  sp->fibroblast.enabled = true;
  sp->fibroblast.spawn_delay = 5;
  sp->fibroblast.spawn_waves = 1;
  sp->ra.enabled = true;
  sp->ra.pannus_fibroblast_boost = 3.0;

  TGFBetaPDE tgfb(sp); tgfb.Init(sim);
  CollagenPDE col(sp); col.Init(sim);

  auto* fr = new FibroblastRecruitment();
  OperationRegistry::GetInstance()->AddOperationImpl(
      "FibroblastRecruitment_pannus", OpComputeTarget::kCpu, fr);
  sim->GetScheduler()->ScheduleOp(
      NewOperation("FibroblastRecruitment_pannus"), OpType::kPreSchedule);

  // Trigger spawn: step 6, wound_age = 6 >= delay 5
  sim->GetScheduler()->Simulate(6);

  // Expected count without boost
  real_t outer_r = sp->wound.radius + sp->dermal_fibroblast_margin;
  int normal = static_cast<int>(
      std::ceil(2.0 * M_PI * outer_r / sp->fibroblast.diameter *
                sp->fibroblast.density_factor));
  if (normal < 6) normal = 6;
  // With 3x boost
  int expected = static_cast<int>(normal * 3.0);
  if (expected < 6) expected = 6;

  EXPECT_EQ(sim->GetResourceManager()->GetNumAgents(),
            static_cast<size_t>(expected));

  delete sim;
}

// ===========================================================================
// M1 prolongation
// ===========================================================================
TEST(RAM1Test, ProlongedM1Duration) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  SanitizeForUnitTest(sp);
  sp->immune.macrophage_m1_duration = 5;
  sp->immune.macrophage_lifespan = 100;
  sp->m1_transition_threshold = 0;  // disable cytokine path
  sp->ra.enabled = true;
  sp->ra.m1_prolongation = 3.0;  // 3x M1 duration

  SetupInflammation(sim);
  sim->GetScheduler()->Simulate(1);

  auto* cell = new ImmuneCell({15, 15, 1.5});
  cell->SetDiameter(3);
  cell->SetImmuneCellType(kMacrophage);
  cell->SetState(kM1Active);
  cell->AddBehavior(new MacrophageBehavior());
  sim->GetResourceManager()->AddAgent(cell);

  // Normal M1 duration=5, RA factor=3 -> effective=15
  // After 8 steps: still M1 (8 < 15)
  sim->GetScheduler()->Simulate(8);
  ImmuneState state = kM1Active;
  sim->GetResourceManager()->ForEachAgent([&](Agent* a) {
    auto* ic = dynamic_cast<ImmuneCell*>(a);
    if (ic) state = ic->GetState();
  });
  EXPECT_EQ(state, kM1Active);

  // After 10 more (total 18): should transition (18 > 15)
  sim->GetScheduler()->Simulate(10);
  sim->GetResourceManager()->ForEachAgent([&](Agent* a) {
    auto* ic = dynamic_cast<ImmuneCell*>(a);
    if (ic) state = ic->GetState();
  });
  EXPECT_EQ(state, kM2Resolving);

  delete sim;
}

// ===========================================================================
// Anti-TNF and anti-IL6R treatment clearance
// ===========================================================================
TEST(RATreatmentTest, AntiTNFReducesTNF) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  SanitizeForUnitTest(sp);
  sp->ra.enabled = true;
  sp->ra.autoimmune_source = 0.0;       // no production
  sp->ra.il6_autoimmune_source = 0.0;
  sp->ra.tnf_m1_rate = 0.0;
  sp->ra.tnf_neutrophil_rate = 0.0;
  sp->ra.il6_m1_rate = 0.0;
  sp->ra.anti_tnf_clearance = 0.5;      // aggressive clearance
  sp->wound.trigger_step = 0;

  SetupInflammation(sim);
  TNFAlphaPDE tnf_pde(sp); tnf_pde.Init(sim);
  IL6PDE il6_pde(sp); il6_pde.Init(sim);
  sim->GetScheduler()->Simulate(1);

  // Seed TNF-alpha manually
  auto* rm = sim->GetResourceManager();
  auto* tnf_grid = rm->GetDiffusionGrid(fields::kTNFAlphaId);
  size_t test_idx = 100;
  tnf_grid->ChangeConcentrationBy(test_idx, 1.0);
  real_t before = tnf_grid->GetConcentration(test_idx);
  EXPECT_NEAR(before, 1.0, 0.01);

  // Apply source hook (which runs clearance)
  GridRegistry reg;
  reg.Fill(sim);
  SignalBoard sig;
  RASourceHook hook;
  hook.Init(reg, sig);

  // Create a wound voxel snapshot
  VoxelSnapshot snap;
  snap.idx = test_idx;
  snap.in_wound = true;
  snap.post_wound = true;
  snap.wound_age = 500;
  snap.stratum = -1.0;

  hook.ApplyDermal(snap, sig);

  real_t after = tnf_grid->GetConcentration(test_idx);
  EXPECT_LT(after, before);  // clearance reduced it

  delete sim;
}

}  // namespace skibidy
}  // namespace bdm
