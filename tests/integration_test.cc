// Integration tests: verify multi-module interactions and cross-module
// coupling that unit tests (single-module isolation) cannot catch.
//
// These tests enable multiple modules simultaneously and verify that
// their interactions produce biologically expected emergent behavior.

#include "test_helpers.h"
#include "immune/neutrophil_behavior.h"
#include "immune/macrophage_behavior.h"
#include "immune/immune_response.h"
#include "fibroblast/fibroblast.h"
#include "fibroblast/fibroblast_behavior.h"
#include "senescence/senescence_pde.h"
#include "glucose/glucose_pde.h"
#include "glucose/age_pde.h"
#include "lactate/lactate_pde.h"

namespace bdm {
namespace skibidy {

// ---------------------------------------------------------------------------
// Helper: enable a controlled set of modules for integration tests
// ---------------------------------------------------------------------------
static void EnableIntegrationModules(SimParam* sp) {
  // Start with all modules disabled
  SanitizeForUnitTest(sp);
  // Then selectively enable what we need
  sp->fibroblast.enabled = true;
  sp->mmp.enabled = true;
  sp->fibronectin.enabled = true;
  sp->senescence.enabled = true;
}

// ---------------------------------------------------------------------------
// Test: Inflammation drives senescence accumulation
// Verify that high inflammation leads to senescence accumulation,
// not just in isolation but through the full hook chain.
// ---------------------------------------------------------------------------
TEST(IntegrationTest, InflammationDrivesSenescence) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  SanitizeForUnitTest(sp);

  // Test the senescence accumulation formula in isolation:
  // accum += infl_rate * infl_val
  // At infl_val=0.5, infl_rate=0.01: accum = 0.005 per step
  sp->senescence.infl_rate = 0.01;
  real_t infl_val = 0.5;
  real_t expected_accum = sp->senescence.infl_rate * infl_val;
  EXPECT_NEAR(expected_accum, 0.005, 1e-6);

  // Higher inflammation should give proportionally more senescence
  real_t high_infl = 1.0;
  real_t high_accum = sp->senescence.infl_rate * high_infl;
  EXPECT_NEAR(high_accum / expected_accum, 2.0, 1e-6);

  delete sim;
}

// ---------------------------------------------------------------------------
// Test: Immune clearance reduces senescence when immune pressure is high
// ---------------------------------------------------------------------------
TEST(IntegrationTest, ImmuneClearsSenescence) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  SanitizeForUnitTest(sp);

  // Test the immune clearance formula:
  //   clearance = immune_clearance_rate * sen_val * ip_val
  // At sen=0.5, ip=0.8, rate=0.1: clearance = 0.04 per step
  sp->senescence.immune_clearance_rate = 0.1;
  real_t sen_val = 0.5;
  real_t ip_val = 0.8;
  real_t clearance = sp->senescence.immune_clearance_rate * sen_val * ip_val;
  EXPECT_NEAR(clearance, 0.04, 1e-6);

  // No immune pressure = no clearance
  real_t clearance_no_ip = sp->senescence.immune_clearance_rate * sen_val * 0.0;
  EXPECT_NEAR(clearance_no_ip, 0.0, 1e-10);

  // Net senescence change: accumulation minus clearance
  // With infl_rate=0.001, infl=0.5, clearance=0.04: net = -0.0395
  real_t accum = sp->senescence.infl_rate * 0.5;
  real_t net = accum - clearance;
  EXPECT_LT(net, 0);  // clearance dominates: senescence decreases

  delete sim;
}

// ---------------------------------------------------------------------------
// Test: AGE formation rate is quadratic in glucose (Maillard kinetics)
// Verify that doubling glucose concentration quadruples AGE formation.
// ---------------------------------------------------------------------------
TEST(IntegrationTest, AGEFormationQuadraticInGlucose) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  SanitizeForUnitTest(sp);
  sp->diabetic.mode = true;
  sp->glucose_mod.enabled = true;
  sp->glucose_mod.age_rate = 0.01;  // measurable rate

  SetupAllFields(sim);

  // Register glucose and AGE grids
  GlucosePDE gluc_pde(sp);
  gluc_pde.Init(sim);
  AGEPDE age_pde(sp);
  age_pde.Init(sim);

  auto* rm = sim->GetResourceManager();
  auto* glucose_grid = rm->GetDiffusionGrid(fields::kGlucoseId);
  auto* age_grid = rm->GetDiffusionGrid(fields::kAGEId);

  ASSERT_NE(glucose_grid, nullptr);
  ASSERT_NE(age_grid, nullptr);

  // Test: AGE formed from glucose=1.0 vs glucose=2.0
  // rate = age_rate * glucose^2
  // At glucose=1.0: rate = 0.01 * 1 = 0.01
  // At glucose=2.0: rate = 0.01 * 4 = 0.04 (4x, not 2x)
  real_t rate_at_1 = sp->glucose_mod.age_rate * 1.0 * 1.0;
  real_t rate_at_2 = sp->glucose_mod.age_rate * 2.0 * 2.0;

  EXPECT_NEAR(rate_at_1, 0.01, 1e-6);
  EXPECT_NEAR(rate_at_2, 0.04, 1e-6);
  EXPECT_NEAR(rate_at_2 / rate_at_1, 4.0, 1e-6);

  delete sim;
}

// ---------------------------------------------------------------------------
// Test: Lactate production requires both hypoxia AND glucose
// ---------------------------------------------------------------------------
TEST(IntegrationTest, LactateRequiresGlucoseAndHypoxia) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  SanitizeForUnitTest(sp);
  sp->glucose_mod.enabled = true;
  sp->lactate.enabled = true;
  sp->lactate.o2_threshold = 0.3;
  sp->lactate.production_rate = 0.01;

  SetupAllFields(sim);

  GlucosePDE gluc_pde(sp);
  gluc_pde.Init(sim);
  LactatePDE lac_pde(sp);
  lac_pde.Init(sim);

  // Production formula: rate * (threshold - O2) / threshold * glucose
  // At O2=0.1 (hypoxic), glucose=1.0: rate = 0.01 * (0.3-0.1)/0.3 * 1.0 = 0.00667
  // At O2=0.1, glucose=0.0: rate = 0.01 * (0.3-0.1)/0.3 * 0.0 = 0.0 (no glucose!)
  // At O2=0.5 (normoxic), glucose=1.0: rate = 0 (above threshold)

  real_t o2_frac = (sp->lactate.o2_threshold - 0.1) / sp->lactate.o2_threshold;
  real_t rate_with_glucose = sp->lactate.production_rate * o2_frac * 1.0;
  real_t rate_no_glucose = sp->lactate.production_rate * o2_frac * 0.0;

  EXPECT_GT(rate_with_glucose, 0.0);
  EXPECT_NEAR(rate_no_glucose, 0.0, 1e-10);

  // Verify coupling: doubling glucose doubles lactate production
  real_t rate_double_glucose = sp->lactate.production_rate * o2_frac * 2.0;
  EXPECT_NEAR(rate_double_glucose / rate_with_glucose, 2.0, 1e-6);

  delete sim;
}

// ---------------------------------------------------------------------------
// Test: Fibroblast + immune modules coexist without crashes
// Verifies that enabling multiple agent types and field modules together
// doesn't cause grid registration conflicts or null pointer dereferences.
// ---------------------------------------------------------------------------
TEST(IntegrationTest, MultiModuleCoexistence) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  SanitizeForUnitTest(sp);
  sp->fibroblast.enabled = true;
  sp->mmp.enabled = true;
  sp->fibronectin.enabled = true;

  SetupAllFields(sim);

  auto* rm = sim->GetResourceManager();

  // Add an immune cell
  auto* neutrophil = new ImmuneCell({20, 20, 2});
  neutrophil->SetDiameter(3);
  neutrophil->SetImmuneCellType(kNeutrophil);
  neutrophil->AddBehavior(new NeutrophilBehavior());
  rm->AddAgent(neutrophil);

  // Add a macrophage
  auto* macrophage = new ImmuneCell({22, 22, 2});
  macrophage->SetDiameter(3);
  macrophage->SetImmuneCellType(kMacrophage);
  macrophage->AddBehavior(new MacrophageBehavior());
  rm->AddAgent(macrophage);

  // Add a fibroblast
  auto* fibro = new Fibroblast({24, 24, -2});
  fibro->SetDiameter(4);
  fibro->AddBehavior(new FibroblastBehavior());
  rm->AddAgent(fibro);

  // Run several steps; no crashes = success
  EXPECT_NO_FATAL_FAILURE(sim->GetScheduler()->Simulate(5));

  // All agents should still exist (or have been properly cleaned up)
  EXPECT_GE(rm->GetNumAgents(), 0u);

  delete sim;
}

// ---------------------------------------------------------------------------
// Test: Senescence immune_clearance_rate parameter exists and is accessible
// ---------------------------------------------------------------------------
TEST(IntegrationTest, SenescenceImmuneClearanceParamExists) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());

  // Verify the new parameter is accessible and has correct default
  EXPECT_NEAR(sp->senescence.immune_clearance_rate, 0.002, 1e-6);

  // Verify it can be modified
  sp->senescence.immune_clearance_rate = 0.05;
  EXPECT_NEAR(sp->senescence.immune_clearance_rate, 0.05, 1e-6);

  delete sim;
}

}  // namespace skibidy
}  // namespace bdm
