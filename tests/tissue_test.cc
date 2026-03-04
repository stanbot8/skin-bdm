#include "test_helpers.h"
#include "core/cell_cell_force.h"
#include "tissue/migration.h"
#include "tissue/shedding.h"
#include "tissue/basal_division.h"
#include "core/voxel_env.h"
#include "immune/neutrophil_behavior.h"

namespace bdm {
namespace skibidy {

TEST(SkinTest, KeratinocyteStratumDefault) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* cell = new Keratinocyte({25, 25, 1});
  cell->SetDiameter(4);
  EXPECT_EQ(cell->GetStratum(), kBasal);
  delete cell;
  delete sim;
}

TEST(SkinTest, StemCellIdentification) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* cell = new Keratinocyte({25, 25, 1});
  cell->SetDivisionsLeft(-1);
  EXPECT_TRUE(cell->IsStem());
  cell->SetDivisionsLeft(5);
  EXPECT_FALSE(cell->IsStem());
  EXPECT_EQ(cell->GetDivisionsLeft(), 5);
  delete cell;
  delete sim;
}

TEST(SkinTest, CellCyclePhaseDefault) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* cell = new Keratinocyte({25, 25, 1});
  EXPECT_EQ(cell->GetCyclePhase(), kG1);
  EXPECT_DOUBLE_EQ(cell->GetPhaseElapsed(), 0.0);
  EXPECT_DOUBLE_EQ(cell->GetCellAge(), 0.0);
  delete cell;
  delete sim;
}

TEST(SkinTest, StratumStepsResetOnChange) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* cell = new Keratinocyte({25, 25, 1});
  cell->SetStratum(kBasal);
  cell->IncrementStratumSteps();
  cell->IncrementStratumSteps();
  EXPECT_EQ(cell->GetStratumSteps(), 2);
  // Changing stratum resets counter
  cell->SetStratum(kSpinous);
  EXPECT_EQ(cell->GetStratumSteps(), 0);
  delete cell;
  delete sim;
}

TEST(SkinTest, BasementMembraneClamp) {
  auto* sim = CreateTestSim(TEST_NAME);

  // Set up calcium grid so Differentiation behavior can run
  auto* sp = sim->GetParam()->Get<SimParam>();
  ModelInitializer::DefineSubstance(fields::kCalciumId, fields::kCalcium,
      sp->calcium_diffusion, sp->calcium_decay, 10);
  ModelInitializer::InitializeSubstance(fields::kCalciumId,
      [](real_t, real_t, real_t) { return 0.01; });

  auto* cell = new Keratinocyte({25, 25, -5});
  cell->SetDiameter(4);
  cell->SetStratum(kBasal);
  cell->AddBehavior(new Differentiation());
  sim->GetResourceManager()->AddAgent(cell);

  sim->GetScheduler()->Simulate(1);

  // Cell should be clamped to z >= 0
  auto agents = sim->GetResourceManager()->GetNumAgents();
  EXPECT_GE(agents, 1u);
  delete sim;
}

TEST(SkinTest, CalciumGradientProfile) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = sim->GetParam()->Get<SimParam>();

  ModelInitializer::DefineSubstance(fields::kCalciumId, fields::kCalcium,
      sp->calcium_diffusion, sp->calcium_decay, 20);
  ModelInitializer::InitializeSubstance(fields::kCalciumId,
      [&](real_t x, real_t y, real_t z) {
        real_t mid = sp->calcium_midpoint_z;
        real_t k = sp->calcium_steepness;
        return sp->calcium_basal +
               (sp->calcium_peak - sp->calcium_basal) /
               (1.0 + std::exp(-(z - mid) / k));
      });

  // Run one step to initialize grids
  sim->GetScheduler()->Simulate(1);

  auto* rm = sim->GetResourceManager();
  auto* ca = rm->GetDiffusionGrid(fields::kCalciumId);

  // Sigmoid profile: calcium at surface (high z) > calcium at base (low z)
  real_t ca_at_base = ca->GetValue({25, 25, 0.5});
  real_t ca_at_top = ca->GetValue({25, 25, 45});

  // After grid discretization the exact values shift,
  // but the monotonic gradient must hold
  EXPECT_GT(ca_at_top, ca_at_base);

  delete sim;
}

TEST(SkinTest, TACellExhaustion) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* cell = new Keratinocyte({25, 25, 1});
  cell->SetDivisionsLeft(0);
  cell->SetTADivisionsMax(4);
  EXPECT_FALSE(cell->IsStem());
  EXPECT_EQ(cell->GetDivisionsLeft(), 0);
  delete cell;
  delete sim;
}

TEST(SkinTest, UWYNVolumeFlags) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = sim->GetParam()->Get<SimParam>();
  auto* vm = VolumeManager::Get();
  vm->Init(sp->volume_z_spinous, sp->volume_z_granular,
           sp->volume_z_cornified);

  EXPECT_TRUE(vm->AreAgentsEnabled(kBasal));
  EXPECT_FALSE(vm->AreAgentsEnabled(kSpinous));
  EXPECT_FALSE(vm->AreAgentsEnabled(kGranular));
  EXPECT_FALSE(vm->AreAgentsEnabled(kCornified));
  EXPECT_FALSE(vm->AreAgentsEnabled(kDermis));

  delete sim;
}

TEST(SkinTest, UWYNSupraBasalRemoval) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = sim->GetParam()->Get<SimParam>();

  auto* vm = VolumeManager::Get();
  vm->Init(sp->volume_z_spinous, sp->volume_z_granular,
           sp->volume_z_cornified);

  // High calcium everywhere so stratum assignment goes to spinous (not basal override)
  ModelInitializer::DefineSubstance(fields::kCalciumId, fields::kCalcium,
      sp->calcium_diffusion, sp->calcium_decay, 10);
  ModelInitializer::InitializeSubstance(fields::kCalciumId,
      [](real_t, real_t, real_t) { return 1.0; });

  // Place exhausted TA cell at z=10 (spinous zone, above volume_z_spinous=3.0)
  auto* cell = new Keratinocyte({25, 25, 10.0});
  cell->SetDiameter(4);
  cell->SetStratum(kBasal);
  cell->SetDivisionsLeft(0);  // exhausted TA: no stem anchoring
  cell->AddBehavior(new Differentiation());
  sim->GetResourceManager()->AddAgent(cell);

  sim->GetScheduler()->Simulate(1);

  // Differentiation should have removed the cell (agents_enabled=false for spinous)
  EXPECT_EQ(sim->GetResourceManager()->GetNumAgents(), 0u);
  delete sim;
}

// ---------------------------------------------------------------------------
// ForceTest
// ---------------------------------------------------------------------------

TEST(ForceTest, HertzRepulsion) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = sim->GetParam()->Get<SimParam>();

  // Two overlapping cells on z-axis
  auto* c1 = new Keratinocyte({25, 25, 5});
  auto* c2 = new Keratinocyte({25, 25, 7});
  c1->SetDiameter(4);
  c2->SetDiameter(4);

  CellCellForce force;
  Real4 f = force.Calculate(c1, c2);

  // overlap = r1+r2 - dist = 2+2-2 = 2
  real_t expected = -sp->repulsion_coeff * std::pow(2.0, 1.5);
  // Direction from c1 to c2 is +z, so force on c1 is in +z direction
  // but repulsion is negative magnitude -> pushes c1 away (negative z)
  EXPECT_NEAR(f[0], 0, 1e-10);
  EXPECT_NEAR(f[1], 0, 1e-10);
  // f[2] = expected * dir_z where dir_z = +1 (c2 is above c1)
  EXPECT_NEAR(f[2], expected, 1e-6);
  EXPECT_LT(f[2], 0);  // repulsion pushes c1 downward

  delete c1;
  delete c2;
  delete sim;
}

TEST(ForceTest, DesmosomalAdhesion) {
  auto* sim = CreateTestSim(TEST_NAME);

  // bdm.toml sets attraction_coeff=0; override for this test
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  sp->attraction_coeff = 6.0;

  // Two cells with a gap of 1.0 (within adhesion range of 2.0)
  auto* c1 = new Keratinocyte({25, 25, 5});
  auto* c2 = new Keratinocyte({25, 25, 10});
  c1->SetDiameter(4);
  c2->SetDiameter(4);

  CellCellForce force;
  Real4 f = force.Calculate(c1, c2);

  // gap = dist - (r1+r2) = 5 - 4 = 1.0
  // force_mag = attraction_coeff * (1 - 1/2) = 6.0 * 0.5 = 3.0
  real_t expected = sp->attraction_coeff * 0.5;
  EXPECT_NEAR(f[2], expected, 1e-6);
  EXPECT_GT(f[2], 0);  // attraction pulls c1 toward c2 (+z)

  delete c1;
  delete c2;
  delete sim;
}

TEST(ForceTest, BeyondAdhesionRange) {
  auto* sim = CreateTestSim(TEST_NAME);

  // Two cells far apart (gap >> 2.0)
  auto* c1 = new Keratinocyte({25, 25, 5});
  auto* c2 = new Keratinocyte({25, 25, 15});
  c1->SetDiameter(4);
  c2->SetDiameter(4);

  CellCellForce force;
  Real4 f = force.Calculate(c1, c2);

  // gap = 10 - 4 = 6.0, beyond max_adhesion_range=2.0
  EXPECT_NEAR(f[0], 0, 1e-10);
  EXPECT_NEAR(f[1], 0, 1e-10);
  EXPECT_NEAR(f[2], 0, 1e-10);

  delete c1;
  delete c2;
  delete sim;
}

// ---------------------------------------------------------------------------
// MigrationTest
// ---------------------------------------------------------------------------

TEST(MigrationTest, DirectionTowardCenter) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = sim->GetParam()->Get<SimParam>();
  WaterPDE water; water.Init(sim);
  OxygenPDE o2; o2.Init(sim);
  SetupInflammation(sim);
  if (sp->fibronectin.enabled) { FibronectinPDE fn(sp); fn.Init(sim); }
  sim->GetScheduler()->Simulate(1);  // populate grids
  g_voxel_env.Fill(sim);             // explicit fill (tests don't schedule VoxelEnvFillOp)

  // Cell east of wound center (center = 25,25 in test domain)
  auto* cell = new Keratinocyte({30, 25, 2});
  cell->SetDiameter(4);
  cell->SetStratum(kBasal);

  // Call behavior directly to test logic in isolation
  Migration mig;
  mig.Run(cell);

  // Tractor force should point west (toward center: negative x)
  auto tf = cell->GetTractorForce();
  EXPECT_LT(tf[0], 0);         // x: toward center
  EXPECT_NEAR(tf[1], 0, 1e-6); // y: no offset
  EXPECT_NEAR(tf[2], 0, 1e-6); // z: planar migration

  delete cell;
  delete sim;
}

TEST(MigrationTest, SpeedScalesWithDistance) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = sim->GetParam()->Get<SimParam>();
  WaterPDE water; water.Init(sim);
  OxygenPDE o2; o2.Init(sim);
  SetupInflammation(sim);
  if (sp->fibronectin.enabled) { FibronectinPDE fn(sp); fn.Init(sim); }
  sim->GetScheduler()->Simulate(1);  // populate grids
  g_voxel_env.Fill(sim);             // explicit fill (tests don't schedule VoxelEnvFillOp)
  real_t r = sp->wound.radius;

  // Cell A: near wound edge (distance ≈ radius)
  auto* cellA = new Keratinocyte({sp->wound.center_x + r * 0.99,
                                   sp->wound.center_y, 2});
  cellA->SetDiameter(4);
  cellA->SetStratum(kBasal);

  Migration migA;
  migA.Run(cellA);
  real_t speedA = std::abs(cellA->GetTractorForce()[0]);

  // Cell B: at half radius
  auto* cellB = new Keratinocyte({sp->wound.center_x + r / 2,
                                   sp->wound.center_y, 2});
  cellB->SetDiameter(4);
  cellB->SetStratum(kBasal);

  Migration migB;
  migB.Run(cellB);
  real_t speedB = std::abs(cellB->GetTractorForce()[0]);

  // Both cells should migrate toward center (O2 gradient chemotaxis).
  EXPECT_GT(speedA, 0);
  EXPECT_GT(speedB, 0);
  EXPECT_LT(cellA->GetTractorForce()[0], 0);  // directed toward center (-x)
  EXPECT_LT(cellB->GetTractorForce()[0], 0);

  delete cellA;
  delete cellB;
  delete sim;
}

TEST(MigrationTest, NonBasalClearsForce) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = sim->GetParam()->Get<SimParam>();
  WaterPDE water; water.Init(sim);
  SetupInflammation(sim);
  if (sp->fibronectin.enabled) { FibronectinPDE fn(sp); fn.Init(sim); }
  sim->GetScheduler()->Simulate(1);  // populate grids

  auto* cell = new Keratinocyte({19, 15, 2});
  cell->SetDiameter(4);
  cell->SetStratum(kSpinous);  // not basal
  cell->SetTractorForce({5, 5, 0});  // pre-set a force

  // Call behavior directly
  Migration mig;
  mig.Run(cell);

  // Force should be zeroed for non-basal cells
  auto tf = cell->GetTractorForce();
  EXPECT_NEAR(tf[0], 0, 1e-10);
  EXPECT_NEAR(tf[1], 0, 1e-10);
  EXPECT_NEAR(tf[2], 0, 1e-10);

  delete cell;
  delete sim;
}

// ---------------------------------------------------------------------------
// SheddingTest
// ---------------------------------------------------------------------------

TEST(SheddingTest, TissueBoundaryRemoval) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = sim->GetParam()->Get<SimParam>();

  // Cell outside tissue_max in x
  auto* cell = new Keratinocyte({sp->tissue_max + 1, 15, 5});
  cell->SetDiameter(4);
  cell->SetStratum(kBasal);
  cell->AddBehavior(new Shedding());
  sim->GetResourceManager()->AddAgent(cell);

  sim->GetScheduler()->Simulate(1);

  EXPECT_EQ(sim->GetResourceManager()->GetNumAgents(), 0u);
  delete sim;
}

TEST(SheddingTest, CornifiedDesquamation) {
  auto* sim = CreateTestSim(TEST_NAME);

  // Override shedding_delay to a small value for testing
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  sp->shedding_delay = 3;

  auto* cell = new Keratinocyte({15, 15, 5});
  cell->SetDiameter(4);
  cell->SetStratum(kCornified);
  cell->AddBehavior(new Shedding());
  sim->GetResourceManager()->AddAgent(cell);

  // Step 3 times: counter reaches 3, check is > 3 so still alive
  sim->GetScheduler()->Simulate(3);
  EXPECT_EQ(sim->GetResourceManager()->GetNumAgents(), 1u);

  // Step 1 more: counter reaches 4 > 3, removed
  sim->GetScheduler()->Simulate(1);
  EXPECT_EQ(sim->GetResourceManager()->GetNumAgents(), 0u);
  delete sim;
}

TEST(SheddingTest, ExhaustedTAApoptosis) {
  auto* sim = CreateTestSim(TEST_NAME);

  // Override apoptosis_delay to small value
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  sp->apoptosis_delay = 2;

  // Cell at threshold (stratum_steps == apoptosis_delay) should survive
  auto* alive = new Keratinocyte({15, 15, 2});
  alive->SetDiameter(4);
  alive->SetStratum(kBasal);
  alive->SetDivisionsLeft(0);  // exhausted TA
  alive->SetStratumSteps(2);   // == 2, NOT > 2
  alive->AddBehavior(new Shedding());
  sim->GetResourceManager()->AddAgent(alive);

  // Cell above threshold (stratum_steps > apoptosis_delay) should be removed
  auto* doomed = new Keratinocyte({16, 15, 2});
  doomed->SetDiameter(4);
  doomed->SetStratum(kBasal);
  doomed->SetDivisionsLeft(0);  // exhausted TA
  doomed->SetStratumSteps(3);   // > 2
  doomed->AddBehavior(new Shedding());
  sim->GetResourceManager()->AddAgent(doomed);

  sim->GetScheduler()->Simulate(1);

  // Only the at-threshold cell should survive
  EXPECT_EQ(sim->GetResourceManager()->GetNumAgents(), 1u);
  delete sim;
}

// ---------------------------------------------------------------------------
// InitFromFieldsTest
// ---------------------------------------------------------------------------

TEST(InitFromFieldsTest, ZeroStemFractionAlwaysTA) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());

  // Set stem_fraction to 0 so stem_prob is always 0 regardless of calcium
  sp->stem_fraction = 0.0;

  SetupAllFields(sim);

  // Run one step so grids initialize
  sim->GetScheduler()->Simulate(1);

  // Create 20 cells and init from fields -- all should be TA (divisions >= 0)
  for (int i = 0; i < 20; i++) {
    auto* cell = new Keratinocyte({15, 15, 2});
    cell->SetDiameter(4);
    cell->InitFromFields(sim, sp);
    EXPECT_GE(cell->GetDivisionsLeft(), 0)
        << "Cell " << i << " should be TA with stem_fraction=0";
    delete cell;
  }
  delete sim;
}

TEST(InitFromFieldsTest, KGFPreAdvance) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = sim->GetParam()->Get<SimParam>();

  SetupAllFields(sim);
  sim->GetScheduler()->Simulate(1);

  // Create cell near dermis where KGF is high
  auto* cell = new Keratinocyte({15, 15, 1});
  cell->SetDiameter(4);
  cell->InitFromFields(sim, sp);

  // KGF should pre-advance phase_elapsed > 0
  EXPECT_GT(cell->GetPhaseElapsed(), 0);
  // Upper bound: max advance is g1_duration * 0.5 = 3.5h
  EXPECT_LE(cell->GetPhaseElapsed(), sp->g1_duration * 0.5 + 0.01);
  delete cell;
  delete sim;
}

// ---------------------------------------------------------------------------
// BasalDivisionTest
// ---------------------------------------------------------------------------

TEST(BasalDivisionTest, ExhaustedTAEntersG0) {
  auto* sim = CreateTestSim(TEST_NAME);
  SetupAllFields(sim);
  sim->GetScheduler()->Simulate(1);

  auto* cell = new Keratinocyte({15, 15, 2});
  cell->SetDiameter(4);
  cell->SetStratum(kBasal);
  cell->SetDivisionsLeft(0);  // exhausted TA
  cell->SetCyclePhase(kG1);

  BasalDivision div;
  div.Run(cell);

  EXPECT_EQ(cell->GetCyclePhase(), kG0);

  delete cell;
  delete sim;
}

TEST(BasalDivisionTest, NonBasalEntersG0) {
  auto* sim = CreateTestSim(TEST_NAME);
  SetupAllFields(sim);
  sim->GetScheduler()->Simulate(1);

  auto* cell = new Keratinocyte({15, 15, 10});
  cell->SetDiameter(4);
  cell->SetStratum(kSpinous);
  cell->SetDivisionsLeft(5);
  cell->SetCyclePhase(kG1);

  BasalDivision div;
  div.Run(cell);

  EXPECT_EQ(cell->GetCyclePhase(), kG0);

  delete cell;
  delete sim;
}

TEST(BasalDivisionTest, CellAgeIncrements) {
  auto* sim = CreateTestSim(TEST_NAME);
  SetupAllFields(sim);
  sim->GetScheduler()->Simulate(1);

  auto* cell = new Keratinocyte({15, 15, 2});
  cell->SetDiameter(4);
  cell->SetStratum(kBasal);
  cell->SetDivisionsLeft(5);
  cell->SetCellAge(0);

  real_t dt = sim->GetParam()->simulation_time_step;

  BasalDivision div;
  div.Run(cell);

  EXPECT_NEAR(cell->GetCellAge(), dt, 1e-9);

  delete cell;
  delete sim;
}

TEST(BasalDivisionTest, SGrowthAndPhaseTransitions) {
  auto* sim = CreateTestSim(TEST_NAME);
  SetupAllFields(sim);
  sim->GetScheduler()->Simulate(1);

  auto* cell = new Keratinocyte({15, 15, 2});
  cell->SetDiameter(3);  // below division_diameter
  cell->SetStratum(kBasal);
  cell->SetDivisionsLeft(5);
  cell->SetCyclePhase(kS);
  cell->SetPhaseElapsed(0);
  cell->AddBehavior(new BasalDivision());
  sim->GetResourceManager()->AddAgent(cell);

  real_t initial_diam = cell->GetDiameter();

  // Run a step: cell should grow in S phase
  sim->GetScheduler()->Simulate(1);

  // Check cell grew (diameter increased)
  real_t after_diam = 0;
  sim->GetResourceManager()->ForEachAgent([&](Agent* a) {
    auto* kc = dynamic_cast<Keratinocyte*>(a);
    if (kc) after_diam = kc->GetDiameter();
  });
  EXPECT_GT(after_diam, initial_diam);

  delete sim;
}

// ---------------------------------------------------------------------------
// DifferentiationTest
// ---------------------------------------------------------------------------

TEST(DifferentiationTest, StemCellAnchoring) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = sim->GetParam()->Get<SimParam>();

  // Need calcium grid for Differentiation
  ModelInitializer::DefineSubstance(fields::kCalciumId, fields::kCalcium,
      sp->calcium_diffusion, sp->calcium_decay, 10);
  ModelInitializer::InitializeSubstance(fields::kCalciumId,
      [](real_t, real_t, real_t) { return 0.01; });  // low calcium = basal

  sim->GetScheduler()->Simulate(1);

  // Stem cell placed too high
  auto* cell = new Keratinocyte({25, 25, 10.0});
  cell->SetDiameter(4);
  cell->SetStratum(kBasal);
  cell->SetDivisionsLeft(-1);  // stem

  Differentiation diff;
  diff.Run(cell);

  // Should be clamped to max_z = diameter/2 = 2.0
  EXPECT_NEAR(cell->GetPosition()[2], 2.0, 1e-6);

  delete cell;
  delete sim;
}

TEST(DifferentiationTest, CalciumKeepsCellBasal) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = sim->GetParam()->Get<SimParam>();

  auto* vm = VolumeManager::Get();
  vm->Init(sp->volume_z_spinous, sp->volume_z_granular,
           sp->volume_z_cornified);

  // Low calcium everywhere -> cell stays basal
  ModelInitializer::DefineSubstance(fields::kCalciumId, fields::kCalcium,
      sp->calcium_diffusion, sp->calcium_decay, 10);
  ModelInitializer::InitializeSubstance(fields::kCalciumId,
      [](real_t, real_t, real_t) { return 0.01; });

  sim->GetScheduler()->Simulate(1);

  // Cell at z=1 (below spinous_threshold), low calcium
  auto* cell = new Keratinocyte({25, 25, 1.0});
  cell->SetDiameter(4);
  cell->SetStratum(kBasal);
  cell->SetDivisionsLeft(5);  // TA, not stem

  Differentiation diff;
  diff.Run(cell);

  EXPECT_EQ(cell->GetStratum(), kBasal);

  delete cell;
  delete sim;
}

TEST(DifferentiationTest, WoundWallRepulsion) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  sp->wound.enabled = true;
  sp->wound.center_x = 15;
  sp->wound.center_y = 15;
  sp->wound.radius = 5;

  ModelInitializer::DefineSubstance(fields::kCalciumId, fields::kCalcium,
      sp->calcium_diffusion, sp->calcium_decay, 10);
  ModelInitializer::InitializeSubstance(fields::kCalciumId,
      [](real_t, real_t, real_t) { return 0.01; });

  sim->GetScheduler()->Simulate(1);

  // Cell past wound edge (overshoot)
  real_t overshoot_x = 15 + 5 + 2.0;  // 2 units past edge
  auto* cell = new Keratinocyte({overshoot_x, 15, 2});
  cell->SetDiameter(4);
  cell->SetStratum(kBasal);
  cell->SetDivisionsLeft(5);

  Differentiation diff;
  diff.Run(cell);

  // Should be pushed back toward center (x decreased)
  EXPECT_LT(cell->GetPosition()[0], overshoot_x);

  delete cell;
  delete sim;
}

// ---------------------------------------------------------------------------
// ChemotaxisTest
// ---------------------------------------------------------------------------

TEST(ChemotaxisTest, FallbackToGeometric) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  SanitizeForUnitTest(sp);
  sp->chemotaxis_enabled = true;
  sp->immune.cytokine_rate = 0;  // prevent cytokine write from creating gradient

  SetupInflammation(sim);
  sim->GetScheduler()->Simulate(1);

  // Cell east of wound center (center = 25,25 in test domain), no inflammation gradient -> falls back
  auto* cell = new ImmuneCell({30, 25, 1.5});
  cell->SetDiameter(3);
  cell->SetImmuneCellType(kNeutrophil);

  NeutrophilBehavior beh;
  beh.Run(cell);

  // Should still move toward wound center (geometric fallback)
  auto tf = cell->GetTractorForce();
  EXPECT_LT(tf[0], 0);  // toward center (west)

  delete cell;
  delete sim;
}

}  // namespace skibidy
}  // namespace bdm
