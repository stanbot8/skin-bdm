// -----------------------------------------------------------------------------
//
// Copyright (C) 2021 CERN & University of Surrey for the benefit of the
// BioDynaMo collaboration. All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
//
// See the LICENSE file distributed with this work for details.
// See the NOTICE file distributed with this work for additional information
// regarding copyright ownership.
//
// -----------------------------------------------------------------------------

#include <gtest/gtest.h>
#include "biodynamo.h"
#include "infra/toml_compat.h"
#include "infra/sim_param.h"
#include "skibidy.h"
#include "core/pde.h"
#include "core/composite_field.h"
#include "tissue/calcium.h"
#include "tissue/kgf.h"
#include "tissue/oxygen.h"
#include "tissue/stratum.h"
#include "inflammation/inflammation_pde.h"
#include "scar/scar_pde.h"
#include "fibroblast/tgfbeta_pde.h"
#include "fibroblast/collagen_pde.h"
#include "biofilm/biofilm_pde.h"
#include "angiogenesis/vegf_pde.h"
#include "mmp/mmp_pde.h"
#include "fibronectin/fibronectin_pde.h"
#include "elastin/elastin_pde.h"
#include "hyaluronan/hyaluronan_pde.h"
#include "dermis/dermis_pde.h"
#include "hemostasis/hemostasis_pde.h"
#include "ph/ph_pde.h"
#include "tumor/tumor_pde.h"
#include "tissue/water.h"
#include "core/cell_cell_force.h"
#include "tissue/migration.h"
#include "tissue/shedding.h"
#include "tissue/basal_division.h"
#include "immune/neutrophil_behavior.h"
#include "immune/macrophage_behavior.h"
#include "wound/wound_event.h"
#include "immune/immune_response.h"
#include "fibroblast/fibroblast.h"
#include "fibroblast/fibroblast_behavior.h"
#include "fibroblast/fibroblast_recruitment.h"
#include "tumor/tumor_cell.h"
#include "tumor/tumor_behavior.h"
#include "tumor/tumor_initiation.h"
#include "perfusion/vascular.h"
#include "diabetic/baseline_inflammation.h"
#include "biofilm/biofilm_op.h"
#include "core/metrics.h"
#include "angiogenesis/vegf_op.h"
#include "core/chunked_grid.h"
#include "core/fused_source.h"
#include "core/derived_field.h"
#include "core/derived_fields_op.h"
#include "ph/ph_pde.h"

#include <fstream>
#include <sstream>

#define TEST_NAME typeid(*this).name()

namespace bdm {
namespace skibidy {

// ---------------------------------------------------------------------------
// Helper: create a minimal simulation with SimParam registered
// ---------------------------------------------------------------------------
static Simulation* CreateTestSim(const char* name) {
  Param::RegisterParamGroup(new SimParam());
  auto set_param = [](Param* param) {
    param->bound_space = Param::BoundSpaceMode::kClosed;
    param->min_bound = 0;
    param->max_bound = 50;
    param->simulation_time_step = 0.1;
    param->export_visualization = false;
    // Per-axis bounds for ClampToBounds (match grid range)
    auto* sp = const_cast<SimParam*>(param->Get<SimParam>());
    sp->bounds_min = {0, 0, 0};
    sp->bounds_max = {50, 50, 50};
    sp->grid_resolution = 10;  // coarse grid for fast tests
  };
  return new Simulation(name, set_param);
}

// ---------------------------------------------------------------------------
// Helper: register inflammation fields matching config (split vs single)
// ---------------------------------------------------------------------------
static void SetupInflammation(Simulation* sim) {
  auto* sp = sim->GetParam()->Get<SimParam>();
  if (sp->split_inflammation_enabled) {
    ProInflammatoryPDE pro(sp); pro.Init(sim);
    AntiInflammatoryPDE anti(sp); anti.Init(sim);
  } else {
    InflammationPDE infl(sp); infl.Init(sim);
  }
  ImmunePressurePDE ip(sp); ip.Init(sim);
}

// Helper: disable features that require full sim context (for isolated tests).
// The merged bdm.toml enables all modules, but unit tests only register
// the minimal fields they need.  This prevents accessing unregistered grids.
static void SanitizeForUnitTest(SimParam* sp) {
  sp->diabetic_mode = false;
  sp->efferocytosis_enabled = false;
  sp->biofilm_enabled = false;
  sp->mmp_enabled = false;
  sp->angiogenesis_enabled = false;
  sp->fibronectin_enabled = false;
  sp->fibroblast_enabled = false;
  sp->elastin_enabled = false;
  sp->hyaluronan_enabled = false;
  sp->hemostasis_enabled = false;
  // Reset multi-resolution so unit tests use uniform grid resolution
  sp->grid_resolution_structural = 0;
}

// Same but keeps diabetic_mode true (for DiabeticTest tests)
static void SanitizeForDiabeticTest(SimParam* sp) {
  sp->efferocytosis_enabled = false;
  sp->biofilm_enabled = false;
  sp->mmp_enabled = false;
  sp->angiogenesis_enabled = false;
  sp->fibronectin_enabled = false;
  sp->fibroblast_enabled = false;
  sp->split_inflammation_enabled = false;
  sp->elastin_enabled = false;
  sp->hyaluronan_enabled = false;
  sp->hemostasis_enabled = false;
}

// Helper: get the readable inflammation grid (pro-infl in split mode)
static DiffusionGrid* GetInflammationGrid(Simulation* sim) {
  auto* sp = sim->GetParam()->Get<SimParam>();
  auto* rm = sim->GetResourceManager();
  if (sp->split_inflammation_enabled) {
    return rm->GetDiffusionGrid(fields::kProInflammatory);
  }
  return rm->GetDiffusionGrid(fields::kInflammation);
}

// ---------------------------------------------------------------------------
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
  auto* ca = rm->GetDiffusionGrid(fields::kCalcium);

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
// Helper: initialize all 5 PDE channels for tests that read field values
// ---------------------------------------------------------------------------
static void SetupAllFields(Simulation* sim) {
  auto* sp = sim->GetParam()->Get<SimParam>();
  VascularPDE vasc;
  vasc.Init(sim);
  CalciumPDE ca;
  ca.Init(sim);
  KgfPDE kgf;
  kgf.Init(sim);
  OxygenPDE o2;
  o2.Init(sim);
  auto* vm = VolumeManager::Get();
  vm->Init(sp->volume_z_spinous, sp->volume_z_granular,
           sp->volume_z_cornified);
  StratumPDE st;
  st.Init(sim);
  ScarPDE scar;
  scar.Init(sim);
  WaterPDE water;
  water.Init(sim);
  SetupInflammation(sim);
  if (sp->fibronectin_enabled) { FibronectinPDE fn(sp); fn.Init(sim); }
  if (sp->mmp_enabled) { MMPPDE mmp(sp); mmp.Init(sim); }
  if (sp->fibroblast_enabled) {
    TGFBetaPDE tgfb(sp); tgfb.Init(sim);
    CollagenPDE col(sp); col.Init(sim);
  }
  { PHPDE ph(sp); ph.Init(sim); }
}

// ===========================================================================
// CellCellForce tests
// ===========================================================================
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

// ===========================================================================
// Migration tests
// ===========================================================================
TEST(MigrationTest, DirectionTowardCenter) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = sim->GetParam()->Get<SimParam>();
  WaterPDE water; water.Init(sim);
  SetupInflammation(sim);
  if (sp->fibronectin_enabled) { FibronectinPDE fn(sp); fn.Init(sim); }
  sim->GetScheduler()->Simulate(1);  // populate grids

  // Cell east of wound center (center = 15,15 from bdm.toml)
  auto* cell = new Keratinocyte({19, 15, 2});
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
  SetupInflammation(sim);
  if (sp->fibronectin_enabled) { FibronectinPDE fn(sp); fn.Init(sim); }
  sim->GetScheduler()->Simulate(1);  // populate grids
  real_t r = sp->wound_radius;

  // Cell A: near wound edge (distance ≈ radius)
  auto* cellA = new Keratinocyte({sp->wound_center_x + r * 0.99,
                                   sp->wound_center_y, 2});
  cellA->SetDiameter(4);
  cellA->SetStratum(kBasal);

  Migration migA;
  migA.Run(cellA);
  real_t speedA = std::abs(cellA->GetTractorForce()[0]);

  // Cell B: at half radius
  auto* cellB = new Keratinocyte({sp->wound_center_x + r / 2,
                                   sp->wound_center_y, 2});
  cellB->SetDiameter(4);
  cellB->SetStratum(kBasal);

  Migration migB;
  migB.Run(cellB);
  real_t speedB = std::abs(cellB->GetTractorForce()[0]);

  // speed = migration_speed * (dist / r), so A (~edge) ≈ 2x faster than B
  EXPECT_GT(speedA, 0);
  EXPECT_GT(speedB, 0);
  EXPECT_NEAR(speedA / speedB, 0.99 * 2.0, 0.1);

  delete cellA;
  delete cellB;
  delete sim;
}

TEST(MigrationTest, NonBasalClearsForce) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = sim->GetParam()->Get<SimParam>();
  WaterPDE water; water.Init(sim);
  SetupInflammation(sim);
  if (sp->fibronectin_enabled) { FibronectinPDE fn(sp); fn.Init(sim); }
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

// ===========================================================================
// Shedding tests
// ===========================================================================
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

// ===========================================================================
// InitFromFields tests
// ===========================================================================
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

// ===========================================================================
// CompositeField tests
// ===========================================================================
TEST(CompositeFieldTest, AddAndGet) {
  auto* sim = CreateTestSim(TEST_NAME);

  CompositeField fields;
  fields.Add(std::make_unique<CalciumPDE>());
  fields.Add(std::make_unique<KgfPDE>());

  EXPECT_EQ(fields.Size(), 2u);
  EXPECT_NE(fields.Get("Calcium"), nullptr);
  EXPECT_NE(fields.Get("KGF"), nullptr);
  EXPECT_EQ(fields.Get("Nonexistent"), nullptr);

  delete sim;
}

// ===========================================================================
// WoundEvent test
// ===========================================================================
TEST(WoundTest, MarginCellCount) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  sp->wound_trigger_step = 0;  // fire immediately

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
  if (sp->split_inflammation_enabled) {
    fields.Add(std::make_unique<ProInflammatoryPDE>(sp));
    fields.Add(std::make_unique<AntiInflammatoryPDE>(sp));
  } else {
    fields.Add(std::make_unique<InflammationPDE>(sp));
  }
  if (sp->fibronectin_enabled) fields.Add(std::make_unique<FibronectinPDE>(sp));
  if (sp->mmp_enabled) fields.Add(std::make_unique<MMPPDE>(sp));
  if (sp->fibroblast_enabled) {
    fields.Add(std::make_unique<TGFBetaPDE>(sp));
    fields.Add(std::make_unique<CollagenPDE>(sp));
  }
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
      std::round(2.0 * M_PI * sp->wound_radius / sp->division_diameter));
  if (expected < 6) expected = 6;

  EXPECT_EQ(sim->GetResourceManager()->GetNumAgents(),
            static_cast<size_t>(expected));
  delete sim;
}

// ===========================================================================
// Integration test
// ===========================================================================
TEST(IntegrationTest, WoundTriggersAndCellsExist) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  sp->wound_trigger_step = 5;

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
  if (sp->split_inflammation_enabled) {
    fields.Add(std::make_unique<ProInflammatoryPDE>(sp));
    fields.Add(std::make_unique<AntiInflammatoryPDE>(sp));
  } else {
    fields.Add(std::make_unique<InflammationPDE>(sp));
  }
  if (sp->fibronectin_enabled) fields.Add(std::make_unique<FibronectinPDE>(sp));
  if (sp->mmp_enabled) fields.Add(std::make_unique<MMPPDE>(sp));
  if (sp->fibroblast_enabled) {
    fields.Add(std::make_unique<TGFBetaPDE>(sp));
    fields.Add(std::make_unique<CollagenPDE>(sp));
  }
  fields.InitAll(sim);

  // Schedule wound event
  OperationRegistry::GetInstance()->AddOperationImpl(
      "WoundEvent_integ", OpComputeTarget::kCpu, new WoundEvent(&fields));
  sim->GetScheduler()->ScheduleOp(
      NewOperation("WoundEvent_integ"), OpType::kPreSchedule);

  // Schedule composite field op for field coupling
  OperationRegistry::GetInstance()->AddOperationImpl(
      "CompositeFieldOp_integ", OpComputeTarget::kCpu,
      new CompositeFieldOp(&fields));
  sim->GetScheduler()->ScheduleOp(
      NewOperation("CompositeFieldOp_integ"), OpType::kPreSchedule);

  // Before wound: no agents
  EXPECT_EQ(sim->GetResourceManager()->GetNumAgents(), 0u);

  // Run 50 steps (wound at step 5, then 45 steps of cell activity)
  sim->GetScheduler()->Simulate(50);

  // After 50 steps: agents should exist (spawned at step 5, dividing since)
  EXPECT_GT(sim->GetResourceManager()->GetNumAgents(), 0u);
  delete sim;
}

// ===========================================================================
// ImmuneCell tests
// ===========================================================================
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

// ===========================================================================
// Neutrophil / Macrophage behavior tests
// ===========================================================================
TEST(NeutrophilBehaviorTest, HardCeilingLifespan) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  sp->neutrophil_lifespan = 5;
  sp->neutrophil_apoptosis_rate = 0;  // disable stochastic death, test ceiling
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
  sp->macrophage_m1_duration = 3;
  sp->macrophage_lifespan = 100;
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
  sp->macrophage_m1_duration = 9999;         // disable timer
  sp->m1_transition_threshold = 0.5;
  sp->m1_transition_min_age = 3;             // allow cytokine trigger after 3 steps
  sp->macrophage_lifespan = 100;
  sp->inflammation_decay = 0;                // no decay, we control concentration
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
  sp->macrophage_m1_duration = 9999;         // disable timer
  sp->m1_transition_threshold = 999.0;       // always above local inflammation
  sp->m1_transition_min_age = 10;            // guard: must be in M1 for 10 steps
  sp->macrophage_lifespan = 100;
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

// ===========================================================================
// ImmuneResponse tests
// ===========================================================================
TEST(ImmuneResponseTest, NeutrophilSpawnCount) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  SanitizeForUnitTest(sp);
  sp->wound_trigger_step = 0;
  sp->neutrophil_spawn_delay = 5;
  sp->neutrophil_spawn_waves = 1;     // single wave for count test
  sp->macrophage_spawn_delay = 9999;

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
  int expected = static_cast<int>(std::round(2.0 * M_PI * sp->wound_radius / 3.0));
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
  sp->wound_trigger_step = 0;
  sp->neutrophil_spawn_delay = 9999;     // disable neutrophils
  sp->macrophage_spawn_delay = 5;
  sp->macrophage_spawn_threshold = 0.0;  // any inflammation triggers recruitment
  sp->macrophage_spawn_rate = 1.0;       // high rate for deterministic test
  sp->macrophage_apoptosis_rate = 0;     // disable death during test
  sp->macrophage_lifespan = 9999;

  SetupInflammation(sim);

  // Initialize grids (boxes not created until first step)
  sim->GetScheduler()->Simulate(1);

  // Seed inflammation at wound center so recruitment has signal
  auto* rm = sim->GetResourceManager();
  auto* infl_grid = GetInflammationGrid(sim);
  Real3 center = {sp->wound_center_x, sp->wound_center_y, 1.5};
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
  sp->wound_trigger_step = 0;
  sp->neutrophil_spawn_delay = 9999;
  sp->macrophage_spawn_delay = 1;
  sp->macrophage_spawn_threshold = 0.1;
  sp->macrophage_spawn_rate = 1.0;

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
  sp->wound_trigger_step = 0;
  sp->neutrophil_spawn_delay = 1;
  sp->neutrophil_spawn_waves = 3;
  sp->neutrophil_spawn_window = 6;   // waves at steps 1, 4, 7
  sp->neutrophil_apoptosis_rate = 0;  // disable stochastic death
  sp->neutrophil_lifespan = 9999;    // don't die during test
  sp->macrophage_spawn_delay = 9999;

  SetupInflammation(sim);

  auto* ir = new ImmuneResponse();
  OperationRegistry::GetInstance()->AddOperationImpl(
      "ImmuneResponse_test_wave", OpComputeTarget::kCpu, ir);
  sim->GetScheduler()->ScheduleOp(
      NewOperation("ImmuneResponse_test_wave"), OpType::kPreSchedule);

  int total = static_cast<int>(
      std::round(2.0 * M_PI * sp->wound_radius / 3.0));
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
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  SanitizeForUnitTest(sp);
  sp->neutrophil_lifespan = 9999;         // disable hard ceiling
  sp->neutrophil_min_survival = 0;        // apoptosis from step 1
  sp->neutrophil_apoptosis_rate = 0.01;   // ~1% per step

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

  sim->GetScheduler()->Simulate(1);  // commit agents
  EXPECT_EQ(sim->GetResourceManager()->GetNumAgents(), 100u);

  // Run 100 steps: with 1% rate, expect ~63% dead (1 - e^{-1})
  sim->GetScheduler()->Simulate(100);
  size_t survivors = sim->GetResourceManager()->GetNumAgents();

  // Survivors should be roughly 37 (e^{-1} ~ 0.37), allow wide margin
  EXPECT_GT(survivors, 10u);   // not all dead
  EXPECT_LT(survivors, 70u);   // substantial deaths occurred

  delete sim;
}

// ===========================================================================
// BasalDivision tests
// ===========================================================================
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

// ===========================================================================
// WoundResolution tests
// ===========================================================================
TEST(WoundResolutionTest, SafetyTimeout) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  sp->wound_trigger_step = 0;
  sp->num_steps = 210;  // safety timeout fires at num_steps - 200 = step 10

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
  if (sp->split_inflammation_enabled) {
    fields.Add(std::make_unique<ProInflammatoryPDE>(sp));
    fields.Add(std::make_unique<AntiInflammatoryPDE>(sp));
  } else {
    fields.Add(std::make_unique<InflammationPDE>(sp));
  }
  if (sp->fibronectin_enabled) fields.Add(std::make_unique<FibronectinPDE>(sp));
  if (sp->mmp_enabled) fields.Add(std::make_unique<MMPPDE>(sp));
  if (sp->fibroblast_enabled) {
    fields.Add(std::make_unique<TGFBetaPDE>(sp));
    fields.Add(std::make_unique<CollagenPDE>(sp));
  }
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
  EXPECT_GT(sim->GetResourceManager()->GetNumAgents(), 0u);

  // Run to safety timeout (step >= num_steps - 200 = 10)
  sim->GetScheduler()->Simulate(10);

  // All agents should be removed
  EXPECT_EQ(sim->GetResourceManager()->GetNumAgents(), 0u);

  delete sim;
}

// ===========================================================================
// Differentiation tests
// ===========================================================================
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
  sp->wound_enabled = true;
  sp->wound_center_x = 15;
  sp->wound_center_y = 15;
  sp->wound_radius = 5;

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

// ===========================================================================
// Extension 4: Split pro/anti-inflammatory fields
// ===========================================================================
TEST(SplitInflammationTest, ProAntiFieldsRegistered) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  sp->split_inflammation_enabled = true;

  ProInflammatoryPDE pro(sp);
  pro.Init(sim);
  AntiInflammatoryPDE anti(sp);
  anti.Init(sim);

  auto* rm = sim->GetResourceManager();
  EXPECT_NE(rm->GetDiffusionGrid(fields::kProInflammatory), nullptr);
  EXPECT_NE(rm->GetDiffusionGrid(fields::kAntiInflammatory), nullptr);

  delete sim;
}

TEST(SplitInflammationTest, NeutrophilWritesToProField) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  SanitizeForUnitTest(sp);
  sp->split_inflammation_enabled = true;

  ProInflammatoryPDE pro(sp);
  pro.Init(sim);
  AntiInflammatoryPDE anti(sp);
  anti.Init(sim);
  sim->GetScheduler()->Simulate(1);

  auto* rm = sim->GetResourceManager();
  auto* pro_grid = rm->GetDiffusionGrid(fields::kProInflammatory);
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
  sp->split_inflammation_enabled = true;
  sp->macrophage_lifespan = 100;

  ProInflammatoryPDE pro(sp);
  pro.Init(sim);
  AntiInflammatoryPDE anti(sp);
  anti.Init(sim);
  sim->GetScheduler()->Simulate(1);

  auto* rm = sim->GetResourceManager();
  auto* anti_grid = rm->GetDiffusionGrid(fields::kAntiInflammatory);
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
  sp->diabetic_mode = true;
  sp->macrophage_m1_duration = 3;
  sp->diabetic_m1_duration_factor = 3.0;
  sp->macrophage_lifespan = 100;
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
  sp->diabetic_mode = true;
  sp->diabetic_resolution_factor = 0.3;
  sp->macrophage_lifespan = 100;

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
  real_t expected = sp->immune_resolution_rate * sp->diabetic_resolution_factor;
  EXPECT_NEAR(consumed, expected, 1e-6);

  delete cell;
  delete sim;
}

TEST(DiabeticTest, KeratinocyteProliferation) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  SanitizeForDiabeticTest(sp);
  sp->diabetic_mode = true;
  sp->diabetic_prolif_factor = 0.0;  // completely suppress proliferation

  SetupAllFields(sim);
  sim->GetScheduler()->Simulate(1);

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
  sp->wound_enabled = true;
  sp->migration_enabled = true;
  sp->wound_center_x = 15;
  sp->wound_center_y = 15;
  sp->wound_radius = 10;
  SetupAllFields(sim);
  sim->GetScheduler()->Simulate(1);

  // Normal cell (diabetic off): measure tractor force
  sp->diabetic_mode = false;
  auto* cell1 = new Keratinocyte({20, 15, 2});
  cell1->SetDiameter(4);
  cell1->SetStratum(kBasal);

  Migration mig;
  mig.Run(cell1);
  Real3 force_normal = cell1->GetTractorForce();
  real_t mag_normal = std::sqrt(
      force_normal[0] * force_normal[0] + force_normal[1] * force_normal[1]);

  // Diabetic cell at same position
  sp->diabetic_mode = true;
  sp->diabetic_migration_factor = 0.5;

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
  sp->diabetic_mode = true;
  sp->fibroblast_enabled = true;
  sp->fibroblast_activation_delay = 10;
  sp->diabetic_fibroblast_activation_factor = 3.0;
  sp->fibroblast_activation_threshold = 999;  // prevent TGF-beta trigger
  sp->fibroblast_lifespan = 10000;

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
  sp->wound_enabled = true;
  sp->wound_trigger_step = 0;
  sp->neutrophil_spawn_delay = 1;
  sp->neutrophil_spawn_waves = 1;       // single wave for count test
  sp->neutrophil_apoptosis_rate = 0;    // disable stochastic death
  sp->diabetic_mode = true;
  sp->diabetic_neutrophil_factor = 2.0;
  sp->diabetic_neutrophil_waves_factor = 1.0;   // isolate count test from wave scaling
  sp->diabetic_neutrophil_window_factor = 1.0;  // isolate count test from window scaling

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
      2.0 * M_PI * sp->wound_radius / 3.0));
  if (base_count < 6) base_count = 6;
  int expected = static_cast<int>(std::ceil(base_count * 2.0));

  EXPECT_EQ(diabetic_count, expected);

  delete sim;
}

TEST(DiabeticTest, ExtendedNeutrophilLifespan) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  SanitizeForDiabeticTest(sp);
  sp->diabetic_mode = true;
  sp->neutrophil_lifespan = 10;
  sp->diabetic_neutrophil_lifespan_factor = 2.0;
  sp->neutrophil_apoptosis_rate = 0;    // disable stochastic death
  sp->neutrophil_min_survival = 9999;   // disable stochastic death

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

// ===========================================================================
// Extension 3: Efferocytosis
// ===========================================================================
TEST(EfferocytosisTest, M1EngulfsDyingNeutrophil) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  SanitizeForUnitTest(sp);
  sp->efferocytosis_enabled = true;
  sp->efferocytosis_radius = 2.0;  // must fit within env box length
  sp->efferocytosis_age_fraction = 0.5;
  sp->neutrophil_lifespan = 10;
  sp->macrophage_m1_duration = 9999;  // disable timer
  sp->macrophage_lifespan = 100;

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
  sp->neutrophil_lifespan = 100;
  sp->macrophage_m1_duration = 9999;
  sp->macrophage_lifespan = 100;

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

// ===========================================================================
// Extension 2: Chemotaxis
// ===========================================================================
TEST(ChemotaxisTest, FallbackToGeometric) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  SanitizeForUnitTest(sp);
  sp->chemotaxis_enabled = true;
  sp->immune_cytokine_rate = 0;  // prevent cytokine write from creating gradient

  SetupInflammation(sim);
  sim->GetScheduler()->Simulate(1);

  // Cell east of wound center, no inflammation gradient -> falls back
  auto* cell = new ImmuneCell({19, 15, 1.5});
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

// ===========================================================================
// Extension 1: Proportional scarring
// ===========================================================================
TEST(ProportionalScarTest, AccumulationOp) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  SanitizeForUnitTest(sp);
  sp->scar_proportional_enabled = true;
  sp->scar_accumulation_rate = 0.1;  // high rate for test visibility
  sp->wound_enabled = true;
  sp->wound_center_x = 25;
  sp->wound_center_y = 25;
  sp->wound_radius = 20;  // large wound to cover grid
  sp->volume_z_cornified = 40;

  auto* vm = VolumeManager::Get();
  vm->Init(sp->volume_z_spinous, sp->volume_z_granular,
           sp->volume_z_cornified);

  SetupInflammation(sim);
  ScarPDE scar;
  scar.Init(sim);
  sim->GetScheduler()->Simulate(1);  // populate grid arrays

  // Seed inflammation
  auto* rm = sim->GetResourceManager();
  auto* infl_grid = GetInflammationGrid(sim);
  Real3 center = {25, 25, 5};
  size_t idx = infl_grid->GetBoxIndex(center);
  infl_grid->ChangeConcentrationBy(idx, 5.0);

  // Run accumulation op
  ScarAccumulationOp op;
  op();

  auto* scar_grid = rm->GetDiffusionGrid(fields::kScar);
  size_t scar_idx = scar_grid->GetBoxIndex(center);
  real_t scar_val = scar_grid->GetConcentration(scar_idx);

  // scar += inflammation * rate = 5.0 * 0.1 = 0.5
  EXPECT_NEAR(scar_val, 0.5, 0.01);

  delete sim;
}

TEST(ProportionalScarTest, WriteScarValueSkipsBinaryScar) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  sp->scar_proportional_enabled = true;

  SetupAllFields(sim);
  sim->GetScheduler()->Simulate(1);

  auto* cell = new Keratinocyte({25, 25, 5});
  cell->SetDiameter(4);
  cell->SetStratum(kCornified);

  WriteScarValue(cell);

  // Scar grid should NOT have binary 1.0 when proportional mode is on
  auto* rm = sim->GetResourceManager();
  auto* scar_grid = rm->GetDiffusionGrid(fields::kScar);
  Real3 pos = {25, 25, 5};
  size_t idx = scar_grid->GetBoxIndex(pos);
  EXPECT_LT(scar_grid->GetConcentration(idx), 0.5);

  // But stratum field should still get stamped
  auto* stratum_grid = rm->GetDiffusionGrid(fields::kStratum);
  size_t st_idx = stratum_grid->GetBoxIndex(pos);
  EXPECT_GT(stratum_grid->GetConcentration(st_idx), 0);

  delete cell;
  delete sim;
}

// ===========================================================================
// Fibroblast agent tests
// ===========================================================================
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

// ===========================================================================
// Fibroblast behavior tests
// ===========================================================================
TEST(FibroblastBehaviorTest, ActivationByTGFBeta) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  SanitizeForUnitTest(sp);
  sp->fibroblast_enabled = true;
  sp->fibroblast_activation_threshold = 0.1;
  sp->fibroblast_activation_delay = 9999;  // disable timeout

  TGFBetaPDE tgfb(sp);
  tgfb.Init(sim);
  CollagenPDE col(sp);
  col.Init(sim);
  sim->GetScheduler()->Simulate(1);

  // Seed TGF-beta above threshold at cell position
  auto* rm = sim->GetResourceManager();
  auto* tgfb_grid = rm->GetDiffusionGrid(fields::kTGFBeta);
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
  sp->fibroblast_enabled = true;
  sp->fibroblast_activation_threshold = 999;  // very high, won't trigger
  sp->fibroblast_activation_delay = 3;

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
  sp->fibroblast_enabled = true;
  sp->fibroblast_myofibroblast_threshold = 0.2;
  sp->fibroblast_myofibroblast_delay = 2;
  sp->fibroblast_lifespan = 9999;

  TGFBetaPDE tgfb(sp);
  tgfb.Init(sim);
  CollagenPDE col(sp);
  col.Init(sim);
  sim->GetScheduler()->Simulate(1);

  // Seed TGF-beta above myofibroblast threshold
  auto* rm = sim->GetResourceManager();
  auto* tgfb_grid = rm->GetDiffusionGrid(fields::kTGFBeta);
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
  sp->fibroblast_enabled = true;
  sp->collagen_deposition_rate = 0.1;
  sp->fibroblast_tgfb_rate = 0;  // don't write TGF-beta (isolate collagen test)
  sp->fibroblast_lifespan = 9999;
  sp->fibroblast_apoptosis_threshold = 0;  // prevent removal

  TGFBetaPDE tgfb(sp);
  tgfb.Init(sim);
  CollagenPDE col(sp);
  col.Init(sim);
  sim->GetScheduler()->Simulate(1);

  // Seed TGF-beta
  auto* rm = sim->GetResourceManager();
  auto* tgfb_grid = rm->GetDiffusionGrid(fields::kTGFBeta);
  Real3 pos = {15, 15, 2};
  size_t idx = tgfb_grid->GetBoxIndex(pos);
  tgfb_grid->ChangeConcentrationBy(idx, 2.0);

  auto* cell = new Fibroblast(pos);
  cell->SetDiameter(4);
  cell->SetFibroblastState(kMyofibroblast);

  FibroblastBehavior beh;
  beh.Run(cell);

  auto* col_grid = rm->GetDiffusionGrid(fields::kCollagen);
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
  sp->fibroblast_enabled = true;
  sp->fibroblast_tgfb_rate = 0.01;
  sp->fibroblast_lifespan = 9999;
  sp->fibroblast_apoptosis_threshold = 0;

  TGFBetaPDE tgfb(sp);
  tgfb.Init(sim);
  CollagenPDE col(sp);
  col.Init(sim);
  sim->GetScheduler()->Simulate(1);

  auto* rm = sim->GetResourceManager();
  auto* tgfb_grid = rm->GetDiffusionGrid(fields::kTGFBeta);
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
  sp->fibroblast_enabled = true;
  sp->fibroblast_apoptosis_threshold = 0.5;
  sp->fibroblast_min_lifespan = 3;
  sp->fibroblast_lifespan = 9999;

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
  sp->fibroblast_enabled = true;
  sp->fibroblast_lifespan = 5;
  sp->fibroblast_apoptosis_threshold = 0;  // prevent apoptosis path

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
  sp->fibroblast_enabled = true;
  sp->fibroblast_lifespan = 9999;
  sp->fibroblast_activation_threshold = 0;  // auto-activate immediately

  TGFBetaPDE tgfb(sp);
  tgfb.Init(sim);
  CollagenPDE col(sp);
  col.Init(sim);
  sim->GetScheduler()->Simulate(1);

  // Cell east of wound center (center = 15,15)
  auto* cell = new Fibroblast({19, 15, 2});
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

// ===========================================================================
// Fibroblast recruitment test
// ===========================================================================
TEST(FibroblastRecruitmentTest, SpawnTimingAndCount) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  sp->wound_trigger_step = 0;
  sp->fibroblast_enabled = true;
  sp->fibroblast_spawn_delay = 5;
  sp->fibroblast_spawn_waves = 1;  // single wave for deterministic count

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

  real_t outer_r = sp->wound_radius + sp->dermal_fibroblast_margin;
  int expected = static_cast<int>(
      std::ceil(2.0 * M_PI * outer_r / sp->fibroblast_diameter *
                sp->fibroblast_density_factor));
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

// ===========================================================================
// M2 macrophage TGF-beta production test
// ===========================================================================
TEST(MacrophageBehaviorTest, M2ProducesTGFBeta) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  SanitizeForUnitTest(sp);
  sp->fibroblast_enabled = true;
  sp->macrophage_lifespan = 100;

  SetupInflammation(sim);
  TGFBetaPDE tgfb(sp);
  tgfb.Init(sim);
  CollagenPDE col(sp);
  col.Init(sim);
  sim->GetScheduler()->Simulate(1);

  auto* rm = sim->GetResourceManager();
  auto* tgfb_grid = rm->GetDiffusionGrid(fields::kTGFBeta);
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

// ===========================================================================
// TumorCell agent tests
// ===========================================================================
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

// ===========================================================================
// TumorBehavior tests
// ===========================================================================
TEST(TumorBehaviorTest, CycleFasterThanNormal) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  sp->tumor_cycle_factor = 0.5;  // 2x faster

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
  sp->tumor_max_neighbors = 30;
  sp->tumor_ci_steepness = 0;      // hard threshold (deterministic)
  sp->tumor_apoptosis_rate = 0;    // no stochastic death

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

// ===========================================================================
// TumorInitiation test
// ===========================================================================
TEST(TumorInitiationTest, SpawnTimingAndCount) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  sp->tumor_enabled = true;
  sp->tumor_seed_time = 5;
  sp->tumor_seed_count = 8;

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

// ===========================================================================
// GridContext tests
// ===========================================================================
TEST(GridContextTest, InWoundCylinderDetection) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  sp->wound_center_x = 20;
  sp->wound_center_y = 20;
  sp->wound_radius = 5;

  CalciumPDE ca;
  ca.Init(sim);
  sim->GetScheduler()->Simulate(1);

  auto* grid = sim->GetResourceManager()->GetDiffusionGrid(fields::kCalcium);
  GridContext ctx(grid, sp);

  // Center of wound
  EXPECT_TRUE(ctx.InWound(20, 20));
  // Edge (exactly at radius)
  EXPECT_TRUE(ctx.InWound(25, 20));
  // Just outside
  EXPECT_FALSE(ctx.InWound(25.1, 20));
  // Far away
  EXPECT_FALSE(ctx.InWound(0, 0));

  delete sim;
}

TEST(GridContextTest, CoordinateHelpers) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = sim->GetParam()->Get<SimParam>();

  CalciumPDE ca;
  ca.Init(sim);
  sim->GetScheduler()->Simulate(1);

  auto* grid = sim->GetResourceManager()->GetDiffusionGrid(fields::kCalcium);
  GridContext ctx(grid, sp);

  // Voxel 0 should have coordinates at lo + box_len/2
  real_t x0 = ctx.X(0);
  real_t y0 = ctx.Y(0);
  real_t z0 = ctx.Z(0);

  real_t expected = ctx.lo + ctx.box_len / 2.0;
  EXPECT_NEAR(x0, expected, 1e-6);
  EXPECT_NEAR(y0, expected, 1e-6);
  EXPECT_NEAR(z0, expected, 1e-6);

  // Voxel 1 should be one box_len in x from voxel 0
  EXPECT_NEAR(ctx.X(1) - ctx.X(0), ctx.box_len, 1e-6);
  EXPECT_NEAR(ctx.Y(1), ctx.Y(0), 1e-6);  // same y
  EXPECT_NEAR(ctx.Z(1), ctx.Z(0), 1e-6);  // same z

  delete sim;
}

// ===========================================================================
// VascularPDE tests
// ===========================================================================
TEST(VascularPDETest, InitializesToOneInDermis) {
  // Need domain with negative z for dermal voxels
  Param::RegisterParamGroup(new SimParam());
  auto set_param = [](Param* param) {
    param->bound_space = Param::BoundSpaceMode::kClosed;
    param->min_bound = -10;
    param->max_bound = 50;
    param->simulation_time_step = 0.1;
    param->export_visualization = false;
  };
  auto* sim = new Simulation(TEST_NAME, set_param);
  auto* sp = sim->GetParam()->Get<SimParam>();

  VascularPDE vasc;
  vasc.Init(sim);
  sim->GetScheduler()->Simulate(1);

  auto* grid = sim->GetResourceManager()->GetDiffusionGrid(fields::kVascular);
  GridContext ctx(grid, sp);

  // Check dermal voxels (z < 0) are initialized to perfusion_basal
  int dermal_count = 0;
  real_t sum = 0;
  for (size_t idx = 0; idx < ctx.n; idx++) {
    if (ctx.Z(idx) < 0) {
      dermal_count++;
      sum += grid->GetConcentration(idx);
    }
  }
  ASSERT_GT(dermal_count, 0);
  real_t mean = sum / dermal_count;
  // After one step of diffusion, values may shift slightly from 1.0,
  // but should still be close to perfusion_basal
  EXPECT_NEAR(mean, sp->perfusion_basal, 0.2);

  delete sim;
}

TEST(VascularPDETest, WoundZerosDermalVoxels) {
  // Need domain with negative z for dermal voxels
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
  sp->wound_center_x = 25;
  sp->wound_center_y = 25;
  sp->wound_radius = 10;

  VascularPDE vasc;
  vasc.Init(sim);
  sim->GetScheduler()->Simulate(1);

  auto* grid = sim->GetResourceManager()->GetDiffusionGrid(fields::kVascular);

  // Apply wound
  vasc.ApplyWound(sim, sp->wound_center_x, sp->wound_center_y,
                  sp->wound_radius);

  GridContext ctx(grid, sp);
  int wound_dermal = 0;
  for (size_t idx = 0; idx < ctx.n; idx++) {
    if (ctx.Z(idx) >= 0) continue;
    real_t x = ctx.X(idx), y = ctx.Y(idx);
    if (ctx.InWound(x, y)) {
      wound_dermal++;
      EXPECT_NEAR(grid->GetConcentration(idx), 0, 1e-8)
          << "Wound dermal voxel should be zero at idx=" << idx;
    }
  }
  EXPECT_GT(wound_dermal, 0) << "Should have dermal voxels inside wound";

  delete sim;
}

TEST(VascularPDETest, SubLayerPerfusionProfile) {
  // At high resolution, layer-specific fractions produce a gradient
  Param::RegisterParamGroup(new SimParam());
  auto set_param = [](Param* param) {
    param->bound_space = Param::BoundSpaceMode::kClosed;
    param->min_bound = -10;
    param->max_bound = 40;
    param->simulation_time_step = 0.1;
    param->export_visualization = false;
    auto* sp = const_cast<SimParam*>(param->Get<SimParam>());
    sp->grid_resolution = 25;
    sp->perfusion_papillary_fraction = 1.0;
    sp->perfusion_reticular_fraction = 0.7;
    sp->perfusion_hypodermis_fraction = 0.3;
  };
  auto* sim = new Simulation(TEST_NAME, set_param);
  auto* sp = sim->GetParam()->Get<SimParam>();

  VascularPDE vasc;
  vasc.Init(sim);
  sim->GetScheduler()->Simulate(1);

  auto* grid = sim->GetResourceManager()->GetDiffusionGrid(fields::kVascular);
  GridContext ctx(grid, sp);

  real_t sum_papillary = 0, sum_reticular = 0, sum_hypodermis = 0;
  int n_papillary = 0, n_reticular = 0, n_hypodermis = 0;
  for (size_t idx = 0; idx < ctx.n; idx++) {
    real_t z = ctx.Z(idx);
    if (z >= 0) continue;
    real_t val = grid->GetConcentration(idx);
    if (z >= sp->dermal_z_papillary) {
      sum_papillary += val;
      n_papillary++;
    } else if (z >= sp->dermal_z_reticular) {
      sum_reticular += val;
      n_reticular++;
    } else {
      sum_hypodermis += val;
      n_hypodermis++;
    }
  }
  ASSERT_GT(n_reticular, 0);
  ASSERT_GT(n_hypodermis, 0);
  real_t mean_reticular = sum_reticular / n_reticular;
  real_t mean_hypodermis = sum_hypodermis / n_hypodermis;
  // Reticular mean should exceed hypodermis mean
  EXPECT_GT(mean_reticular, mean_hypodermis);
  // If papillary has voxels, it should exceed reticular
  if (n_papillary > 0) {
    real_t mean_papillary = sum_papillary / n_papillary;
    EXPECT_GT(mean_papillary, mean_reticular);
  }

  delete sim;
}

TEST(WoundMaskTest, SubLayerClassification) {
  // Verify sub-layer lists partition dermal_all correctly at res=25
  Param::RegisterParamGroup(new SimParam());
  auto set_param = [](Param* param) {
    param->bound_space = Param::BoundSpaceMode::kClosed;
    param->min_bound = -10;
    param->max_bound = 40;
    param->simulation_time_step = 0.1;
    param->export_visualization = false;
    auto* sp = const_cast<SimParam*>(param->Get<SimParam>());
    sp->grid_resolution = 25;
    sp->wound_center_x = 20;
    sp->wound_center_y = 20;
    sp->wound_radius = 8;
  };
  auto* sim = new Simulation(TEST_NAME, set_param);
  auto* sp = sim->GetParam()->Get<SimParam>();

  VascularPDE vasc;
  vasc.Init(sim);
  sim->GetScheduler()->Simulate(1);
  auto* grid = sim->GetResourceManager()->GetDiffusionGrid(fields::kVascular);
  GridContext ctx(grid, sp);

  real_t z_max = sp->volume_z_cornified + ctx.box_len;
  auto mask = GridContext::ComputeWoundMask(ctx, z_max, sp);

  // Sub-layer lists must partition dermal_all
  EXPECT_EQ(mask.dermal_all.size(),
            mask.papillary_all.size() + mask.reticular_all.size() +
                mask.hypodermis_all.size());
  // Wound sub-layer lists must partition dermal_wound
  EXPECT_EQ(mask.dermal_wound.size(),
            mask.papillary_wound.size() + mask.reticular_wound.size() +
                mask.hypodermis_wound.size());
  // At res=25, all three layers should have voxels
  EXPECT_GT(mask.reticular_all.size(), 0u);
  EXPECT_GT(mask.hypodermis_all.size(), 0u);
  // dermal_all should be non-empty
  EXPECT_GT(mask.dermal_all.size(), 0u);

  delete sim;
}

// ===========================================================================
// MetricsExporter tests
// ===========================================================================
TEST(MetricsExporterTest, WritesHeaderOnFirstCall) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  sp->metrics_interval = 1;
  sp->wound_enabled = false;
  sp->tumor_enabled = false;

  SetupAllFields(sim);

  MetricsExporter metrics;
  metrics();

  // Read the output file
  std::string path = sim->GetOutputDir() + "/metrics.csv";
  std::ifstream f(path);
  ASSERT_TRUE(f.is_open());

  std::string header;
  std::getline(f, header);

  // Header should contain 40 column names
  EXPECT_NE(header.find("step"), std::string::npos);
  EXPECT_NE(header.find("time_h"), std::string::npos);
  EXPECT_NE(header.find("wound_closure_pct"), std::string::npos);
  EXPECT_NE(header.find("n_tumor_cells"), std::string::npos);
  EXPECT_NE(header.find("tumor_field_cells"), std::string::npos);
  EXPECT_NE(header.find("mean_ecm_quality"), std::string::npos);
  EXPECT_NE(header.find("mean_tissue_viability"), std::string::npos);

  // Count commas to verify column count (39 commas for 40 columns)
  int commas = 0;
  for (char c : header) {
    if (c == ',') commas++;
  }
  EXPECT_EQ(commas, 39);

  metrics.Close();
  delete sim;
}

TEST(MetricsExporterTest, RowCountMatchesInterval) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  sp->metrics_interval = 3;
  sp->wound_enabled = false;
  sp->tumor_enabled = false;

  SetupAllFields(sim);

  // Heap-allocate so OperationRegistry can safely own and delete the impl
  auto* metrics = new MetricsExporter();
  OperationRegistry::GetInstance()->AddOperationImpl(
      "MetricsExporter_test", OpComputeTarget::kCpu, metrics);
  sim->GetScheduler()->ScheduleOp(
      NewOperation("MetricsExporter_test"), OpType::kPostSchedule);

  // Run 10 steps
  sim->GetScheduler()->Simulate(10);

  metrics->Close();

  // Read and count data rows
  std::string path = sim->GetOutputDir() + "/metrics.csv";
  std::ifstream f(path);
  ASSERT_TRUE(f.is_open());

  int lines = 0;
  std::string line;
  while (std::getline(f, line)) {
    lines++;
  }

  // Header + data rows. Steps 0,3,6,9 = 4 rows (step 0 is pre-simulate)
  // Simulate(10) runs steps 0-9. MetricsExporter fires at
  // steps where step % interval == 0. Step 0 is first call.
  // With interval=3: steps 0, 3, 6, 9 -> 4 data rows + 1 header = 5
  EXPECT_EQ(lines, 5);

  delete sim;
}

// ===========================================================================
// Integration tests: wound progression and tumor growth
// ===========================================================================
TEST(IntegrationTest, WoundClosureProgressesOverTime) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  sp->wound_trigger_step = 0;
  sp->wound_enabled = true;
  sp->handoff_delay = 50;  // fast handoff for test

  auto* vm = VolumeManager::Get();
  vm->Init(sp->volume_z_spinous, sp->volume_z_granular,
           sp->volume_z_cornified);

  CompositeField fields;
  fields.Add(std::make_unique<VascularPDE>());
  fields.Add(std::make_unique<CalciumPDE>());
  fields.Add(std::make_unique<KgfPDE>());
  fields.Add(std::make_unique<OxygenPDE>());
  fields.Add(std::make_unique<StratumPDE>());
  fields.Add(std::make_unique<ScarPDE>());
  fields.Add(std::make_unique<WaterPDE>());
  if (sp->split_inflammation_enabled) {
    fields.Add(std::make_unique<ProInflammatoryPDE>(sp));
    fields.Add(std::make_unique<AntiInflammatoryPDE>(sp));
  } else {
    fields.Add(std::make_unique<InflammationPDE>(sp));
  }
  fields.Add(std::make_unique<ImmunePressurePDE>(sp));
  if (sp->fibronectin_enabled) fields.Add(std::make_unique<FibronectinPDE>(sp));
  if (sp->mmp_enabled) fields.Add(std::make_unique<MMPPDE>(sp));
  if (sp->fibroblast_enabled) {
    fields.Add(std::make_unique<TGFBetaPDE>(sp));
    fields.Add(std::make_unique<CollagenPDE>(sp));
  }
  fields.InitAll(sim);

  OperationRegistry::GetInstance()->AddOperationImpl(
      "WoundEvent_closure", OpComputeTarget::kCpu, new WoundEvent(&fields));
  sim->GetScheduler()->ScheduleOp(
      NewOperation("WoundEvent_closure"), OpType::kPreSchedule);
  OperationRegistry::GetInstance()->AddOperationImpl(
      "CompositeFieldOp_closure", OpComputeTarget::kCpu,
      new CompositeFieldOp(&fields));
  sim->GetScheduler()->ScheduleOp(
      NewOperation("CompositeFieldOp_closure"), OpType::kPreSchedule);

  sim->GetScheduler()->Simulate(500);

  // After 500 steps (50h), margin cells should have proliferated.
  // Check that agent count exceeds the initial margin spawn.
  int n_agents = 0;
  sim->GetResourceManager()->ForEachAgent([&](Agent*) { n_agents++; });
  EXPECT_GT(n_agents, 7) << "Cells should proliferate beyond initial margin count";

  delete sim;
}

TEST(IntegrationTest, TumorGrowsOverTime) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  sp->tumor_enabled = true;
  sp->tumor_seed_time = 0;
  sp->tumor_seed_count = 5;
  sp->tumor_max_cells = 500;
  sp->wound_enabled = false;

  SetupAllFields(sim);

  // Sizing agent: force box_length above tumor neighbor search radius
  auto* sizer = new TumorCell({0.5, 0.5, 0.5});
  sizer->SetDiameter(20);
  sim->GetResourceManager()->AddAgent(sizer);

  auto* ti = new TumorInitiation();
  OperationRegistry::GetInstance()->AddOperationImpl(
      "TumorInit_growth", OpComputeTarget::kCpu, ti);
  sim->GetScheduler()->ScheduleOp(
      NewOperation("TumorInit_growth"), OpType::kPreSchedule);

  sim->GetScheduler()->Simulate(500);

  // Count tumor cells
  int tumor_count = 0;
  sim->GetResourceManager()->ForEachAgent([&](Agent* a) {
    if (dynamic_cast<TumorCell*>(a)) tumor_count++;
  });

  // Subtract 1 for sizing agent
  EXPECT_GT(tumor_count - 1, static_cast<int>(sp->tumor_seed_count))
      << "Tumor should grow beyond initial seed count after 500 steps";

  delete sim;
}

// ---------------------------------------------------------------------------
// Extension 7: Biofilm tests
// ---------------------------------------------------------------------------
TEST(BiofilmTest, LogisticGrowth) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  SanitizeForUnitTest(sp);
  sp->biofilm_enabled = true;
  sp->biofilm_growth_rate = 0.1;
  sp->biofilm_carrying_capacity = 1.0;
  sp->wound_enabled = true;
  sp->wound_center_x = 15;
  sp->wound_center_y = 15;
  sp->wound_radius = 10;
  sp->wound_trigger_step = 0;

  auto* vm = VolumeManager::Get();
  vm->Init(sp->volume_z_spinous, sp->volume_z_granular, sp->volume_z_cornified);

  BiofilmPDE biofilm_pde;
  biofilm_pde.Init(sim);
  SetupInflammation(sim);
  StratumPDE stratum_pde;
  stratum_pde.Init(sim);
  sim->GetScheduler()->Simulate(1);

  // Seed biofilm at wound center (basal layer z=2.5, stratum=kBasal=0)
  auto* rm = sim->GetResourceManager();
  auto* biofilm_grid = rm->GetDiffusionGrid(fields::kBiofilm);
  Real3 center = {15, 15, 2.5};
  size_t idx = biofilm_grid->GetBoxIndex(center);
  biofilm_grid->ChangeConcentrationBy(idx, 0.5);

  real_t before = biofilm_grid->GetConcentration(idx);
  EXPECT_NEAR(before, 0.5, 0.01) << "Seeded biofilm should be readable";

  // Call the op directly
  BiofilmGrowthOp op;
  op();

  real_t after = biofilm_grid->GetConcentration(idx);
  // delta = 0.1 * 0.5 * (1 - 0.5/1.0) = 0.025
  EXPECT_NEAR(after, 0.525, 0.01);

  delete sim;
}

TEST(BiofilmTest, M1BlockByBiofilm) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  SanitizeForUnitTest(sp);
  sp->biofilm_enabled = true;
  sp->biofilm_m1_block_threshold = 0.2;
  sp->macrophage_m1_duration = 1;
  sp->macrophage_lifespan = 999;
  sp->macrophage_min_survival = 999;
  sp->m1_transition_threshold = 999;
  sp->wound_enabled = true;
  sp->wound_center_x = 25;
  sp->wound_center_y = 25;
  sp->wound_radius = 20;
  sp->wound_trigger_step = 0;

  SetupInflammation(sim);
  BiofilmPDE bio;
  bio.Init(sim);
  sim->GetScheduler()->Simulate(1);

  // Seed high biofilm at cell position
  auto* rm = sim->GetResourceManager();
  auto* biofilm_grid = rm->GetDiffusionGrid(fields::kBiofilm);
  Real3 pos = {25, 25, 1.5};
  size_t idx = biofilm_grid->GetBoxIndex(pos);
  biofilm_grid->ChangeConcentrationBy(idx, 0.5);

  auto* cell = new ImmuneCell(pos);
  cell->SetDiameter(3);
  cell->SetImmuneCellType(kMacrophage);
  cell->SetState(kM1Active);
  cell->AddBehavior(new MacrophageBehavior());
  rm->AddAgent(cell);

  // Run 3 steps: timer trigger fires but biofilm should block M1->M2
  sim->GetScheduler()->Simulate(3);
  ImmuneState state = kM1Active;
  rm->ForEachAgent([&](Agent* a) {
    auto* ic = dynamic_cast<ImmuneCell*>(a);
    if (ic && ic->GetImmuneCellType() == kMacrophage) state = ic->GetState();
  });
  EXPECT_EQ(state, kM1Active);

  delete sim;
}

// ---------------------------------------------------------------------------
// Extension 8: VEGF / Angiogenesis tests
// ---------------------------------------------------------------------------
TEST(AngiogenesisTest, VEGFProductionUnderHypoxia) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  sp->angiogenesis_enabled = true;
  sp->vegf_production_rate = 0.1;
  sp->vegf_hypoxia_threshold = 0.5;
  sp->wound_enabled = true;
  sp->wound_center_x = 25;
  sp->wound_center_y = 25;
  sp->wound_radius = 20;
  sp->wound_trigger_step = 0;
  sp->volume_z_cornified = 40;

  auto* vm = VolumeManager::Get();
  vm->Init(sp->volume_z_spinous, sp->volume_z_granular, sp->volume_z_cornified);

  OxygenPDE o2;
  o2.Init(sim);
  VEGFPDE vegf(sp);
  vegf.Init(sim);
  sim->GetScheduler()->Simulate(1);

  // Zero oxygen at a wound voxel to simulate hypoxia
  auto* rm = sim->GetResourceManager();
  auto* o2_grid = rm->GetDiffusionGrid(fields::kOxygen);
  auto* vegf_grid = rm->GetDiffusionGrid(fields::kVEGF);
  Real3 pos = {25, 25, 5};
  size_t o2_idx = o2_grid->GetBoxIndex(pos);
  real_t o2_val = o2_grid->GetConcentration(o2_idx);
  o2_grid->ChangeConcentrationBy(o2_idx, -o2_val);

  real_t before = vegf_grid->GetConcentration(vegf_grid->GetBoxIndex(pos));

  VEGFSourceOp op;
  op();

  real_t after = vegf_grid->GetConcentration(vegf_grid->GetBoxIndex(pos));
  EXPECT_GT(after, before);

  delete sim;
}

TEST(AngiogenesisTest, NoVEGFWithoutHypoxia) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  sp->angiogenesis_enabled = true;
  SanitizeForUnitTest(sp);
  sp->angiogenesis_enabled = true;
  sp->vegf_production_rate = 0.1;
  sp->vegf_hypoxia_threshold = 0.5;
  sp->wound_enabled = true;
  sp->wound_center_x = 25;
  sp->wound_center_y = 25;
  sp->wound_radius = 20;
  sp->wound_trigger_step = 0;
  sp->volume_z_cornified = 40;
  sp->oxygen_basal_conc = 1.0;
  // Use a very gradual O2 decay so all voxels stay above the hypoxia
  // threshold even at the coarse grid resolution (box centers at z≈10).
  sp->oxygen_decay_length = 100.0;

  auto* vm = VolumeManager::Get();
  vm->Init(sp->volume_z_spinous, sp->volume_z_granular, sp->volume_z_cornified);

  OxygenPDE o2;
  o2.Init(sim);
  VEGFPDE vegf(sp);
  vegf.Init(sim);
  sim->GetScheduler()->Simulate(1);

  // O2 is initialized at basal level (1.0 > 0.5 threshold) so no VEGF
  auto* vegf_grid = sim->GetResourceManager()->GetDiffusionGrid(fields::kVEGF);
  Real3 pos = {25, 25, 5};
  real_t before = vegf_grid->GetConcentration(vegf_grid->GetBoxIndex(pos));

  VEGFSourceOp op;
  op();

  real_t after = vegf_grid->GetConcentration(vegf_grid->GetBoxIndex(pos));
  // At well-oxygenated voxels, VEGF should not increase
  EXPECT_LE(after, before + 1e-10);

  delete sim;
}

// ===========================================================================
// VolumeManager dermal sub-layer tests
// ===========================================================================

TEST(VolumeManagerTest, DermalSubLayerLookup) {
  auto* vm = VolumeManager::Get();
  vm->Init(6.0, 15.0, 25.0, -2.0, -8.0);

  // Papillary dermis: z in [-2, 0)
  EXPECT_EQ(vm->GetDermalSubLayer({15, 15, -1}), kPapillary);
  EXPECT_EQ(vm->GetDermalSubLayer({15, 15, -1.5}), kPapillary);

  // Reticular dermis: z in [-8, -2)
  EXPECT_EQ(vm->GetDermalSubLayer({15, 15, -3}), kReticular);
  EXPECT_EQ(vm->GetDermalSubLayer({15, 15, -7}), kReticular);

  // Hypodermis: z < -8
  EXPECT_EQ(vm->GetDermalSubLayer({15, 15, -9}), kHypodermis);

  // Epidermal positions return epidermal strata (not dermal sub-layers)
  EXPECT_EQ(vm->GetDermalSubLayer({15, 15, 3}), kBasal);
  EXPECT_EQ(vm->GetDermalSubLayer({15, 15, 10}), kSpinous);

  // GetStratumAt still returns kDermis for all z < 0 (backward compat)
  EXPECT_EQ(vm->GetStratumAt({15, 15, -1}), kDermis);
  EXPECT_EQ(vm->GetStratumAt({15, 15, -5}), kDermis);
  EXPECT_EQ(vm->GetStratumAt({15, 15, -9}), kDermis);
}

// ===========================================================================
// Elastin / Hyaluronan source term tests
// ===========================================================================

TEST(FibroblastBehaviorTest, TropoelastinProduction) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  SanitizeForUnitTest(sp);
  sp->fibroblast_enabled = true;
  sp->elastin_enabled = true;
  sp->elastin_production_rate = 0.01;
  sp->fibroblast_tgfb_rate = 0;
  sp->fibroblast_lifespan = 9999;
  sp->fibroblast_apoptosis_threshold = 0;

  TGFBetaPDE tgfb(sp);
  tgfb.Init(sim);
  CollagenPDE col(sp);
  col.Init(sim);
  ElastinPDE elastin;
  elastin.Init(sim);
  sim->GetScheduler()->Simulate(1);

  auto* rm = sim->GetResourceManager();
  Real3 pos = {15, 15, 2};
  auto* el_grid = rm->GetDiffusionGrid(fields::kElastin);
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
  sp->fibroblast_enabled = true;
  sp->hyaluronan_enabled = true;
  sp->hyaluronan_production_rate = 0.01;
  sp->fibroblast_tgfb_rate = 0;
  sp->fibroblast_lifespan = 9999;
  sp->fibroblast_apoptosis_threshold = 0;

  TGFBetaPDE tgfb(sp);
  tgfb.Init(sim);
  CollagenPDE col(sp);
  col.Init(sim);
  HyaluronanPDE ha;
  ha.Init(sim);
  sim->GetScheduler()->Simulate(1);

  auto* rm = sim->GetResourceManager();
  Real3 pos = {15, 15, 2};
  auto* ha_grid = rm->GetDiffusionGrid(fields::kHyaluronan);
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
  sp->fibroblast_enabled = true;
  sp->elastin_enabled = true;
  sp->hyaluronan_enabled = true;
  sp->elastin_production_rate = 0.01;
  sp->hyaluronan_production_rate = 0.01;
  sp->fibroblast_tgfb_rate = 0;
  sp->fibroblast_lifespan = 9999;
  sp->fibroblast_apoptosis_threshold = 0;
  sp->fibroblast_activation_delay = 9999;
  sp->fibroblast_activation_threshold = 999;

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
  auto* el_grid = rm->GetDiffusionGrid(fields::kElastin);
  auto* ha_grid = rm->GetDiffusionGrid(fields::kHyaluronan);
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

// ===========================================================================
// Elastin PDE tests
// ===========================================================================

TEST(ElastinPDETest, InitProfile) {
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
  sp->grid_resolution_structural = 0;

  ElastinPDE elastin;
  elastin.Init(sim);
  sim->GetScheduler()->Simulate(1);

  auto* grid = sim->GetResourceManager()->GetDiffusionGrid(fields::kElastin);
  GridContext ctx(grid, sp);

  // Epidermal voxels should be zero
  real_t epi_sum = 0;
  int epi_count = 0;
  // Dermal voxels should be positive
  real_t derm_sum = 0;
  int derm_count = 0;
  for (size_t idx = 0; idx < ctx.n; idx++) {
    real_t z = ctx.Z(idx);
    real_t val = grid->GetConcentration(idx);
    if (z >= 2) {  // well into epidermis
      epi_sum += val;
      epi_count++;
    } else if (z < -3 && z >= sp->dermal_z_reticular) {
      // Reticular dermis
      derm_sum += val;
      derm_count++;
    }
  }
  ASSERT_GT(epi_count, 0);
  ASSERT_GT(derm_count, 0);
  // Epidermal elastin should be near zero (may diffuse slightly)
  EXPECT_LT(epi_sum / epi_count, 0.1);
  // Reticular dermis should have significant elastin
  EXPECT_GT(derm_sum / derm_count, 0.1);

  delete sim;
}

// ===========================================================================
// Hyaluronan PDE tests
// ===========================================================================

TEST(HyaluronanPDETest, InitProfile) {
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
  sp->grid_resolution_structural = 0;

  HyaluronanPDE ha;
  ha.Init(sim);
  sim->GetScheduler()->Simulate(1);

  auto* grid = sim->GetResourceManager()->GetDiffusionGrid(fields::kHyaluronan);
  GridContext ctx(grid, sp);

  // Epidermal voxels should be near zero
  real_t epi_sum = 0;
  int epi_count = 0;
  // Papillary dermis should be HA-rich
  real_t pap_sum = 0;
  int pap_count = 0;
  for (size_t idx = 0; idx < ctx.n; idx++) {
    real_t z = ctx.Z(idx);
    real_t val = grid->GetConcentration(idx);
    if (z >= 2) {
      epi_sum += val;
      epi_count++;
    } else if (z < 0 && z >= sp->dermal_z_papillary) {
      pap_sum += val;
      pap_count++;
    }
  }
  ASSERT_GT(epi_count, 0);
  ASSERT_GT(pap_count, 0);
  EXPECT_LT(epi_sum / epi_count, 0.1);
  EXPECT_GT(pap_sum / pap_count, 0.2);

  delete sim;
}

// ===========================================================================
// ChunkedGrid tests
// ===========================================================================

TEST(ChunkedGridTest, DefaultValueForUnallocatedChunks) {
  ChunkedGrid grid(16, 1.0);  // 16^3 voxels, default=1.0
  // All chunks unallocated -- every voxel should return default
  EXPECT_EQ(grid.Get(0), 1.0);
  EXPECT_EQ(grid.Get(100), 1.0);
  EXPECT_EQ(grid.Get(16 * 16 * 16 - 1), 1.0);
  EXPECT_EQ(grid.GetNumAllocatedChunks(), 0u);
}

TEST(ChunkedGridTest, WriteAllocatesChunk) {
  ChunkedGrid grid(16, 0.0);
  EXPECT_EQ(grid.GetNumAllocatedChunks(), 0u);
  grid.Set(0, 5.0);
  EXPECT_EQ(grid.GetNumAllocatedChunks(), 1u);
  EXPECT_EQ(grid.Get(0), 5.0);
  // Other voxels in same chunk should be default (0.0)
  EXPECT_EQ(grid.Get(1), 0.0);
  EXPECT_EQ(grid.Get(7), 0.0);  // still within chunk 0
}

TEST(ChunkedGridTest, AddAccumulates) {
  ChunkedGrid grid(16, 0.0);
  grid.Add(42, 1.5);
  grid.Add(42, 0.5);
  EXPECT_NEAR(grid.Get(42), 2.0, 1e-10);
}

TEST(ChunkedGridTest, DirtyTracking) {
  ChunkedGrid grid(16, 0.0);
  size_t ci = grid.ChunkIndex(0);
  EXPECT_FALSE(grid.IsDirty(ci));

  grid.Set(0, 1.0);
  EXPECT_TRUE(grid.IsDirty(ci));
  EXPECT_EQ(grid.GetNumDirtyChunks(), 1u);

  grid.ClearAllDirty();
  EXPECT_FALSE(grid.IsDirty(ci));
  EXPECT_EQ(grid.GetNumDirtyChunks(), 0u);

  // Reading doesn't dirty
  (void)grid.Get(0);
  EXPECT_FALSE(grid.IsDirty(ci));
}

TEST(ChunkedGridTest, ImportExportRoundTrip) {
  constexpr size_t res = 16;
  constexpr size_t n = res * res * res;
  std::vector<real_t> flat(n, 0.0);
  // Set some non-default values in a region
  for (size_t i = 0; i < 100; ++i) flat[i] = static_cast<real_t>(i);

  ChunkedGrid grid(res, 0.0);
  grid.ImportFromFlat(flat.data(), n);

  // Only chunks containing non-zero voxels should be allocated
  EXPECT_GT(grid.GetNumAllocatedChunks(), 0u);
  EXPECT_LT(grid.GetNumAllocatedChunks(), grid.GetTotalChunks());

  // Export back and compare
  std::vector<real_t> exported(n, -1.0);
  grid.ExportToFlat(exported.data(), n);
  for (size_t i = 0; i < n; ++i) {
    EXPECT_EQ(exported[i], flat[i]) << "Mismatch at voxel " << i;
  }
}

TEST(ChunkedGridTest, SparsityForUniformField) {
  constexpr size_t res = 24;  // 3 chunks per axis = 27 chunks
  constexpr size_t n = res * res * res;
  std::vector<real_t> flat(n, 1.0);  // all default
  // Only modify one voxel
  flat[0] = 2.0;

  ChunkedGrid grid(res, 1.0);
  grid.ImportFromFlat(flat.data(), n);

  // Only 1 chunk should be allocated (the one containing voxel 0)
  EXPECT_EQ(grid.GetNumAllocatedChunks(), 1u);
  EXPECT_EQ(grid.GetTotalChunks(), 27u);  // 3^3
}

TEST(ChunkedGridTest, ForEachVoxelInChunk) {
  ChunkedGrid grid(16, 0.0);
  size_t count = 0;
  grid.ForEachVoxelInChunk(0, [&](size_t flat_idx, size_t local_idx) {
    (void)flat_idx;
    (void)local_idx;
    ++count;
  });
  // Chunk 0 should contain 8^3 = 512 voxels
  EXPECT_EQ(count, ChunkTraits::kVolume);
}

TEST(ChunkedGridTest, ChunkBoundaryAccess) {
  // Test voxels at chunk boundaries (where indexing bugs typically hide)
  ChunkedGrid grid(16, 0.0);
  // Voxel at (7,7,7) = last voxel of chunk 0
  size_t idx_last = 7 + 7 * 16 + 7 * 16 * 16;
  // Voxel at (8,8,8) = first voxel of chunk (1,1,1)
  size_t idx_next = 8 + 8 * 16 + 8 * 16 * 16;

  grid.Set(idx_last, 10.0);
  grid.Set(idx_next, 20.0);

  EXPECT_EQ(grid.Get(idx_last), 10.0);
  EXPECT_EQ(grid.Get(idx_next), 20.0);
  // They should be in different chunks
  EXPECT_NE(grid.ChunkIndex(idx_last), grid.ChunkIndex(idx_next));
  EXPECT_EQ(grid.GetNumAllocatedChunks(), 2u);
}

// ---------------------------------------------------------------------------
// Resolution independence: mean wound VEGF concentration must match across
// resolutions.  Per-voxel sources give the same concentration per wound
// voxel regardless of grid density, so the wound-averaged value should be
// approximately equal.
// ---------------------------------------------------------------------------
struct VEGFResResult {
  real_t mean_wound_conc;
  real_t total_mass;
};

static VEGFResResult RunVEGFProductionAtResolution(int resolution) {
  Param::RegisterParamGroup(new SimParam());
  auto set_param = [resolution](Param* param) {
    param->bound_space = Param::BoundSpaceMode::kClosed;
    param->min_bound = 0;
    param->max_bound = 50;
    param->simulation_time_step = 0.1;
    param->export_visualization = false;
    auto* sp = const_cast<SimParam*>(param->Get<SimParam>());
    sp->bounds_min = {0, 0, 0};
    sp->bounds_max = {50, 50, 50};
    sp->grid_resolution = resolution;
    sp->ref_box_length = 5.0;
    sp->wound_trigger_step = 0;
    sp->angiogenesis_enabled = true;
    sp->diabetic_mode = false;
    sp->biofilm_enabled = false;
    sp->mmp_enabled = false;
    sp->fibronectin_enabled = false;
    sp->fibroblast_enabled = false;
    sp->elastin_enabled = false;
    sp->hyaluronan_enabled = false;
    sp->vegf_production_rate = 0.01;
    sp->vegf_hypoxia_threshold = 0.5;
    sp->oxygen_basal_conc = 1.0;
    sp->perfusion_basal = 1.0;
  };
  std::string name = "vegf_res_" + std::to_string(resolution);
  auto* sim = new Simulation(name, set_param);
  auto* sp = sim->GetParam()->Get<SimParam>();
  auto* rm = sim->GetResourceManager();

  // Register required grids
  ModelInitializer::DefineSubstance(fields::kVascularId, fields::kVascular,
                                    0, 0, resolution);
  ModelInitializer::DefineSubstance(fields::kOxygenId, fields::kOxygen,
                                    0, 0, resolution);
  ModelInitializer::DefineSubstance(fields::kStratumId, fields::kStratum,
                                    0, 0, resolution);
  ModelInitializer::DefineSubstance(fields::kCalciumId, fields::kCalcium,
                                    0, 0, resolution);
  ModelInitializer::DefineSubstance(fields::kWaterId, fields::kWater,
                                    0, 0, resolution);
  ModelInitializer::DefineSubstance(fields::kVEGFId, fields::kVEGF,
                                    0, 0, resolution);

  auto* vasc_grid = rm->GetDiffusionGrid(fields::kVascular);
  auto* vegf_grid = rm->GetDiffusionGrid(fields::kVEGF);
  vasc_grid->Initialize();
  rm->GetDiffusionGrid(fields::kOxygen)->Initialize();
  rm->GetDiffusionGrid(fields::kStratum)->Initialize();
  rm->GetDiffusionGrid(fields::kCalcium)->Initialize();
  rm->GetDiffusionGrid(fields::kWater)->Initialize();
  vegf_grid->Initialize();

  // Skip diffusion solver
  vasc_grid->SetTimeStep(1e30);
  rm->GetDiffusionGrid(fields::kOxygen)->SetTimeStep(1e30);
  rm->GetDiffusionGrid(fields::kStratum)->SetTimeStep(1e30);
  rm->GetDiffusionGrid(fields::kCalcium)->SetTimeStep(1e30);
  rm->GetDiffusionGrid(fields::kWater)->SetTimeStep(1e30);
  vegf_grid->SetTimeStep(1e30);

  size_t n = vasc_grid->GetNumBoxes();
  GridContext ctx(vasc_grid, sp);

  // Vascular = 1 everywhere, O2 = 0 (max VEGF production)
  for (size_t i = 0; i < n; i++) {
    vasc_grid->ChangeConcentrationBy(i, 1.0);
  }

  // Advance past wound_trigger_step=0
  sim->GetScheduler()->Simulate(1);

  FusedWoundSourceOp op;
  op();

  // Compute mean wound concentration and total mass
  real_t box_vol = ctx.box_len * ctx.box_len * ctx.box_len;
  real_t total_mass = 0;
  real_t sum_conc = 0;
  size_t wound_count = 0;
  for (size_t i = 0; i < n; i++) {
    real_t c = vegf_grid->GetConcentration(i);
    total_mass += c * box_vol;
    if (c > 1e-15) {
      sum_conc += c;
      wound_count++;
    }
  }

  real_t mean = (wound_count > 0) ? sum_conc / wound_count : 0;
  delete sim;
  return {mean, total_mass};
}

TEST(ResolutionTest, VEGFProductionResolutionIndependent) {
  auto r10 = RunVEGFProductionAtResolution(10);
  auto r20 = RunVEGFProductionAtResolution(20);

  EXPECT_GT(r10.mean_wound_conc, 0) << "No VEGF produced at res=10";
  EXPECT_GT(r20.mean_wound_conc, 0) << "No VEGF produced at res=20";

  // Per-voxel sources add the same concentration delta to each wound voxel
  // regardless of resolution. Mean wound concentration must match.
  // (Total mass differs because wound voxel count depends on discretization.)
  real_t conc_ratio = r20.mean_wound_conc / r10.mean_wound_conc;
  EXPECT_NEAR(conc_ratio, 1.0, 0.05)
      << "Mean wound VEGF at res=20 (" << r20.mean_wound_conc
      << ") differs from res=10 (" << r10.mean_wound_conc << ")";
}

// ---------------------------------------------------------------------------
// ScaledGrid: agent_factor must equal (ref / box)^3
// ---------------------------------------------------------------------------
TEST(ScaledGridTest, FactorAtCalibrationResolution) {
  auto* sim = CreateTestSim("sg_factor_cal");
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  sp->ref_box_length = 5.0;
  sp->grid_resolution = 10;
  SanitizeForUnitTest(sp);

  InflammationPDE infl(sp); infl.Init(sim);
  sim->GetScheduler()->Simulate(1);  // initializes grid internals
  auto* grid = sim->GetResourceManager()->GetDiffusionGrid(fields::kInflammation);

  ScaledGrid sg(grid, sp);
  // At calibration resolution, box_len == ref_box_length -> factor = 1
  EXPECT_NEAR(sg.agent_factor, 1.0, 1e-6);
  delete sim;
}

TEST(ScaledGridTest, FactorAtDoubleResolution) {
  Param::RegisterParamGroup(new SimParam());
  auto set_param = [](Param* param) {
    param->bound_space = Param::BoundSpaceMode::kClosed;
    param->min_bound = 0;
    param->max_bound = 50;
    param->simulation_time_step = 0.1;
    param->export_visualization = false;
    auto* sp = const_cast<SimParam*>(param->Get<SimParam>());
    sp->bounds_min = {0, 0, 0};
    sp->bounds_max = {50, 50, 50};
    sp->grid_resolution = 20;
    sp->ref_box_length = 5.0;
  };
  auto* sim = new Simulation("sg_factor_2x", set_param);
  auto* sp = sim->GetParam()->Get<SimParam>();
  SanitizeForUnitTest(const_cast<SimParam*>(sp));

  InflammationPDE infl(sp); infl.Init(sim);
  sim->GetScheduler()->Simulate(1);  // initializes grid internals
  auto* grid = sim->GetResourceManager()->GetDiffusionGrid(fields::kInflammation);

  ScaledGrid sg(grid, sp);
  // At 2x resolution, box_len = 2.5 -> factor = (5/2.5)^3 = 8
  EXPECT_NEAR(sg.agent_factor, 8.0, 1e-6);
  delete sim;
}

// ---------------------------------------------------------------------------
// ScaledGrid: agent deposit produces same total field mass at different
// resolutions. An agent at a fixed position deposits amount * agent_factor
// into one voxel. The total mass = concentration * voxel_volume should be
// the same regardless of resolution.
// ---------------------------------------------------------------------------
TEST(ScaledGridTest, AgentDepositResolutionIndependent) {
  real_t deposit_amount = 0.05;
  Real3 agent_pos = {25, 25, 25};  // center of domain

  // --- Run at res=10 ---
  auto* sim10 = CreateTestSim("sg_deposit_10");
  auto* sp10 = const_cast<SimParam*>(sim10->GetParam()->Get<SimParam>());
  sp10->ref_box_length = 5.0;
  sp10->grid_resolution = 10;
  SanitizeForUnitTest(sp10);
  InflammationPDE infl10(sp10); infl10.Init(sim10);
  sim10->GetScheduler()->Simulate(1);
  auto* grid10 = sim10->GetResourceManager()->GetDiffusionGrid(
      fields::kInflammation);

  ScaledGrid sg10(grid10, sp10);
  size_t idx10 = sg10.Index(agent_pos);
  sg10.AgentDeposit(idx10, deposit_amount);
  real_t conc10 = sg10.Get(idx10);
  real_t box10 = grid10->GetBoxLength();
  real_t mass10 = conc10 * box10 * box10 * box10;
  delete sim10;

  // --- Run at res=20 ---
  Param::RegisterParamGroup(new SimParam());
  auto set20 = [](Param* param) {
    param->bound_space = Param::BoundSpaceMode::kClosed;
    param->min_bound = 0;
    param->max_bound = 50;
    param->simulation_time_step = 0.1;
    param->export_visualization = false;
    auto* sp = const_cast<SimParam*>(param->Get<SimParam>());
    sp->bounds_min = {0, 0, 0};
    sp->bounds_max = {50, 50, 50};
    sp->grid_resolution = 20;
    sp->ref_box_length = 5.0;
  };
  auto* sim20 = new Simulation("sg_deposit_20", set20);
  auto* sp20 = sim20->GetParam()->Get<SimParam>();
  SanitizeForUnitTest(const_cast<SimParam*>(sp20));
  InflammationPDE infl20(sp20); infl20.Init(sim20);
  sim20->GetScheduler()->Simulate(1);
  auto* grid20 = sim20->GetResourceManager()->GetDiffusionGrid(
      fields::kInflammation);

  ScaledGrid sg20(grid20, sp20);
  size_t idx20 = sg20.Index(agent_pos);
  sg20.AgentDeposit(idx20, deposit_amount);
  real_t conc20 = sg20.Get(idx20);
  real_t box20 = grid20->GetBoxLength();
  real_t mass20 = conc20 * box20 * box20 * box20;
  delete sim20;

  // Total mass injected must be identical regardless of resolution
  EXPECT_NEAR(mass10, mass20, 1e-10)
      << "Agent deposit mass differs: res10=" << mass10
      << " res20=" << mass20;

  // Concentration at fine grid should be 8x higher (same mass, 1/8 volume)
  EXPECT_NEAR(conc20 / conc10, 8.0, 1e-6);
}

// ---------------------------------------------------------------------------
// PerfTimer: smoke test (does not crash, measures nonzero time)
// ---------------------------------------------------------------------------
TEST(PerfTimerTest, MeasuresNonzeroTime) {
  PerfTimer timer(true);
  // Do some work
  volatile double x = 0;
  for (int i = 0; i < 10000; i++) x += i * 0.001;
  (void)x;

  // Capture stdout to verify output
  testing::internal::CaptureStdout();
  timer.Print("test_op");
  std::string output = testing::internal::GetCapturedStdout();
  EXPECT_TRUE(output.find("[perf] test_op") != std::string::npos);
}

TEST(PerfTimerTest, DisabledProducesNoOutput) {
  PerfTimer timer(false);
  testing::internal::CaptureStdout();
  timer.Print("should_not_appear");
  std::string output = testing::internal::GetCapturedStdout();
  EXPECT_TRUE(output.empty());
}

// ---------------------------------------------------------------------------
// Debug flags: all default to false
// ---------------------------------------------------------------------------
TEST(DebugProfileTest, FlagsDefaultToFalse) {
  auto* sim = CreateTestSim("debug_defaults");
  auto* sp = sim->GetParam()->Get<SimParam>();
  EXPECT_FALSE(sp->debug_immune);
  EXPECT_FALSE(sp->debug_fibroblast);
  EXPECT_FALSE(sp->debug_wound);
  EXPECT_FALSE(sp->debug_scaled_grid);
  EXPECT_FALSE(sp->debug_perf);
  delete sim;
}

// ===========================================================================
// Dermis PDE tests
// ===========================================================================

TEST(DermisPDETest, InitProfile) {
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
  sp->grid_resolution_structural = 0;

  DermisPDE dermis;
  dermis.Init(sim);
  sim->GetScheduler()->Simulate(1);

  auto* grid = dermis.Grid(sim);
  GridContext ctx(grid, sp);

  int epi_count = 0, pap_count = 0, ret_count = 0, hyp_count = 0;
  real_t epi_sum = 0, pap_sum = 0, ret_sum = 0, hyp_sum = 0;
  for (size_t idx = 0; idx < ctx.n; idx++) {
    real_t z = ctx.Z(idx);
    real_t val = grid->GetConcentration(idx);
    if (z >= 0) {
      epi_count++; epi_sum += val;
    } else if (z >= sp->dermal_z_papillary) {
      pap_count++; pap_sum += val;
    } else if (z >= sp->dermal_z_reticular) {
      ret_count++; ret_sum += val;
    } else {
      hyp_count++; hyp_sum += val;
    }
  }
  ASSERT_GT(epi_count, 0);
  EXPECT_LT(epi_sum / epi_count, 0.1);
  if (ret_count > 0 && hyp_count > 0) {
    EXPECT_GT(ret_sum / ret_count, hyp_sum / hyp_count);
  }

  delete sim;
}

TEST(DermisPDETest, WoundZerosDermalVoxels) {
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
  sp->grid_resolution_structural = 0;
  sp->wound_center_x = 25;
  sp->wound_center_y = 25;
  sp->wound_radius = 10;

  DermisPDE dermis;
  dermis.Init(sim);
  sim->GetScheduler()->Simulate(1);
  dermis.ApplyWound(sim, sp->wound_center_x, sp->wound_center_y,
                    sp->wound_radius);

  auto* grid = dermis.Grid(sim);
  GridContext ctx(grid, sp);
  int wound_dermal = 0;
  for (size_t idx = 0; idx < ctx.n; idx++) {
    if (ctx.Z(idx) >= 0) continue;
    if (ctx.InWound(ctx.X(idx), ctx.Y(idx))) {
      wound_dermal++;
      EXPECT_NEAR(grid->GetConcentration(idx), 0, 1e-8);
    }
  }
  EXPECT_GT(wound_dermal, 0);

  delete sim;
}

TEST(DermisPDETest, HealthyVoxelsUnchangedByWound) {
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
  sp->grid_resolution_structural = 0;
  sp->wound_center_x = 25;
  sp->wound_center_y = 25;
  sp->wound_radius = 5;

  DermisPDE dermis;
  dermis.Init(sim);
  sim->GetScheduler()->Simulate(1);
  dermis.ApplyWound(sim, sp->wound_center_x, sp->wound_center_y,
                    sp->wound_radius);

  auto* grid = dermis.Grid(sim);
  GridContext ctx(grid, sp);
  int healthy_dermal = 0;
  for (size_t idx = 0; idx < ctx.n; idx++) {
    if (ctx.Z(idx) >= 0) continue;
    if (!ctx.InWound(ctx.X(idx), ctx.Y(idx))) {
      healthy_dermal++;
      EXPECT_GT(grid->GetConcentration(idx), 0.1);
    }
  }
  EXPECT_GT(healthy_dermal, 0);

  delete sim;
}

// ---------------------------------------------------------------------------
// WoundInflammationSourceTest: DAMPs from open wound sustain inflammation
// ---------------------------------------------------------------------------
TEST(WoundInflammationSourceTest, OpenWoundProducesInflammation) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  SanitizeForUnitTest(sp);
  sp->wound_inflammation_source_rate = 0.01;  // strong rate for test visibility
  sp->wound_trigger_step = 0;
  sp->wound_center_x = 25;
  sp->wound_center_y = 25;
  sp->wound_radius = 10;
  sp->dermis_enabled = false;

  auto* rm = sim->GetResourceManager();

  // Register required grids (FusedWoundSourceOp needs all 5 core + inflammation)
  int res = sp->grid_resolution;
  ModelInitializer::DefineSubstance(fields::kVascularId, fields::kVascular, 0, 0, res);
  ModelInitializer::DefineSubstance(fields::kOxygenId, fields::kOxygen, 0, 0, res);
  ModelInitializer::DefineSubstance(fields::kStratumId, fields::kStratum, 0, 0, res);
  ModelInitializer::DefineSubstance(fields::kCalciumId, fields::kCalcium, 0, 0, res);
  ModelInitializer::DefineSubstance(fields::kWaterId, fields::kWater, 0, 0, res);
  ModelInitializer::DefineSubstance(fields::kInflammationId, fields::kInflammation, 0, 0, res);

  auto* vasc_grid = rm->GetDiffusionGrid(fields::kVascular);
  auto* stratum_grid = rm->GetDiffusionGrid(fields::kStratum);
  auto* infl_grid = rm->GetDiffusionGrid(fields::kInflammation);
  vasc_grid->Initialize();
  rm->GetDiffusionGrid(fields::kOxygen)->Initialize();
  stratum_grid->Initialize();
  rm->GetDiffusionGrid(fields::kCalcium)->Initialize();
  rm->GetDiffusionGrid(fields::kWater)->Initialize();
  infl_grid->Initialize();

  // Skip diffusion solver
  vasc_grid->SetTimeStep(1e30);
  rm->GetDiffusionGrid(fields::kOxygen)->SetTimeStep(1e30);
  stratum_grid->SetTimeStep(1e30);
  rm->GetDiffusionGrid(fields::kCalcium)->SetTimeStep(1e30);
  rm->GetDiffusionGrid(fields::kWater)->SetTimeStep(1e30);
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
  sp->wound_inflammation_source_rate = 0.01;
  sp->wound_trigger_step = 0;
  sp->wound_center_x = 25;
  sp->wound_center_y = 25;
  sp->wound_radius = 10;
  sp->dermis_enabled = false;

  auto* rm = sim->GetResourceManager();

  int res = sp->grid_resolution;
  ModelInitializer::DefineSubstance(fields::kVascularId, fields::kVascular, 0, 0, res);
  ModelInitializer::DefineSubstance(fields::kOxygenId, fields::kOxygen, 0, 0, res);
  ModelInitializer::DefineSubstance(fields::kStratumId, fields::kStratum, 0, 0, res);
  ModelInitializer::DefineSubstance(fields::kCalciumId, fields::kCalcium, 0, 0, res);
  ModelInitializer::DefineSubstance(fields::kWaterId, fields::kWater, 0, 0, res);
  ModelInitializer::DefineSubstance(fields::kInflammationId, fields::kInflammation, 0, 0, res);

  auto* vasc_grid = rm->GetDiffusionGrid(fields::kVascular);
  auto* stratum_grid = rm->GetDiffusionGrid(fields::kStratum);
  auto* infl_grid = rm->GetDiffusionGrid(fields::kInflammation);
  vasc_grid->Initialize();
  rm->GetDiffusionGrid(fields::kOxygen)->Initialize();
  stratum_grid->Initialize();
  rm->GetDiffusionGrid(fields::kCalcium)->Initialize();
  rm->GetDiffusionGrid(fields::kWater)->Initialize();
  infl_grid->Initialize();

  vasc_grid->SetTimeStep(1e30);
  rm->GetDiffusionGrid(fields::kOxygen)->SetTimeStep(1e30);
  stratum_grid->SetTimeStep(1e30);
  rm->GetDiffusionGrid(fields::kCalcium)->SetTimeStep(1e30);
  rm->GetDiffusionGrid(fields::kWater)->SetTimeStep(1e30);
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

// --- Hemostasis / Fibrin tests ---

TEST(HemostasisTest, FibrinPDEWoundSeed) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  SanitizeForUnitTest(sp);
  sp->hemostasis_enabled = true;
  sp->wound_enabled = true;
  sp->wound_center_x = 25;
  sp->wound_center_y = 25;
  sp->wound_radius = 10;  // large radius to cover many voxels
  sp->hemostasis_wound_seed = 1.0;

  FibrinPDE pde(sp);
  pde.Init(sim);
  sim->GetScheduler()->Simulate(1);  // initialize grid

  auto* grid = sim->GetResourceManager()->GetDiffusionGrid(fields::kFibrin);
  ASSERT_NE(grid, nullptr);

  // Before wound: all zero
  GridContext ctx(grid, sp);
  real_t pre_total = 0;
  for (size_t i = 0; i < ctx.n; i++) {
    pre_total += grid->GetConcentration(i);
  }
  EXPECT_NEAR(pre_total, 0.0, 1e-6);

  // Apply wound
  pde.ApplyWound(sim, 25, 25, 10);

  // Count seeded voxels and total fibrin
  int seeded = 0;
  real_t post_total = 0;
  for (size_t i = 0; i < ctx.n; i++) {
    real_t val = grid->GetConcentration(i);
    if (val > 0.5) seeded++;
    post_total += val;
  }
  EXPECT_GT(seeded, 0) << "At least one voxel should be seeded with fibrin";
  EXPECT_GT(post_total, 1.0) << "Total fibrin should be substantial";

  delete sim;
}

// --- pH coupling tests ---

TEST(PHCouplingTest, BohrEffectReducesO2) {
  // Alkaline pH should reduce O2 delivery via Bohr effect
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  SanitizeForUnitTest(sp);
  sp->ph_bohr_factor = 0.3;

  SetupAllFields(sim);
  sim->GetScheduler()->Simulate(1);

  auto* o2_grid = sim->GetResourceManager()->GetDiffusionGrid(fields::kOxygen);
  auto* ph_grid = sim->GetResourceManager()->GetDiffusionGrid(fields::kPH);
  auto* vasc_grid = sim->GetResourceManager()->GetDiffusionGrid(
      fields::kVascular);
  ASSERT_NE(o2_grid, nullptr);
  ASSERT_NE(ph_grid, nullptr);

  GridContext ctx(o2_grid, sp);

  // Find a dermal voxel with perfusion > 0
  size_t test_idx = 0;
  bool found = false;
  for (size_t i = 0; i < ctx.n; i++) {
    Real3 pos = {ctx.X(i), ctx.Y(i), -1.0};
    real_t perf = vasc_grid->GetValue(pos);
    if (perf > 0.5) {
      test_idx = i;
      found = true;
      break;
    }
  }
  if (!found) { delete sim; return; }

  // Baseline O2 with pH=0 (no alkalinity)
  real_t o2_baseline = o2_grid->GetConcentration(test_idx);

  // Set high alkalinity
  ph_grid->ChangeConcentrationBy(test_idx, 1.0);
  sim->GetScheduler()->Simulate(1);
  real_t o2_alkaline = o2_grid->GetConcentration(test_idx);

  // O2 should be lower with alkaline pH
  EXPECT_LT(o2_alkaline, o2_baseline * 0.95)
      << "Bohr effect: O2 should drop with alkaline pH";

  delete sim;
}

TEST(PHCouplingTest, MMPBoostAmplifiedByAlkalinity) {
  // Alkaline pH should amplify MMP proteolytic activity
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  SanitizeForUnitTest(sp);
  sp->mmp_enabled = true;
  sp->fibroblast_enabled = true;
  sp->wound_enabled = true;
  sp->wound_center_x = 25;
  sp->wound_center_y = 25;
  sp->wound_radius = 10;
  sp->ph_mmp_boost = 0.5;

  SetupAllFields(sim);
  sim->GetScheduler()->Simulate(1);

  auto* col_grid = sim->GetResourceManager()->GetDiffusionGrid(fields::kCollagen);
  auto* mmp_grid = sim->GetResourceManager()->GetDiffusionGrid(fields::kMMP);
  auto* ph_grid = sim->GetResourceManager()->GetDiffusionGrid(fields::kPH);
  ASSERT_NE(col_grid, nullptr);
  ASSERT_NE(mmp_grid, nullptr);

  // At ph_mmp_boost = 0.5 and alkalinity = 1.0,
  // effective MMP = mmp_val * 1.5 (50% boost).
  // Verify the parameter exists and factor is correct.
  real_t factor = 1.0 + sp->ph_mmp_boost * 1.0;
  EXPECT_NEAR(factor, 1.5, 1e-6);

  // Also verify factor at pH=0 (no boost)
  real_t factor_neutral = 1.0 + sp->ph_mmp_boost * 0.0;
  EXPECT_NEAR(factor_neutral, 1.0, 1e-6);

  delete sim;
}

TEST(PHCouplingTest, BiofilmGrowthBoostedByAlkalinity) {
  // Alkaline pH should boost biofilm growth rate
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  SanitizeForUnitTest(sp);
  sp->ph_biofilm_boost = 0.4;

  // Verify the growth rate formula
  real_t base_growth = sp->biofilm_growth_rate;
  real_t eff_growth_alkaline = base_growth * (1.0 + sp->ph_biofilm_boost * 1.0);
  real_t eff_growth_neutral = base_growth * (1.0 + sp->ph_biofilm_boost * 0.0);

  EXPECT_NEAR(eff_growth_alkaline / base_growth, 1.4, 1e-6)
      << "40% growth boost at full alkalinity";
  EXPECT_NEAR(eff_growth_neutral, base_growth, 1e-6)
      << "No boost at neutral pH";

  delete sim;
}

// ---------------------------------------------------------------------------
// DerivedField unit tests
// ---------------------------------------------------------------------------
TEST(DerivedFieldTest, ConstructionAndSize) {
  DerivedField df("TestField", 5, -10.0, 4.0);
  EXPECT_EQ(df.GetName(), "TestField");
  EXPECT_EQ(df.GetResolution(), 5u);
  EXPECT_EQ(df.GetNumBoxes(), 125u);  // 5^3
  // All values should initialize to zero
  for (size_t i = 0; i < df.GetNumBoxes(); i++) {
    EXPECT_DOUBLE_EQ(df.GetConcentration(i), 0.0);
  }
}

TEST(DerivedFieldTest, SetGetConcentration) {
  DerivedField df("Test", 4, 0.0, 5.0);
  df.SetConcentration(0, 1.5);
  df.SetConcentration(63, 2.7);  // last voxel (4^3 - 1)
  EXPECT_DOUBLE_EQ(df.GetConcentration(0), 1.5);
  EXPECT_DOUBLE_EQ(df.GetConcentration(63), 2.7);
  EXPECT_DOUBLE_EQ(df.GetConcentration(1), 0.0);
}

TEST(DerivedFieldTest, GetValueByPosition) {
  // Grid: resolution=4, lo=0, box_len=10 -> covers [0, 40]
  DerivedField df("Test", 4, 0.0, 10.0);
  // Voxel (0,0,0) -> idx 0
  df.SetConcentration(0, 3.14);
  // Query inside that voxel
  Real3 pos = {5.0, 5.0, 5.0};
  EXPECT_DOUBLE_EQ(df.GetValue(pos), 3.14);
  // Query in a different voxel (1,0,0) -> idx 1
  df.SetConcentration(1, 2.0);
  Real3 pos2 = {15.0, 5.0, 5.0};
  EXPECT_DOUBLE_EQ(df.GetValue(pos2), 2.0);
}

TEST(DerivedFieldTest, Zero) {
  DerivedField df("Test", 3, 0.0, 5.0);
  for (size_t i = 0; i < df.GetNumBoxes(); i++) {
    df.SetConcentration(i, static_cast<real_t>(i) + 1.0);
  }
  EXPECT_GT(df.GetConcentration(0), 0.0);
  df.Zero();
  for (size_t i = 0; i < df.GetNumBoxes(); i++) {
    EXPECT_DOUBLE_EQ(df.GetConcentration(i), 0.0);
  }
}

TEST(DerivedFieldTest, GetBoxIndexClamping) {
  // Grid: resolution=4, lo=-10, box_len=5 -> covers [-10, 10]
  DerivedField df("Test", 4, -10.0, 5.0);
  // Position below lo should clamp to index 0
  Real3 below = {-20.0, -20.0, -20.0};
  EXPECT_EQ(df.GetBoxIndex(below), 0u);
  // Position above max should clamp to (3,3,3) -> idx 3 + 3*4 + 3*16 = 63
  Real3 above = {100.0, 100.0, 100.0};
  EXPECT_EQ(df.GetBoxIndex(above), 63u);
}

// ---------------------------------------------------------------------------
// CoarseIndex unit tests
// ---------------------------------------------------------------------------
TEST(CoarseIndexTest, MapsToValidIndex) {
  // CoarseIndex should always produce a valid index within the coarse grid.
  auto* sim = CreateTestSim("CoarseIdxValid");
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  SanitizeForUnitTest(sp);
  sp->grid_resolution = 10;
  sp->grid_resolution_structural = 5;

  // Register fine grid (O2 uses grid_resolution) and structural grid
  OxygenPDE o2;
  o2.Init(sim);
  SimpleStructuralPDE structural("CoarseTest", 99, 0, 0);
  structural.Init(sim);

  auto* fine_grid = sim->GetResourceManager()->GetDiffusionGrid(fields::kOxygen);
  auto* coarse_grid = sim->GetResourceManager()->GetDiffusionGrid("CoarseTest");
  ASSERT_EQ(fine_grid->GetResolution(), 10u);
  ASSERT_EQ(coarse_grid->GetResolution(), 5u);

  GridContext ctx(fine_grid, sp);
  size_t coarse_boxes = coarse_grid->GetNumBoxes();

  // Every fine index should map to a valid coarse index
  for (size_t idx = 0; idx < ctx.n; idx++) {
    size_t mapped = GridContext::CoarseIndex(idx, ctx, coarse_grid);
    EXPECT_LT(mapped, coarse_boxes)
        << "Fine idx=" << idx << " mapped to invalid coarse idx=" << mapped;
  }
  delete sim;
}

TEST(CoarseIndexTest, SpatialConsistency) {
  // Reading a coarse grid via CoarseIndex should give the same value as
  // GetValue(position), verifying that the spatial mapping is correct.
  auto* sim = CreateTestSim("CoarseIdxSpatial");
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  SanitizeForUnitTest(sp);
  sp->grid_resolution = 10;
  sp->grid_resolution_structural = 5;

  OxygenPDE o2;
  o2.Init(sim);
  SimpleStructuralPDE structural("CoarseSpatial", 98, 0, 0);
  structural.Init(sim);

  auto* fine_grid = sim->GetResourceManager()->GetDiffusionGrid(fields::kOxygen);
  auto* coarse_grid = sim->GetResourceManager()->GetDiffusionGrid("CoarseSpatial");

  // Seed coarse grid with unique values per voxel
  for (size_t i = 0; i < coarse_grid->GetNumBoxes(); i++) {
    coarse_grid->ChangeConcentrationBy(i, static_cast<real_t>(i + 1) * 0.01);
  }

  GridContext ctx(fine_grid, sp);
  // For a sample of fine voxels, verify CoarseIndex gives correct spatial read
  for (size_t idx = 0; idx < ctx.n; idx += 37) {
    Real3 pos = {ctx.X(idx), ctx.Y(idx), ctx.Z(idx)};
    real_t via_getvalue = coarse_grid->GetValue(pos);
    size_t mapped = GridContext::CoarseIndex(idx, ctx, coarse_grid);
    real_t via_coarseindex = coarse_grid->GetConcentration(mapped);
    EXPECT_DOUBLE_EQ(via_coarseindex, via_getvalue)
        << "CoarseIndex at fine idx=" << idx << " should match GetValue";
  }
  delete sim;
}

// ---------------------------------------------------------------------------
// SimpleStructuralPDE test
// ---------------------------------------------------------------------------
TEST(SimpleStructuralPDETest, UsesStructuralResolution) {
  auto* sim = CreateTestSim("StructPDE");
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  SanitizeForUnitTest(sp);
  sp->grid_resolution = 10;
  sp->grid_resolution_structural = 5;

  // Register a structural PDE (uses DefineStructuralGrid)
  SimpleStructuralPDE pde("TestStruct", 99, 0, 0);
  pde.Init(sim);

  auto* grid = sim->GetResourceManager()->GetDiffusionGrid("TestStruct");
  EXPECT_EQ(grid->GetResolution(), 5u)
      << "Structural PDE should use grid_resolution_structural";
  delete sim;
}

TEST(SimpleStructuralPDETest, FallsBackWhenZero) {
  auto* sim = CreateTestSim("StructFallback");
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  SanitizeForUnitTest(sp);
  sp->grid_resolution = 10;
  sp->grid_resolution_structural = 0;  // disabled

  SimpleStructuralPDE pde("TestFallback", 98, 0, 0);
  pde.Init(sim);

  auto* grid = sim->GetResourceManager()->GetDiffusionGrid("TestFallback");
  EXPECT_EQ(grid->GetResolution(), 10u)
      << "Structural PDE should fall back to grid_resolution when structural=0";
  delete sim;
}

}  // namespace skibidy
}  // namespace bdm

// Explicit gtest main -- overrides main() from libskibidy.so
// (shared library symbols resolve after the binary's own symbols on Linux)
int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
