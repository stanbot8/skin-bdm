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
//
// Core infrastructure tests extracted from test-suite.cc.
//

#include "test_helpers.h"
#include "core/chunked_grid.h"
#include "core/fused_source.h"
#include "core/derived_field.h"
#include "core/derived_fields_op.h"
#include "core/metrics.h"
#include "wound/wound_event.h"
#include "immune/immune_response.h"
#include "fibroblast/fibroblast.h"
#include "fibroblast/fibroblast_behavior.h"
#include "fibroblast/fibroblast_recruitment.h"
#include "biofilm/biofilm_op.h"
#include "angiogenesis/vegf_op.h"
#include "tumor/tumor_cell.h"
#include "tumor/tumor_behavior.h"
#include "tumor/tumor_initiation.h"
#include "tumor/tumor_pde.h"

namespace bdm {
namespace skibidy {

// ---------------------------------------------------------------------------
// CompositeFieldTest
// ---------------------------------------------------------------------------
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

// ---------------------------------------------------------------------------
// IntegrationTest: WoundTriggersAndCellsExist
// ---------------------------------------------------------------------------
TEST(IntegrationTest, WoundTriggersAndCellsExist) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  sp->wound.trigger_step = 5;

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

// ---------------------------------------------------------------------------
// GridContextTest
// ---------------------------------------------------------------------------
TEST(GridContextTest, InWoundCylinderDetection) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  sp->wound.center_x = 20;
  sp->wound.center_y = 20;
  sp->wound.radius = 5;

  CalciumPDE ca;
  ca.Init(sim);
  sim->GetScheduler()->Simulate(1);

  auto* grid = sim->GetResourceManager()->GetDiffusionGrid(fields::kCalciumId);
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

  auto* grid = sim->GetResourceManager()->GetDiffusionGrid(fields::kCalciumId);
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

// ---------------------------------------------------------------------------
// WoundMaskTest
// ---------------------------------------------------------------------------
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
    sp->wound.center_x = 20;
    sp->wound.center_y = 20;
    sp->wound.radius = 8;
  };
  auto* sim = new Simulation(TEST_NAME, set_param);
  auto* sp = sim->GetParam()->Get<SimParam>();

  VascularPDE vasc;
  vasc.Init(sim);
  sim->GetScheduler()->Simulate(1);
  auto* grid = sim->GetResourceManager()->GetDiffusionGrid(fields::kVascularId);
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

// ---------------------------------------------------------------------------
// MetricsExporterTest
// ---------------------------------------------------------------------------
TEST(MetricsExporterTest, WritesHeaderOnFirstCall) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  sp->metrics_interval = 1;
  sp->wound.enabled = false;
  sp->tumor.enabled = false;

  SetupAllFields(sim);

  MetricsExporter metrics;
  metrics();

  // Read the output file
  std::string path = sim->GetOutputDir() + "/metrics.csv";
  std::ifstream f(path);
  ASSERT_TRUE(f.is_open());

  std::string header;
  std::getline(f, header);

  // Header should contain 45 column names
  EXPECT_NE(header.find("step"), std::string::npos);
  EXPECT_NE(header.find("time_h"), std::string::npos);
  EXPECT_NE(header.find("wound_closure_pct"), std::string::npos);
  EXPECT_NE(header.find("n_tumor_cells"), std::string::npos);
  EXPECT_NE(header.find("tumor_field_cells"), std::string::npos);
  EXPECT_NE(header.find("mean_ecm_quality"), std::string::npos);
  EXPECT_NE(header.find("mean_tissue_viability"), std::string::npos);
  EXPECT_NE(header.find("n_activated_fibroblasts"), std::string::npos);
  EXPECT_NE(header.find("mean_temperature_wound"), std::string::npos);
  EXPECT_NE(header.find("mean_glucose_wound"), std::string::npos);
  EXPECT_NE(header.find("mean_lactate_wound"), std::string::npos);
  EXPECT_NE(header.find("mean_no_wound"), std::string::npos);
  EXPECT_NE(header.find("mean_ros_wound"), std::string::npos);
  EXPECT_NE(header.find("mean_stiffness_wound"), std::string::npos);
  EXPECT_NE(header.find("mean_lymphatic_wound"), std::string::npos);
  EXPECT_NE(header.find("mean_edema_wound"), std::string::npos);
  EXPECT_NE(header.find("mean_voltage_wound"), std::string::npos);
  EXPECT_NE(header.find("mean_tcell_wound"), std::string::npos);
  EXPECT_NE(header.find("mean_bone_wound"), std::string::npos);
  EXPECT_NE(header.find("mean_scab_wound"), std::string::npos);

  // Count commas to verify column count (56 commas for 57 columns)
  int commas = 0;
  for (char c : header) {
    if (c == ',') commas++;
  }
  EXPECT_EQ(commas, 56);

  metrics.Close();
  delete sim;
}

TEST(MetricsExporterTest, RowCountMatchesInterval) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  sp->metrics_interval = 3;
  sp->wound.enabled = false;
  sp->tumor.enabled = false;

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

// ---------------------------------------------------------------------------
// IntegrationTest: WoundClosureProgressesOverTime
// ---------------------------------------------------------------------------
TEST(IntegrationTest, WoundClosureProgressesOverTime) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  sp->wound.trigger_step = 0;
  sp->wound.enabled = true;
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
  if (sp->inflammation.split_inflammation_enabled) {
    fields.Add(std::make_unique<ProInflammatoryPDE>(sp));
    fields.Add(std::make_unique<AntiInflammatoryPDE>(sp));
  } else {
    fields.Add(std::make_unique<InflammationPDE>(sp));
  }
  fields.Add(std::make_unique<ImmunePressurePDE>(sp));
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

// ---------------------------------------------------------------------------
// IntegrationTest: TumorGrowsOverTime
// ---------------------------------------------------------------------------
TEST(IntegrationTest, TumorGrowsOverTime) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  sp->tumor.enabled = true;
  sp->tumor.seed_time = 0;
  sp->tumor.seed_count = 5;
  sp->tumor.max_cells = 500;
  sp->tumor.cycle_factor = 0.5;   // fast cycle so growth is reliable in 500 steps
  sp->tumor.growth_rate = 20;     // fast growth to reach division_diameter reliably
  sp->wound.enabled = false;

  SetupAllFields(sim);

  // Initialize O2 so tumor cells can proliferate (G1->S requires o2 > 0)
  ModelInitializer::InitializeSubstance(fields::kOxygenId,
      [](real_t, real_t, real_t) { return 1.0; });

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
  EXPECT_GT(tumor_count - 1, static_cast<int>(sp->tumor.seed_count))
      << "Tumor should grow beyond initial seed count after 500 steps";

  delete sim;
}

// ---------------------------------------------------------------------------
// ChunkedGridTest
// ---------------------------------------------------------------------------
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
// VolumeManagerTest
// ---------------------------------------------------------------------------
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

// ---------------------------------------------------------------------------
// ResolutionTest (with VEGFResResult struct and helper)
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
    sp->wound.trigger_step = 0;
    sp->angiogenesis.enabled = true;
    sp->diabetic.mode = false;
    sp->biofilm.enabled = false;
    sp->ra.enabled = false;
    sp->mmp.enabled = false;
    sp->fibronectin.enabled = false;
    sp->fibroblast.enabled = false;
    sp->elastin.enabled = false;
    sp->hyaluronan.enabled = false;
    sp->angiogenesis.vegf_production_rate = 0.01;
    sp->angiogenesis.vegf_hypoxia_threshold = 0.5;
    sp->oxygen_basal_conc = 1.0;
    sp->perfusion.basal = 1.0;
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

  // Force environment update so grids get valid bounds before Initialize()
  sim->GetEnvironment()->Update();

  auto* vasc_grid = rm->GetDiffusionGrid(fields::kVascularId);
  auto* vegf_grid = rm->GetDiffusionGrid(fields::kVEGFId);
  vasc_grid->Initialize();
  rm->GetDiffusionGrid(fields::kOxygenId)->Initialize();
  rm->GetDiffusionGrid(fields::kStratumId)->Initialize();
  rm->GetDiffusionGrid(fields::kCalciumId)->Initialize();
  rm->GetDiffusionGrid(fields::kWaterId)->Initialize();
  vegf_grid->Initialize();

  // Skip diffusion solver
  vasc_grid->SetTimeStep(1e30);
  rm->GetDiffusionGrid(fields::kOxygenId)->SetTimeStep(1e30);
  rm->GetDiffusionGrid(fields::kStratumId)->SetTimeStep(1e30);
  rm->GetDiffusionGrid(fields::kCalciumId)->SetTimeStep(1e30);
  rm->GetDiffusionGrid(fields::kWaterId)->SetTimeStep(1e30);
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
// ScaledGridTest
// ---------------------------------------------------------------------------
TEST(ScaledGridTest, FactorAtCalibrationResolution) {
  auto* sim = CreateTestSim("sg_factor_cal");
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  sp->ref_box_length = 5.0;
  sp->grid_resolution = 10;
  SanitizeForUnitTest(sp);

  InflammationPDE infl(sp); infl.Init(sim);
  sim->GetScheduler()->Simulate(1);  // initializes grid internals
  auto* grid = sim->GetResourceManager()->GetDiffusionGrid(fields::kInflammationId);

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
  auto* grid = sim->GetResourceManager()->GetDiffusionGrid(fields::kInflammationId);

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
      fields::kInflammationId);

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
// PerfTimerTest
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

// ---------------------------------------------------------------------------
// DerivedFieldTest
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

  auto* fine_grid = sim->GetResourceManager()->GetDiffusionGrid(fields::kOxygenId);
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

  auto* fine_grid = sim->GetResourceManager()->GetDiffusionGrid(fields::kOxygenId);
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
// SimpleStructuralPDETest
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
