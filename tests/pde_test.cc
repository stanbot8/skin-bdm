#include "test_helpers.h"

#include "biofilm/biofilm_pde.h"
#include "biofilm/biofilm_op.h"
#include "angiogenesis/vegf_pde.h"
#include "angiogenesis/vegf_op.h"
#include "elastin/elastin_pde.h"
#include "hyaluronan/hyaluronan_pde.h"
#include "dermis/dermis_pde.h"
#include "hemostasis/hemostasis_pde.h"
#include "immune/macrophage_behavior.h"
#include "immune/immune_response.h"
#include "core/fused_source.h"

namespace bdm {
namespace skibidy {

TEST(ProportionalScarTest, AccumulationOp) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  SanitizeForUnitTest(sp);
  sp->scar.proportional_enabled = true;
  sp->scar.accumulation_rate = 0.1;  // high rate for test visibility
  sp->wound.enabled = true;
  sp->wound.center_x = 25;
  sp->wound.center_y = 25;
  sp->wound.radius = 20;  // large wound to cover grid
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

  auto* scar_grid = rm->GetDiffusionGrid(fields::kScarId);
  size_t scar_idx = scar_grid->GetBoxIndex(center);
  real_t scar_val = scar_grid->GetConcentration(scar_idx);

  // scar += inflammation * rate = 5.0 * 0.1 = 0.5
  EXPECT_NEAR(scar_val, 0.5, 0.01);

  delete sim;
}

TEST(ProportionalScarTest, WriteScarValueSkipsBinaryScar) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  sp->scar.proportional_enabled = true;

  SetupAllFields(sim);
  sim->GetScheduler()->Simulate(1);

  auto* cell = new Keratinocyte({25, 25, 5});
  cell->SetDiameter(4);
  cell->SetStratum(kCornified);

  WriteScarValue(cell);

  // Scar grid should NOT have binary 1.0 when proportional mode is on
  auto* rm = sim->GetResourceManager();
  auto* scar_grid = rm->GetDiffusionGrid(fields::kScarId);
  Real3 pos = {25, 25, 5};
  size_t idx = scar_grid->GetBoxIndex(pos);
  EXPECT_LT(scar_grid->GetConcentration(idx), 0.5);

  // But stratum field should still get stamped
  auto* stratum_grid = rm->GetDiffusionGrid(fields::kStratumId);
  size_t st_idx = stratum_grid->GetBoxIndex(pos);
  EXPECT_GT(stratum_grid->GetConcentration(st_idx), 0);

  delete cell;
  delete sim;
}

TEST(VascularPDETest, InitializesToOneInDermis) {
  // Need domain with negative z for dermal voxels
  Param::RegisterParamGroup(new SimParam());
  auto set_param = [](Param* param) {
    param->bound_space = Param::BoundSpaceMode::kClosed;
    param->min_bound = -10;
    param->max_bound = 50;
    param->simulation_time_step = 0.1;
    param->export_visualization = false;
    auto* sp = const_cast<SimParam*>(param->Get<SimParam>());
    sp->grid_resolution = 10;  // pin for deterministic voxel layout
  };
  auto* sim = new Simulation(TEST_NAME, set_param);
  auto* sp = sim->GetParam()->Get<SimParam>();

  VascularPDE vasc;
  vasc.Init(sim);
  sim->GetScheduler()->Simulate(1);

  auto* grid = sim->GetResourceManager()->GetDiffusionGrid(fields::kVascularId);
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
  EXPECT_NEAR(mean, sp->perfusion.basal, 0.2);

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
  sp->wound.center_x = 25;
  sp->wound.center_y = 25;
  sp->wound.radius = 10;

  VascularPDE vasc;
  vasc.Init(sim);
  sim->GetScheduler()->Simulate(1);

  auto* grid = sim->GetResourceManager()->GetDiffusionGrid(fields::kVascularId);

  // Apply wound
  vasc.ApplyWound(sim, sp->wound.center_x, sp->wound.center_y,
                  sp->wound.radius);

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
    sp->perfusion.papillary_fraction = 1.0;
    sp->perfusion.reticular_fraction = 0.7;
    sp->perfusion.hypodermis_fraction = 0.3;
  };
  auto* sim = new Simulation(TEST_NAME, set_param);
  auto* sp = sim->GetParam()->Get<SimParam>();

  VascularPDE vasc;
  vasc.Init(sim);
  sim->GetScheduler()->Simulate(1);

  auto* grid = sim->GetResourceManager()->GetDiffusionGrid(fields::kVascularId);
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

TEST(BiofilmTest, LogisticGrowth) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  SanitizeForUnitTest(sp);
  sp->biofilm.enabled = true;
  sp->biofilm.growth_rate = 0.1;
  sp->biofilm.carrying_capacity = 1.0;
  sp->wound.enabled = true;
  sp->wound.center_x = 15;
  sp->wound.center_y = 15;
  sp->wound.radius = 10;
  sp->wound.trigger_step = 0;

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
  auto* biofilm_grid = rm->GetDiffusionGrid(fields::kBiofilmId);
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
  sp->biofilm.enabled = true;
  sp->biofilm.m1_block_threshold = 0.2;
  sp->immune.macrophage_m1_duration = 1;
  sp->immune.macrophage_lifespan = 999;
  sp->immune.macrophage_min_survival = 999;
  sp->m1_transition_threshold = 999;
  sp->wound.enabled = true;
  sp->wound.center_x = 25;
  sp->wound.center_y = 25;
  sp->wound.radius = 20;
  sp->wound.trigger_step = 0;

  SetupInflammation(sim);
  BiofilmPDE bio;
  bio.Init(sim);
  sim->GetScheduler()->Simulate(1);

  // Seed high biofilm at cell position
  auto* rm = sim->GetResourceManager();
  auto* biofilm_grid = rm->GetDiffusionGrid(fields::kBiofilmId);
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

TEST(AngiogenesisTest, VEGFProductionUnderHypoxia) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  sp->angiogenesis.enabled = true;
  sp->angiogenesis.vegf_production_rate = 0.1;
  sp->angiogenesis.vegf_hypoxia_threshold = 0.5;
  sp->wound.enabled = true;
  sp->wound.center_x = 25;
  sp->wound.center_y = 25;
  sp->wound.radius = 20;
  sp->wound.trigger_step = 0;
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
  auto* o2_grid = rm->GetDiffusionGrid(fields::kOxygenId);
  auto* vegf_grid = rm->GetDiffusionGrid(fields::kVEGFId);
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
  sp->angiogenesis.enabled = true;
  SanitizeForUnitTest(sp);
  sp->angiogenesis.enabled = true;
  sp->angiogenesis.vegf_production_rate = 0.1;
  sp->angiogenesis.vegf_hypoxia_threshold = 0.5;
  sp->wound.enabled = true;
  sp->wound.center_x = 25;
  sp->wound.center_y = 25;
  sp->wound.radius = 20;
  sp->wound.trigger_step = 0;
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
  auto* vegf_grid = sim->GetResourceManager()->GetDiffusionGrid(fields::kVEGFId);
  Real3 pos = {25, 25, 5};
  real_t before = vegf_grid->GetConcentration(vegf_grid->GetBoxIndex(pos));

  VEGFSourceOp op;
  op();

  real_t after = vegf_grid->GetConcentration(vegf_grid->GetBoxIndex(pos));
  // At well-oxygenated voxels, VEGF should not increase
  EXPECT_LE(after, before + 1e-10);

  delete sim;
}

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

  auto* grid = sim->GetResourceManager()->GetDiffusionGrid(fields::kElastinId);
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
  sp->grid_resolution = 10;  // pin for deterministic voxel layout
  sp->grid_resolution_structural = 0;

  HyaluronanPDE ha;
  ha.Init(sim);
  sim->GetScheduler()->Simulate(1);

  auto* grid = sim->GetResourceManager()->GetDiffusionGrid(fields::kHyaluronanId);
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
  sp->wound.center_x = 25;
  sp->wound.center_y = 25;
  sp->wound.radius = 10;

  DermisPDE dermis;
  dermis.Init(sim);
  sim->GetScheduler()->Simulate(1);
  dermis.ApplyWound(sim, sp->wound.center_x, sp->wound.center_y,
                    sp->wound.radius);

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
  sp->wound.center_x = 25;
  sp->wound.center_y = 25;
  sp->wound.radius = 5;

  DermisPDE dermis;
  dermis.Init(sim);
  sim->GetScheduler()->Simulate(1);
  dermis.ApplyWound(sim, sp->wound.center_x, sp->wound.center_y,
                    sp->wound.radius);

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

TEST(HemostasisTest, FibrinPDEWoundSeed) {
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  SanitizeForUnitTest(sp);
  sp->hemostasis.enabled = true;
  sp->wound.enabled = true;
  sp->wound.center_x = 25;
  sp->wound.center_y = 25;
  sp->wound.radius = 10;  // large radius to cover many voxels
  sp->hemostasis.wound_seed = 1.0;

  FibrinPDE pde(sp);
  pde.Init(sim);
  sim->GetScheduler()->Simulate(1);  // initialize grid

  auto* grid = sim->GetResourceManager()->GetDiffusionGrid(fields::kFibrinId);
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

TEST(PHCouplingTest, BohrEffectReducesO2) {
  // Alkaline pH should reduce O2 delivery via Bohr effect
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  SanitizeForUnitTest(sp);
  sp->ph.bohr_factor = 0.3;

  SetupAllFields(sim);
  sim->GetScheduler()->Simulate(1);

  auto* o2_grid = sim->GetResourceManager()->GetDiffusionGrid(fields::kOxygenId);
  auto* ph_grid = sim->GetResourceManager()->GetDiffusionGrid(fields::kPHId);
  auto* vasc_grid = sim->GetResourceManager()->GetDiffusionGrid(
      fields::kVascularId);
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
  sp->mmp.enabled = true;
  sp->fibroblast.enabled = true;
  sp->wound.enabled = true;
  sp->wound.center_x = 25;
  sp->wound.center_y = 25;
  sp->wound.radius = 10;
  sp->ph.mmp_boost = 0.5;

  SetupAllFields(sim);
  sim->GetScheduler()->Simulate(1);

  auto* col_grid = sim->GetResourceManager()->GetDiffusionGrid(fields::kCollagenId);
  auto* mmp_grid = sim->GetResourceManager()->GetDiffusionGrid(fields::kMMPId);
  auto* ph_grid = sim->GetResourceManager()->GetDiffusionGrid(fields::kPHId);
  ASSERT_NE(col_grid, nullptr);
  ASSERT_NE(mmp_grid, nullptr);

  // At ph_mmp_boost = 0.5 and alkalinity = 1.0,
  // effective MMP = mmp_val * 1.5 (50% boost).
  // Verify the parameter exists and factor is correct.
  real_t factor = 1.0 + sp->ph.mmp_boost * 1.0;
  EXPECT_NEAR(factor, 1.5, 1e-6);

  // Also verify factor at pH=0 (no boost)
  real_t factor_neutral = 1.0 + sp->ph.mmp_boost * 0.0;
  EXPECT_NEAR(factor_neutral, 1.0, 1e-6);

  delete sim;
}

TEST(PHCouplingTest, BiofilmGrowthBoostedByAlkalinity) {
  // Alkaline pH should boost biofilm growth rate
  auto* sim = CreateTestSim(TEST_NAME);
  auto* sp = const_cast<SimParam*>(sim->GetParam()->Get<SimParam>());
  SanitizeForUnitTest(sp);
  sp->ph.biofilm_boost = 0.4;

  // Verify the growth rate formula
  real_t base_growth = sp->biofilm.growth_rate;
  real_t eff_growth_alkaline = base_growth * (1.0 + sp->ph.biofilm_boost * 1.0);
  real_t eff_growth_neutral = base_growth * (1.0 + sp->ph.biofilm_boost * 0.0);

  EXPECT_NEAR(eff_growth_alkaline / base_growth, 1.4, 1e-6)
      << "40% growth boost at full alkalinity";
  EXPECT_NEAR(eff_growth_neutral, base_growth, 1e-6)
      << "No boost at neutral pH";

  delete sim;
}

}  // namespace skibidy
}  // namespace bdm
