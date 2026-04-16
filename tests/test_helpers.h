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
// Shared test helpers for skibidy unit tests.
// Each test file includes this header for CreateTestSim(), SanitizeForUnitTest(),
// SetupInflammation(), SetupAllFields(), and GetInflammationGrid().
//
#ifndef SKIBIDY_TEST_HELPERS_H_
#define SKIBIDY_TEST_HELPERS_H_

#include <gtest/gtest.h>
#include "biodynamo.h"
#include "infra/toml_compat.h"
#include "infra/sim_param.h"
#include "skibidy.h"
#include "core/pde.h"
#include "core/composite_field.h"
#include "core/field_names.h"
#include "tissue/calcium.h"
#include "tissue/kgf.h"
#include "tissue/oxygen.h"
#include "tissue/stratum.h"
#include "tissue/water.h"
#include "inflammation/inflammation_pde.h"
#include "scar/scar_pde.h"
#include "fibroblast/tgfbeta_pde.h"
#include "fibroblast/collagen_pde.h"
#include "fibronectin/fibronectin_pde.h"
#include "perfusion/vascular.h"
#include "ph/ph_pde.h"
#include "ros/ros_pde.h"
#include "mechanotransduction/stiffness_pde.h"
#include "lymphatic/lymphatic_pde.h"
#include "lymphatic/edema_pde.h"
#include "bioelectric/voltage_pde.h"
#include "mmp/mmp_pde.h"

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
    sp->tissue_min = 0;        // match test domain (with headroom for boundary tests)
    sp->tissue_max = 45;       // leave 5 units inside domain for boundary tests
    sp->wound.center_x = 25;   // center of test domain
    sp->wound.center_y = 25;   // center of test domain
    sp->grid_resolution = 10;  // coarse grid for fast tests
  };
  return new Simulation(name, set_param);
}

// ---------------------------------------------------------------------------
// Helper: register inflammation fields matching config (split vs single)
// ---------------------------------------------------------------------------
static void SetupInflammation(Simulation* sim) {
  auto* sp = sim->GetParam()->Get<SimParam>();
  if (sp->inflammation.split_inflammation_enabled) {
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
  sp->diabetic.mode = false;
  sp->efferocytosis_enabled = false;
  sp->biofilm.enabled = false;
  sp->mmp.enabled = false;
  sp->angiogenesis.enabled = false;
  sp->fibronectin.enabled = false;
  sp->fibroblast.enabled = false;
  sp->elastin.enabled = false;
  sp->hyaluronan.enabled = false;
  sp->hemostasis.enabled = false;
  // New physics modules: disable to avoid costly "grid not found" error logging
  sp->temperature.enabled = false;
  sp->glucose_mod.enabled = false;
  sp->lactate.enabled = false;
  sp->nitric_oxide.enabled = false;
  sp->dermis.enabled = false;
  sp->mechanotransduction.enabled = false;
  sp->lymphatic.enabled = false;
  sp->bioelectric.enabled = false;
  // RA module
  sp->ra.enabled = false;
  // Reset multi-resolution so unit tests use uniform grid resolution
  sp->grid_resolution_structural = 0;
}

// Same but keeps diabetic_mode true (for DiabeticTest tests)
static void SanitizeForDiabeticTest(SimParam* sp) {
  SanitizeForUnitTest(sp);
  sp->diabetic.mode = true;
  sp->inflammation.split_inflammation_enabled = false;
}

// Helper: get the readable inflammation grid (pro-infl in split mode)
static DiffusionGrid* GetInflammationGrid(Simulation* sim) {
  auto* sp = sim->GetParam()->Get<SimParam>();
  auto* rm = sim->GetResourceManager();
  if (sp->inflammation.split_inflammation_enabled) {
    return rm->GetDiffusionGrid(fields::kProInflammatoryId);
  }
  return rm->GetDiffusionGrid(fields::kInflammationId);
}

// ---------------------------------------------------------------------------
// Helper: initialize all PDE channels for tests that read field values
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
  if (sp->fibronectin.enabled) { FibronectinPDE fn(sp); fn.Init(sim); }
  if (sp->mmp.enabled) { MMPPDE mmp(sp); mmp.Init(sim); }
  if (sp->fibroblast.enabled) {
    TGFBetaPDE tgfb(sp); tgfb.Init(sim);
    CollagenPDE col(sp); col.Init(sim);
  }
  { PHPDE ph(sp); ph.Init(sim); }
  if (sp->ros.enabled) { ROSPDE ros(sp); ros.Init(sim); }
  if (sp->mechanotransduction.enabled) { StiffnessPDE stiff(sp); stiff.Init(sim); }
  if (sp->lymphatic.enabled) {
    LymphaticPDE lymph(sp); lymph.Init(sim);
    EdemaPDE edema(sp); edema.Init(sim);
  }
  if (sp->bioelectric.enabled) { VoltagePDE volt(sp); volt.Init(sim); }
}

}  // namespace skibidy
}  // namespace bdm

#endif  // SKIBIDY_TEST_HELPERS_H_
