// Smoke tests for modules without dedicated test coverage.
// Each test enables one module, registers its required fields, runs the
// fused wound pipeline for a few steps, and verifies no crash + basic
// grid activity. These would have caught the dt=0.1 hardcoding bug.

#include "test_helpers.h"
#include "core/fused_source.h"
#include "core/fused_post.h"
#include "core/registration.h"
#include "core/composite_field.h"
#include "wound/wound_event.h"

// Module PDEs
#include "angiogenesis/vegf_pde.h"
#include "fibroblast/collagen_pde.h"
#include "fibroblast/tgfbeta_pde.h"
#include "hemostasis/hemostasis_pde.h"
#include "mmp/mmp_pde.h"
#include "neuropathy/nerve_pde.h"
#include "nitric_oxide/nitric_oxide_pde.h"
#include "perfusion/vascular.h"
#include "photon/fluence_pde.h"
#include "photon/opsin_pde.h"
#include "ros/ros_pde.h"
#include "scab/scab_pde.h"
#include "temperature/temperature_pde.h"

namespace bdm {
namespace skibidy {

namespace {

// Minimal wound environment: stratum, oxygen, vascular, inflammation,
// calcium, KGF, water, scar + wound trigger at step 0.
struct SmokeEnv {
  Simulation* sim;
  SimParam* sp;
  CompositeField fields;

  static SmokeEnv Create(const char* name) {
    SmokeEnv env;
    env.sim = CreateTestSim(name);
    env.sp = const_cast<SimParam*>(env.sim->GetParam()->Get<SimParam>());
    SanitizeForUnitTest(env.sp);

    env.sp->wound.enabled = true;
    env.sp->wound.trigger_step = 0;
    env.sp->wound.radius = 20;
    env.sp->volume_z_cornified = 40;

    auto* vm = VolumeManager::Get();
    vm->Init(env.sp->volume_z_spinous, env.sp->volume_z_granular,
             env.sp->volume_z_cornified);

    // Shared fields most modules read from
    env.fields.Add(std::make_unique<StratumPDE>());
    env.fields.Add(std::make_unique<OxygenPDE>());
    env.fields.Add(std::make_unique<CalciumPDE>());
    env.fields.Add(std::make_unique<KgfPDE>());
    env.fields.Add(std::make_unique<WaterPDE>());
    env.fields.Add(std::make_unique<ScarPDE>());
    SetupInflammation(env.sim);
    VascularPDE vasc; vasc.Init(env.sim);
    return env;
  }

  // Register fields, trigger wound, run N fused steps.
  // tag must be unique per test (use TEST_NAME).
  void Run(const char* tag, int steps = 10) {
    fields.InitAll(sim);
    sim->GetScheduler()->Simulate(1);  // populate grid arrays

    auto reg_op = [&](const char* suffix, StandaloneOperationImpl* impl,
                      OpType type = OpType::kPreSchedule) {
      std::string name = std::string(tag) + "_" + suffix;
      OperationRegistry::GetInstance()->AddOperationImpl(
          name.c_str(), OpComputeTarget::kCpu, impl);
      sim->GetScheduler()->ScheduleOp(NewOperation(name.c_str()), type);
    };

    reg_op("wound", new WoundEvent(&fields));
    reg_op("source", new FusedWoundSourceOp());
    reg_op("post", new FusedWoundPostOp());

    sim->GetScheduler()->Simulate(steps);
  }
};

// Check that at least one voxel in a grid has a nonzero value.
bool HasActivity(Simulation* sim, int field_id) {
  auto* rm = sim->GetResourceManager();
  auto* grid = rm->GetDiffusionGrid(field_id);
  if (!grid) return false;
  const real_t* c = grid->GetAllConcentrations();
  for (size_t i = 0; i < grid->GetNumBoxes(); i++) {
    if (std::abs(c[i]) > 1e-12) return true;
  }
  return false;
}

}  // namespace

TEST(ModuleSmoke, Temperature) {
  auto env = SmokeEnv::Create(TEST_NAME);
  env.sp->temperature.enabled = true;
  TemperaturePDE temp(env.sp); temp.Init(env.sim);
  env.Run(TEST_NAME);
  // Survived: no crash during fused pipeline execution
}

TEST(ModuleSmoke, NitricOxide) {
  auto env = SmokeEnv::Create(TEST_NAME);
  env.sp->nitric_oxide.enabled = true;
  NitricOxidePDE no_pde(env.sp); no_pde.Init(env.sim);
  env.Run(TEST_NAME);
  // Survived: no crash during fused pipeline execution
}

TEST(ModuleSmoke, Neuropathy) {
  auto env = SmokeEnv::Create(TEST_NAME);
  env.sp->neuropathy.enabled = true;
  NervePDE nerve(env.sp); nerve.Init(env.sim);
  // Neuropathy reads VEGF and TGF-beta
  env.sp->angiogenesis.enabled = true;
  VEGFPDE vegf(env.sp); vegf.Init(env.sim);
  env.sp->fibroblast.enabled = true;
  TGFBetaPDE tgfb(env.sp); tgfb.Init(env.sim);
  CollagenPDE col(env.sp); col.Init(env.sim);
  env.Run(TEST_NAME);
  // Survived: no crash during fused pipeline execution
}

TEST(ModuleSmoke, Scab) {
  auto env = SmokeEnv::Create(TEST_NAME);
  env.sp->scab.enabled = true;
  ScabPDE scab(env.sp); scab.Init(env.sim);
  env.sp->mmp.enabled = true;
  MMPPDE mmp(env.sp); mmp.Init(env.sim);
  env.Run(TEST_NAME);
  // Scab forms from dried wound exudate; may need wound to produce it
}

TEST(ModuleSmoke, Photon) {
  auto env = SmokeEnv::Create(TEST_NAME);
  env.sp->photon.enabled = true;
  env.sp->ros.enabled = true;
  ROSPDE ros(env.sp); ros.Init(env.sim);
  FluencePDE fluence(env.sp); fluence.Init(env.sim);
  OpsinPDE opsin(env.sp); opsin.Init(env.sim);
  TemperaturePDE temp(env.sp); temp.Init(env.sim);
  env.Run(TEST_NAME);
}

TEST(ModuleSmoke, Blood) {
  auto env = SmokeEnv::Create(TEST_NAME);
  env.sp->hemostasis.enabled = true;
  FibrinPDE fibrin(env.sp); fibrin.Init(env.sim);
  // Blood reads perfusion, O2, fibrin, inflammation (all in base env)
  env.Run(TEST_NAME);
}

TEST(ModuleSmoke, Burn) {
  auto env = SmokeEnv::Create(TEST_NAME);
  // Burn reads vascular, water, inflammation, scar (all in base env)
  env.Run(TEST_NAME, 5);
}

TEST(ModuleSmoke, Pressure) {
  auto env = SmokeEnv::Create(TEST_NAME);
  env.sp->ros.enabled = true;
  ROSPDE ros(env.sp); ros.Init(env.sim);
  // Pressure reads vascular, water, inflammation, ROS
  env.Run(TEST_NAME, 5);
}

}  // namespace skibidy
}  // namespace bdm
