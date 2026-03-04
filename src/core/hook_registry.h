#ifndef HOOK_REGISTRY_H_
#define HOOK_REGISTRY_H_

#include <tuple>

// Source hooks
#include "angiogenesis/source_hook.h"
#include "bioelectric/source_hook.h"
#include "blood/source_hook.h"
#include "burn/source_hook.h"
#include "dermis/source_hook.h"
#include "elastin/source_hook.h"
#include "fibronectin/source_hook.h"
#include "glucose/source_hook.h"
#include "hemostasis/source_hook.h"
#include "scab/source_hook.h"
#include "inflammation/source_hook.h"
#include "lactate/source_hook.h"
#include "lymphatic/source_hook.h"
#include "mechanotransduction/source_hook.h"
#include "neuropathy/source_hook.h"
#include "perfusion/source_hook.h"
#include "photon/source_hook.h"
#include "ph/source_hook.h"
#include "pressure/source_hook.h"
#include "ros/source_hook.h"
#include "temperature/source_hook.h"

// Post hooks
#include "biofilm/post_hook.h"
#include "diabetic/post_hook.h"
#include "fibroblast/post_hook.h"
#include "fibronectin/post_hook.h"
#include "glucose/post_hook.h"
#include "lymphatic/post_hook.h"
#include "mmp/post_hook.h"
#include "neuropathy/post_hook.h"
#include "perfusion/post_hook.h"
#include "ros/post_hook.h"
#include "scar/post_hook.h"
#include "senescence/post_hook.h"

#include "core/hook_traits.h"
#include "core/grid_registry.h"
#include "core/signal_board.h"
#include "core/study_hooks_gen.h"

namespace bdm {
namespace skibidy {

// ---------------------------------------------------------------------------
// Hook tuple types -- two-tier registration (Unreal Engine/Project plugin pattern).
//
// Engine hooks (shared biology): include above, add to the tuple below.
// Study hooks (disease-specific): create studies/{name}/modules/study_hooks.h,
//   run gen_study_hooks.py, rebuild. No core files need to change.
// ---------------------------------------------------------------------------

// Source hooks: dispatched in dermal and epidermal wound loops.
// Ordering: angiogenesis first (publishes SignalBoard signals for downstream),
// perfusion second (O2/water physics), then remaining modules alphabetically.
using SourceHooks = std::tuple<
    AngioSourceHook,
    PerfusionSourceHook,
    WoundInflSourceHook,
    PHSourceHook,
    ROSSourceHook,
    TemperatureSourceHook,
    GlucoseSourceHook,
    LactateSourceHook,
    MechanoSourceHook,
    LymphaticSourceHook,
    BioelectricSourceHook,
    NeuropathySourceHook,
    BloodSourceHook,
    BurnSourceHook,
    PressureSourceHook,
    PhotonSourceHook,
    DermisSourceHook>;

// Decay hooks: self-contained iteration (called after main loops).
// Each hook owns its own coarse_map and iterates wound voxels internally.
using DecayHooks = std::tuple<
    FibronectinDecayHook,
    ElastinDecayHook,
    HemostasisSourceHook,
    ScabDecayHook>;

// Post hooks: dispatched in epidermal wound loop only.
// Dependency ordering preserved: biofilm -> diabetic inflammation -> MMP ->
// fibroblast clearance -> lymphatic/perfusion -> fibronectin -> scar ->
// AGE -> senescence -> neuropathy -> ROS -> RA.
using PostHooks = std::tuple<
    BiofilmPostHook,
    DiabeticPostHook,
    MMPPostHook,
    FibroblastPostHook,
    LymphaticPostHook,
    PerfusionPostHook,
    FibronectinPostHook,
    ScarPostHook,
    AGEPostHook,
    SenescencePostHook,
    NeuropathyPostHook,
    ROSPostHook>;

// ---------------------------------------------------------------------------
// Tuple helpers -- fold-expression dispatch over hook collections.
// ---------------------------------------------------------------------------

// Initialize all hooks in a tuple.
template<typename Tuple>
inline void InitHooks(Tuple& hooks, const GridRegistry& reg,
                      SignalBoard& sig) {
  std::apply([&](auto&... h) { (h.Init(reg, sig), ...); }, hooks);
}

// Check if any hook in a tuple is active.
template<typename Tuple>
inline bool AnyActive(const Tuple& hooks) {
  bool any = false;
  std::apply([&](const auto&... h) {
    ((any = any || h.active), ...);
  }, hooks);
  return any;
}

}  // namespace skibidy
}  // namespace bdm

#endif  // HOOK_REGISTRY_H_
