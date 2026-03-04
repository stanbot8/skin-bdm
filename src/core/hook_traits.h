#ifndef HOOK_TRAITS_H_
#define HOOK_TRAITS_H_

#include <type_traits>
#include "core/voxel_snapshot.h"
#include "core/signal_board.h"
#include "core/grid_registry.h"

namespace bdm {
namespace skibidy {

// ---------------------------------------------------------------------------
// Detection traits -- compile-time check for which Apply methods a hook
// provides. The engine uses if constexpr to dispatch only to hooks that
// implement the relevant method. Dead code elimination removes the rest.
//
// Inspired by Unreal Engine's UCLASS/UFUNCTION reflection: declare
// capabilities, engine dispatches automatically.
// ---------------------------------------------------------------------------

// Does Hook have ApplyDermal(const VoxelSnapshot&, SignalBoard&)?
template<typename T, typename = void>
struct has_apply_dermal : std::false_type {};

template<typename T>
struct has_apply_dermal<T, std::void_t<decltype(
    std::declval<T&>().ApplyDermal(
        std::declval<const VoxelSnapshot&>(),
        std::declval<SignalBoard&>()))>> : std::true_type {};

template<typename T>
inline constexpr bool has_apply_dermal_v = has_apply_dermal<T>::value;

// Does Hook have ApplyEpiWound(const VoxelSnapshot&, SignalBoard&)?
template<typename T, typename = void>
struct has_apply_epi_wound : std::false_type {};

template<typename T>
struct has_apply_epi_wound<T, std::void_t<decltype(
    std::declval<T&>().ApplyEpiWound(
        std::declval<const VoxelSnapshot&>(),
        std::declval<SignalBoard&>()))>> : std::true_type {};

template<typename T>
inline constexpr bool has_apply_epi_wound_v = has_apply_epi_wound<T>::value;

// Does Hook have ApplyDecay(const WoundMaskData&, SidxFn, real_t, real_t)?
// Decay hooks own their iteration -- engine calls them after the main loops.
template<typename T, typename = void>
struct has_apply_decay : std::false_type {};

template<typename T>
struct has_apply_decay<T, std::void_t<decltype(
    std::declval<T&>().ApplyDecay(
        std::declval<const GridContext::WoundMaskData&>(),
        std::declval<real_t>()))>> : std::true_type {};

template<typename T>
inline constexpr bool has_apply_decay_v = has_apply_decay<T>::value;

// ---------------------------------------------------------------------------
// Dispatch helpers -- call a hook's method only if it exists and is active.
// Compile-time elimination: if a hook lacks ApplyDermal, the entire call
// is removed by the compiler (zero overhead, no vtable, no branch).
// ---------------------------------------------------------------------------

template<typename Hook>
inline void DispatchDermal(Hook& hook, const VoxelSnapshot& snap,
                           SignalBoard& sig) {
  if constexpr (has_apply_dermal_v<Hook>) {
    if (hook.active) hook.ApplyDermal(snap, sig);
  }
}

template<typename Hook>
inline void DispatchEpiWound(Hook& hook, const VoxelSnapshot& snap,
                              SignalBoard& sig) {
  if constexpr (has_apply_epi_wound_v<Hook>) {
    if (hook.active) hook.ApplyEpiWound(snap, sig);
  }
}

template<typename Hook>
inline void DispatchDecay(Hook& hook,
                          const GridContext::WoundMaskData& mask,
                          real_t dt) {
  if constexpr (has_apply_decay_v<Hook>) {
    if (hook.active) hook.ApplyDecay(mask, dt);
  }
}

}  // namespace skibidy
}  // namespace bdm

#endif  // HOOK_TRAITS_H_
