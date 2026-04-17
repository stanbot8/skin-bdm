// Aggregate include header for the hook API surface.
// Every module hook (source_hook.h, post_hook.h) reads from GridRegistry,
// fills or reads from VoxelSnapshot, and can raise flags via SignalBoard.
// Include this instead of repeating the three core headers.

#ifndef CORE_HOOK_API_H_
#define CORE_HOOK_API_H_

#include "core/grid_registry.h"
#include "core/signal_board.h"
#include "core/voxel_snapshot.h"

#endif  // CORE_HOOK_API_H_
