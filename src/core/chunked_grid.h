#ifndef CORE_CHUNKED_GRID_H_
#define CORE_CHUNKED_GRID_H_

// ChunkedGrid: sparse, chunk-based concentration storage with dirty tracking.
//
// Motivation: BioDynaMo's DiffusionGrid allocates a dense res^3 array for the
// entire simulation domain.  For wound healing, most of the domain is healthy
// tissue at equilibrium -- only the wound region and its diffusion front have
// interesting dynamics.  ChunkedGrid divides the domain into small cubic chunks
// (default 8^3 voxels) and only allocates chunks that deviate from a default
// value.  A per-chunk dirty flag lets solvers skip chunks whose concentrations
// haven't changed since the last step.
//
// This is a standalone data structure that wraps a flat voxel domain.  It does
// NOT replace DiffusionGrid -- it sits alongside it as an acceleration layer
// that the fused source/post ops can query for skip decisions.
//
// Design inspired by:
//   - Minecraft chunk system (16^3 sub-volumes, allocated on demand)
//   - OpenVDB sparse voxel trees (active/inactive tiles)
//   - AMR dirty-flagging (only recompute changed regions)
//
// Future: a ChunkedDiffusionGrid subclass could replace the dense solver with
// a chunk-aware FTCS stencil that skips uniform interior chunks entirely.

#include <cassert>
#include <cstring>
#include <memory>
#include <vector>

#include "core/real_t.h"

namespace bdm {
namespace skibidy {

// Compile-time chunk dimensions.  Must be a power of 2 for fast bit-shift
// indexing.  8^3 = 512 voxels per chunk is a sweet spot: small enough for
// fine-grained sparsity, large enough for SIMD-friendly inner loops.
struct ChunkTraits {
  static constexpr size_t kSide = 8;
  static constexpr size_t kSideBits = 3;  // log2(kSide)
  static constexpr size_t kVolume = kSide * kSide * kSide;  // 512
  static constexpr size_t kSideMask = kSide - 1;  // 0b111
};

// A single chunk: a dense 8^3 block of concentration values plus metadata.
struct Chunk {
  real_t data[ChunkTraits::kVolume];
  bool dirty = false;      // true if any voxel was modified since last reset
  bool allocated = false;  // false = all voxels at default_value (not stored)

  void Fill(real_t value) {
    for (size_t i = 0; i < ChunkTraits::kVolume; ++i) {
      data[i] = value;
    }
  }
};

// ChunkedGrid: sparse voxel storage over a cubic domain.
//
// The domain is divided into chunks of kSide^3 voxels.  Unallocated chunks
// are treated as uniform at `default_value_`.  Reading from an unallocated
// chunk returns the default.  Writing to an unallocated chunk triggers
// allocation (fill with default, then apply the write).
//
// Thread safety: same as DiffusionGrid -- callers are responsible for
// synchronization when writing.  The fused ops already serialize per-voxel
// writes within a single thread's chunk of the iteration range.
class ChunkedGrid {
 public:
  // Construct a chunked grid over a cubic domain.
  // resolution: number of voxels per axis (must be divisible by kSide)
  // default_value: concentration for unallocated chunks
  ChunkedGrid(size_t resolution, real_t default_value = 0.0)
      : resolution_(resolution),
        default_value_(default_value) {
    assert(resolution % ChunkTraits::kSide == 0 &&
           "Resolution must be divisible by chunk side length");
    chunks_per_axis_ = resolution / ChunkTraits::kSide;
    total_chunks_ = chunks_per_axis_ * chunks_per_axis_ * chunks_per_axis_;
    chunks_.resize(total_chunks_);
  }

  // --- Voxel access ---

  // Get concentration at flat voxel index.  Returns default if chunk not
  // allocated.  O(1) with decompose + bit-shift addressing.
  real_t Get(size_t voxel_idx) const {
    size_t vx, vy, vz;
    DecomposeVoxel(voxel_idx, vx, vy, vz);
    size_t ci = ChunkIndexFromCoords(vx, vy, vz);
    if (!chunks_[ci]) return default_value_;
    return chunks_[ci]->data[LocalIndex(vx, vy, vz)];
  }

  // Set concentration at flat voxel index.  Allocates chunk on first write.
  void Set(size_t voxel_idx, real_t value) {
    size_t vx, vy, vz;
    DecomposeVoxel(voxel_idx, vx, vy, vz);
    size_t ci = ChunkIndexFromCoords(vx, vy, vz);
    EnsureAllocated(ci);
    chunks_[ci]->data[LocalIndex(vx, vy, vz)] = value;
    chunks_[ci]->dirty = true;
  }

  // Add delta to concentration at flat voxel index.  Allocates on first write.
  void Add(size_t voxel_idx, real_t delta) {
    size_t vx, vy, vz;
    DecomposeVoxel(voxel_idx, vx, vy, vz);
    size_t ci = ChunkIndexFromCoords(vx, vy, vz);
    EnsureAllocated(ci);
    chunks_[ci]->data[LocalIndex(vx, vy, vz)] += delta;
    chunks_[ci]->dirty = true;
  }

  // --- Chunk-level queries ---

  // Is this chunk allocated (has any voxel ever been written)?
  bool IsAllocated(size_t chunk_idx) const {
    return chunks_[chunk_idx] != nullptr;
  }

  // Has this chunk been modified since the last ClearDirty()?
  bool IsDirty(size_t chunk_idx) const {
    return chunks_[chunk_idx] && chunks_[chunk_idx]->dirty;
  }

  // Reset dirty flags on all chunks.  Call after the diffusion solver step.
  void ClearAllDirty() {
    for (auto& c : chunks_) {
      if (c) c->dirty = false;
    }
  }

  // Reset dirty flag on a single chunk.
  void ClearDirty(size_t chunk_idx) {
    if (chunks_[chunk_idx]) chunks_[chunk_idx]->dirty = false;
  }

  // --- Bulk operations ---

  // Copy all concentrations from a flat array (e.g. DiffusionGrid::c1_).
  // Only allocates chunks that differ from default_value_.
  void ImportFromFlat(const real_t* flat, size_t total_voxels) {
    assert(total_voxels == resolution_ * resolution_ * resolution_);
    for (size_t ci = 0; ci < total_chunks_; ++ci) {
      // Check if chunk has any non-default voxels
      size_t base = ChunkBaseVoxel(ci);
      bool all_default = true;
      for (size_t li = 0; li < ChunkTraits::kVolume && all_default; ++li) {
        size_t vi = FlatFromChunkLocal(ci, li);
        if (vi < total_voxels && flat[vi] != default_value_) {
          all_default = false;
        }
      }
      if (all_default) {
        chunks_[ci].reset();  // Free memory for uniform chunks
        continue;
      }
      EnsureAllocated(ci);
      for (size_t li = 0; li < ChunkTraits::kVolume; ++li) {
        size_t vi = FlatFromChunkLocal(ci, li);
        chunks_[ci]->data[li] = (vi < total_voxels) ? flat[vi] : default_value_;
      }
      chunks_[ci]->dirty = true;
    }
  }

  // Export all concentrations to a flat array.
  void ExportToFlat(real_t* flat, size_t total_voxels) const {
    assert(total_voxels == resolution_ * resolution_ * resolution_);
    for (size_t ci = 0; ci < total_chunks_; ++ci) {
      for (size_t li = 0; li < ChunkTraits::kVolume; ++li) {
        size_t vi = FlatFromChunkLocal(ci, li);
        if (vi < total_voxels) {
          flat[vi] = chunks_[ci] ? chunks_[ci]->data[li] : default_value_;
        }
      }
    }
  }

  // --- Statistics ---

  size_t GetNumAllocatedChunks() const {
    size_t count = 0;
    for (auto& c : chunks_) {
      if (c) ++count;
    }
    return count;
  }

  size_t GetNumDirtyChunks() const {
    size_t count = 0;
    for (auto& c : chunks_) {
      if (c && c->dirty) ++count;
    }
    return count;
  }

  size_t GetTotalChunks() const { return total_chunks_; }
  size_t GetChunksPerAxis() const { return chunks_per_axis_; }
  size_t GetResolution() const { return resolution_; }
  real_t GetDefaultValue() const { return default_value_; }

  // --- Index math ---

  // Convert flat voxel index to chunk index (decomposes internally).
  size_t ChunkIndex(size_t voxel_idx) const {
    size_t vx, vy, vz;
    DecomposeVoxel(voxel_idx, vx, vy, vz);
    return ChunkIndexFromCoords(vx, vy, vz);
  }

  // Convert voxel coordinates to chunk index (fast bit-shift path).
  size_t ChunkIndexFromCoords(size_t vx, size_t vy, size_t vz) const {
    size_t cx = vx >> ChunkTraits::kSideBits;
    size_t cy = vy >> ChunkTraits::kSideBits;
    size_t cz = vz >> ChunkTraits::kSideBits;
    return cx + cy * chunks_per_axis_ + cz * chunks_per_axis_ * chunks_per_axis_;
  }

  // Local index within chunk from voxel coordinates (bit-mask extraction).
  static size_t LocalIndex(size_t vx, size_t vy, size_t vz) {
    size_t lx = vx & ChunkTraits::kSideMask;
    size_t ly = vy & ChunkTraits::kSideMask;
    size_t lz = vz & ChunkTraits::kSideMask;
    return lx + ly * ChunkTraits::kSide + lz * ChunkTraits::kSide * ChunkTraits::kSide;
  }

  // Decompose flat voxel index into (vx, vy, vz) coordinates.
  void DecomposeVoxel(size_t voxel_idx, size_t& vx, size_t& vy,
                      size_t& vz) const {
    vx = voxel_idx % resolution_;
    vy = (voxel_idx / resolution_) % resolution_;
    vz = voxel_idx / (resolution_ * resolution_);
  }

  // Get chunk 3D coordinates from chunk index.
  void DecomposeChunk(size_t chunk_idx, size_t& cx, size_t& cy,
                      size_t& cz) const {
    cx = chunk_idx % chunks_per_axis_;
    cy = (chunk_idx / chunks_per_axis_) % chunks_per_axis_;
    cz = chunk_idx / (chunks_per_axis_ * chunks_per_axis_);
  }

  // Iterate over voxel indices belonging to a chunk.
  // Calls f(flat_voxel_idx, local_idx) for each voxel in the chunk.
  template <typename F>
  void ForEachVoxelInChunk(size_t chunk_idx, F&& f) const {
    size_t cx, cy, cz;
    DecomposeChunk(chunk_idx, cx, cy, cz);
    size_t bx = cx * ChunkTraits::kSide;
    size_t by = cy * ChunkTraits::kSide;
    size_t bz = cz * ChunkTraits::kSide;
    size_t li = 0;
    for (size_t lz = 0; lz < ChunkTraits::kSide; ++lz) {
      for (size_t ly = 0; ly < ChunkTraits::kSide; ++ly) {
        for (size_t lx = 0; lx < ChunkTraits::kSide; ++lx, ++li) {
          size_t vx = bx + lx;
          size_t vy = by + ly;
          size_t vz = bz + lz;
          if (vx < resolution_ && vy < resolution_ && vz < resolution_) {
            size_t flat = vx + vy * resolution_ + vz * resolution_ * resolution_;
            f(flat, li);
          }
        }
      }
    }
  }

 private:
  size_t resolution_;        // voxels per axis
  size_t chunks_per_axis_;   // chunks per axis
  size_t total_chunks_;      // total number of chunk slots
  real_t default_value_;     // concentration for unallocated chunks

  std::vector<std::unique_ptr<Chunk>> chunks_;

  void EnsureAllocated(size_t chunk_idx) {
    if (!chunks_[chunk_idx]) {
      chunks_[chunk_idx] = std::make_unique<Chunk>();
      chunks_[chunk_idx]->Fill(default_value_);
      chunks_[chunk_idx]->allocated = true;
    }
  }

  // Base flat voxel index for the first voxel in a chunk (not sequential
  // in flat order due to 3D layout, but useful for iteration).
  size_t ChunkBaseVoxel(size_t chunk_idx) const {
    size_t cx, cy, cz;
    DecomposeChunk(chunk_idx, cx, cy, cz);
    size_t bx = cx * ChunkTraits::kSide;
    size_t by = cy * ChunkTraits::kSide;
    size_t bz = cz * ChunkTraits::kSide;
    return bx + by * resolution_ + bz * resolution_ * resolution_;
  }

  // Convert chunk index + local index to flat voxel index.
  size_t FlatFromChunkLocal(size_t chunk_idx, size_t local_idx) const {
    size_t cx, cy, cz;
    DecomposeChunk(chunk_idx, cx, cy, cz);
    size_t lx = local_idx % ChunkTraits::kSide;
    size_t ly = (local_idx / ChunkTraits::kSide) % ChunkTraits::kSide;
    size_t lz = local_idx / (ChunkTraits::kSide * ChunkTraits::kSide);
    size_t vx = cx * ChunkTraits::kSide + lx;
    size_t vy = cy * ChunkTraits::kSide + ly;
    size_t vz = cz * ChunkTraits::kSide + lz;
    return vx + vy * resolution_ + vz * resolution_ * resolution_;
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // CORE_CHUNKED_GRID_H_
