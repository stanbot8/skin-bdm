#ifndef VOXEL_SNAPSHOT_H_
#define VOXEL_SNAPSHOT_H_

#include <vector>
#include "core/field_names.h"
#include "core/material.h"
#include "core/pde.h"
#include "core/grid_registry.h"

namespace bdm {
namespace skibidy {

// ---------------------------------------------------------------------------
// VoxelSnapshot -- per-voxel read-once context.
// Populated once at the start of each voxel iteration. All hooks read common
// field values from here instead of calling GetConcentration() independently.
// Eliminates ~100-200k redundant grid reads per step.
//
// Design: reads are many-to-one (10 hooks read stratum), writes are
// one-to-few. VoxelSnapshot handles the read side. Hooks write directly
// to their cached grid pointers.
// ---------------------------------------------------------------------------
struct VoxelSnapshot {
  // Spatial context
  size_t idx;          // fine-grid flat index
  real_t x, y, z;     // world-space voxel center

  // Wound state
  bool in_wound;       // inside wound cylinder
  bool post_wound;     // step > wound_trigger_step
  uint64_t wound_age;  // step - wound_trigger_step (0 if pre-wound)

  // Commonly-read field values (shared by many hooks)
  real_t stratum;      // epidermal coverage (0..1, or -1 wound void)
  real_t barrier;      // StratumGate(stratum): re-epithelialization gate
  real_t vasc;         // vascular density [0..1]
  real_t o2;           // oxygen concentration
  real_t ph;           // pH alkalinity (structural grid, coarse-indexed)

  // Multi-resolution structural grid mapping
  size_t coarse_si;    // structural grid index for this fine voxel
  real_t coarse_w;     // volume correction factor for structural writes

  // Material layer
  uint8_t material_id = 0;                    // tissue type at this voxel
  const MaterialProperties* mat = nullptr;    // pointer into MaterialRegistry (never null after fill)
};

// ---------------------------------------------------------------------------
// SnapshotFiller -- fills VoxelSnapshot from registry + mask.
// Constructed once per fused-op invocation, then called per-voxel.
// ---------------------------------------------------------------------------
struct SnapshotFiller {
  DiffusionGrid* vasc_g = nullptr;
  DiffusionGrid* o2_g = nullptr;
  DiffusionGrid* strat_g = nullptr;
  DiffusionGrid* ph_g = nullptr;
  const GridContext* ctx = nullptr;
  const GridContext::WoundMaskData* mask = nullptr;
  const std::vector<size_t>* coarse_map = nullptr;
  bool coarse = false;
  real_t coarse_w = 1.0;
  bool post_wound = false;
  uint64_t wound_age = 0;

  // Material layer
  const MaterialRegistry* mat_reg = nullptr;
  real_t z_epi_top = 25.0;     // volume_z_cornified (epidermal boundary)
  real_t z_derm_bottom = -8.0; // dermal_z_reticular (dermis/hypodermis)

  void Init(const GridRegistry& reg, const GridContext& ctx_,
            const GridContext::WoundMaskData& mask_,
            const std::vector<size_t>& cm,
            bool c, real_t cw) {
    vasc_g = reg.Get(fields::kVascularId);
    o2_g = reg.Get(fields::kOxygenId);
    strat_g = reg.Get(fields::kStratumId);
    ph_g = reg.Get(fields::kPHId);
    ctx = &ctx_;
    mask = &mask_;
    coarse_map = &cm;
    coarse = c;
    coarse_w = cw;

    auto* sp = reg.Params();
    uint64_t wound_step = static_cast<uint64_t>(sp->wound.trigger_step);
    post_wound = (reg.Step() > wound_step);
    wound_age = post_wound ? reg.Step() - wound_step : 0;

    z_epi_top = sp->volume_z_cornified;
    z_derm_bottom = sp->dermal_z_reticular;

    mat_reg = &reg.Materials();
  }

  // Classify voxel into material by z-depth (skin default mapping).
  // Override this for non-skin tissue types.
  uint8_t ClassifyMaterial(real_t z) const {
    if (z >= 0) return material::kSkinEpidermis;
    if (z >= z_derm_bottom) return material::kSkinDermis;
    return material::kSkinHypodermis;
  }

  size_t CoarseSI(size_t fine) const {
    return (coarse && !coarse_map->empty()) ? (*coarse_map)[fine] : fine;
  }

  void FillDermal(VoxelSnapshot& s, size_t idx) const {
    s.idx = idx;
    s.x = ctx->X(idx);
    s.y = ctx->Y(idx);
    s.z = ctx->Z(idx);
    s.in_wound = mask->mask[idx];
    s.post_wound = post_wound;
    s.wound_age = wound_age;
    s.vasc = vasc_g ? vasc_g->GetConcentration(idx) : 0;
    s.o2 = o2_g ? o2_g->GetConcentration(idx) : 0;
    s.coarse_si = CoarseSI(idx);
    s.coarse_w = coarse_w;
    s.ph = ph_g ? ph_g->GetConcentration(s.coarse_si) : 0;
    s.stratum = 0;
    s.barrier = 0;
    s.material_id = ClassifyMaterial(s.z);
    s.mat = mat_reg ? &mat_reg->Get(s.material_id) : nullptr;
  }

  void FillEpi(VoxelSnapshot& s, size_t idx) const {
    s.idx = idx;
    s.x = ctx->X(idx);
    s.y = ctx->Y(idx);
    s.z = ctx->Z(idx);
    s.in_wound = true;  // always inside wound in epi_wound list
    s.post_wound = post_wound;
    s.wound_age = wound_age;
    s.stratum = strat_g ? strat_g->GetConcentration(idx) : 0;
    s.barrier = GridContext::StratumGate(s.stratum);
    s.vasc = vasc_g ? vasc_g->GetConcentration(idx) : 0;
    s.o2 = o2_g ? o2_g->GetConcentration(idx) : 0;
    s.coarse_si = CoarseSI(idx);
    s.coarse_w = coarse_w;
    s.ph = ph_g ? ph_g->GetConcentration(s.coarse_si) : 0;
    s.material_id = material::kSkinEpidermis;  // epi loop is always epidermis
    s.mat = mat_reg ? &mat_reg->Get(s.material_id) : nullptr;
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // VOXEL_SNAPSHOT_H_
