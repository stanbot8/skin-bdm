#ifndef FIELDS_PDE_H_
#define FIELDS_PDE_H_

#include <chrono>
#include <cmath>
#include <iostream>

#include "biodynamo.h"
#include "core/field_names.h"
#include "infra/sim_param.h"

namespace bdm {
namespace skibidy {

class CompositeField;  // forward declaration

// Precomputed grid metadata + wound geometry for voxel iteration.
// Constructed once per loop from any DiffusionGrid (all share the same mesh).
struct GridContext {
  size_t res;
  size_t n;
  real_t lo;
  real_t box_len;
  real_t cx, cy, r2;

  GridContext(DiffusionGrid* grid, const SimParam* sp)
      : res(grid->GetResolution()),
        n(grid->GetNumBoxes()),
        box_len(grid->GetBoxLength()),
        cx(sp->wound_center_x),
        cy(sp->wound_center_y),
        r2(sp->wound_radius * sp->wound_radius) {
    auto dims = grid->GetDimensions();
    lo = static_cast<real_t>(dims[0]);
  }

  real_t X(size_t idx) const {
    return lo + (idx % res) * box_len + box_len / 2.0;
  }
  real_t Y(size_t idx) const {
    return lo + ((idx / res) % res) * box_len + box_len / 2.0;
  }
  real_t Z(size_t idx) const {
    return lo + (idx / (res * res)) * box_len + box_len / 2.0;
  }

  bool InWound(real_t x, real_t y) const {
    real_t dx = x - cx, dy = y - cy;
    return dx * dx + dy * dy <= r2;
  }

  // Map a fine-grid voxel index to its corresponding coarse-grid voxel index
  // via world coordinates. Used when fused ops on the fine grid need to read
  // or write structural fields on the coarser grid.
  static size_t CoarseIndex(size_t fine_idx, const GridContext& fine_ctx,
                            DiffusionGrid* coarse_grid) {
    Real3 pos = {fine_ctx.X(fine_idx), fine_ctx.Y(fine_idx),
                 fine_ctx.Z(fine_idx)};
    return coarse_grid->GetBoxIndex(pos);
  }

  // Precomputed wound mask + layer map for fused voxel passes.
  // Computed once at first use; eliminates per-voxel InWound() + Z() checks.
  enum Layer : uint8_t { kDermal = 0, kEpidermal = 1, kAbove = 2 };
  enum DermalSubLayer : uint8_t {
    kSubPapillary = 0, kSubReticular = 1, kSubHypodermis = 2
  };
  struct WoundMaskData {
    std::vector<uint8_t> mask;              // 1 = inside wound cylinder
    std::vector<uint8_t> layer;             // kDermal / kEpidermal / kAbove
    std::vector<uint8_t> dermal_sub_layer;  // per-voxel sub-layer tag
    // Compact index lists: iterate only relevant voxels instead of all n.
    // At 10^3 grid this is ~3x fewer iterations; at 50^3 it's ~30x fewer.
    std::vector<size_t> dermal_all;       // all dermal voxels (O2/water pin)
    std::vector<size_t> dermal_wound;     // dermal voxels inside wound cylinder
    std::vector<size_t> epi_wound;        // epidermal wound voxels
    // Dermal sub-layer lists (papillary may be empty at low resolution)
    std::vector<size_t> papillary_all;
    std::vector<size_t> papillary_wound;
    std::vector<size_t> reticular_all;
    std::vector<size_t> reticular_wound;
    std::vector<size_t> hypodermis_all;
    std::vector<size_t> hypodermis_wound;
  };
  static WoundMaskData ComputeWoundMask(const GridContext& ctx, real_t z_max,
                                         const SimParam* sp = nullptr) {
    WoundMaskData data;
    data.mask.resize(ctx.n);
    data.layer.resize(ctx.n);
    data.dermal_sub_layer.resize(ctx.n, 0);

    real_t z_papillary = sp ? sp->dermal_z_papillary : -2.0;
    real_t z_reticular = sp ? sp->dermal_z_reticular : -8.0;

    for (size_t idx = 0; idx < ctx.n; idx++) {
      bool in_wound = ctx.InWound(ctx.X(idx), ctx.Y(idx));
      data.mask[idx] = in_wound ? 1 : 0;
      real_t z = ctx.Z(idx);
      if (z < 0) {
        data.layer[idx] = kDermal;
        data.dermal_all.push_back(idx);
        if (in_wound) data.dermal_wound.push_back(idx);

        // Sub-layer classification
        if (z >= z_papillary) {
          data.dermal_sub_layer[idx] = kSubPapillary;
          data.papillary_all.push_back(idx);
          if (in_wound) data.papillary_wound.push_back(idx);
        } else if (z >= z_reticular) {
          data.dermal_sub_layer[idx] = kSubReticular;
          data.reticular_all.push_back(idx);
          if (in_wound) data.reticular_wound.push_back(idx);
        } else {
          data.dermal_sub_layer[idx] = kSubHypodermis;
          data.hypodermis_all.push_back(idx);
          if (in_wound) data.hypodermis_wound.push_back(idx);
        }
      } else if (z <= z_max) {
        data.layer[idx] = kEpidermal;
        if (in_wound) data.epi_wound.push_back(idx);
      } else {
        data.layer[idx] = kAbove;
      }
    }
    return data;
  }

  // Profile helpers (pure math).
  static real_t ExpDecay(real_t z, real_t base, real_t length) {
    if (z < 0) return base;  // dermal voxels: no supraphysiological values
    return base * std::exp(-z / length);
  }
  static real_t Sigmoid(real_t z, real_t lo, real_t hi, real_t mid, real_t k) {
    return lo + (hi - lo) / (1.0 + std::exp(-(z - mid) / k));
  }
  static real_t StratumGate(real_t stratum_val) {
    return std::min(static_cast<real_t>(1.0), stratum_val / 3.0);
  }

  // Pin dermal voxels (z < 0) to base_conc * local vascular perfusion.
  static void PinDermal(DiffusionGrid* grid, DiffusionGrid* vasc_grid,
                        const GridContext& ctx, real_t base_conc) {
    for (size_t idx = 0; idx < ctx.n; idx++) {
      if (ctx.Z(idx) >= 0) continue;
      real_t target = base_conc * vasc_grid->GetConcentration(idx);
      real_t delta = target - grid->GetConcentration(idx);
      if (std::abs(delta) > 1e-10) {
        grid->ChangeConcentrationBy(idx, delta);
      }
    }
  }

  static void ZeroWoundDermis(DiffusionGrid* grid, const GridContext& ctx) {
    for (size_t idx = 0; idx < ctx.n; idx++) {
      if (ctx.Z(idx) >= 0) continue;
      if (ctx.InWound(ctx.X(idx), ctx.Y(idx))) Zero(grid, idx);
    }
  }

  static void ZeroWoundEpidermis(DiffusionGrid* grid, const GridContext& ctx) {
    for (size_t idx = 0; idx < ctx.n; idx++) {
      if (ctx.Z(idx) < 0) continue;
      if (ctx.InWound(ctx.X(idx), ctx.Y(idx))) Zero(grid, idx);
    }
  }

  static void Zero(DiffusionGrid* grid, size_t idx) {
    real_t val = grid->GetConcentration(idx);
    if (std::abs(val) > 1e-10) {
      grid->ChangeConcentrationBy(idx, -val);
    }
  }
};

// ---------------------------------------------------------------------------
// ScaledGrid -- resolution-aware write wrapper for DIFFUSING fields.
//
// Use ONLY for fields with diffusion > 0 (cytokines, signaling molecules).
// Do NOT use for non-diffusing structural fields (collagen, fibronectin,
// elastin, biofilm) where concentration is a local density property.
//
// Agent-level writes inject a fixed physical mass of substance into one
// voxel.  At finer grids the voxel volume shrinks, so the concentration
// delta must increase proportionally to maintain the same mass:
//
//   mass = concentration * voxel_vol
//   → concentration = mass / voxel_vol = rate * (ref_vol / actual_vol)
//
// Per-voxel sources (fused loops) do NOT need scaling because the number
// of wound voxels grows as 1/voxel_vol, keeping total mass constant.
//
// Usage:
//   ScaledGrid infl(rm->GetDiffusionGrid(fields::kInflammation), sp);
//   infl.AgentDeposit(idx, sp->immune_cytokine_rate);   // auto-scaled
// ---------------------------------------------------------------------------
struct ScaledGrid {
  DiffusionGrid* grid;
  real_t agent_factor;   // ref_vol / actual_vol  (>= 1 at fine grids)

  ScaledGrid() : grid(nullptr), agent_factor(1) {}

  ScaledGrid(DiffusionGrid* g, const SimParam* sp) : grid(g) {
    real_t box = g->GetBoxLength();
    real_t ref = sp->ref_box_length;
    agent_factor = (ref * ref * ref) / (box * box * box);
    if (sp->debug_scaled_grid) {
      static bool printed = false;
      if (!printed) {
        std::cout << "[scaled_grid] box_len=" << box << " ref=" << ref
                  << " agent_factor=" << agent_factor << std::endl;
        printed = true;
      }
    }
  }

  // Per-agent deposit: fixed mass injection scaled to voxel volume.
  void AgentDeposit(size_t idx, real_t amount) const {
    if (amount > 1e-12)
      grid->ChangeConcentrationBy(idx, amount * agent_factor);
  }

  // Per-agent removal: fixed mass removal, capped by current concentration.
  void AgentRemove(size_t idx, real_t amount) const {
    real_t current = grid->GetConcentration(idx);
    real_t rem = std::min(current, amount * agent_factor);
    if (rem > 1e-10)
      grid->ChangeConcentrationBy(idx, -rem);
  }

  // Raw read (no scaling needed for concentration queries).
  real_t Get(size_t idx) const { return grid->GetConcentration(idx); }
  size_t Index(const Real3& pos) const { return grid->GetBoxIndex(pos); }
};

// Abstract base for a composable PDE channel.
// Each subclass wraps one BioDynaMo DiffusionGrid and encapsulates its
// initialization, per-step source term, and wound disruption behavior.
struct PDE {
  virtual ~PDE() = default;

  virtual const char* GetName() const = 0;
  virtual int GetId() const = 0;

  // Called once at setup. Define + initialize the BioDynaMo DiffusionGrid.
  virtual void Init(Simulation* sim) = 0;

  // Per-step source term. Reads other channels via `fields` for coupling.
  // Default: no-op (static/prescribed channels).
  virtual void ApplySource(Simulation* sim, const CompositeField& fields) {
    (void)sim;
    (void)fields;
  }

  // Wound disruption. Default: no-op (most fields start at 0).
  // Override for fields that need zeroing or custom wound response.
  virtual void ApplyWound(Simulation* sim, real_t cx, real_t cy, real_t r) {
    (void)sim; (void)cx; (void)cy; (void)r;
  }

  // Helper: zero all voxels inside wound cylinder.
  void ZeroInWound(Simulation* sim) {
    auto* grid = Grid(sim);
    auto* sp = sim->GetParam()->Get<SimParam>();
    GridContext ctx(grid, sp);
    for (size_t idx = 0; idx < ctx.n; idx++) {
      if (ctx.InWound(ctx.X(idx), ctx.Y(idx))) {
        GridContext::Zero(grid, idx);
      }
    }
  }

  // Cached grid accessor.
  DiffusionGrid* Grid(Simulation* sim) const {
    if (!grid_) {
      grid_ = sim->GetResourceManager()->GetDiffusionGrid(GetName());
    }
    return grid_;
  }

  // Register this field's DiffusionGrid with BioDynaMo. Resolution is always
  // read from SimParam::grid_resolution so every grid uses the same mesh.
  void DefineGrid(Simulation* sim, real_t diffusion, real_t decay) {
    auto* sp = sim->GetParam()->Get<SimParam>();
    ModelInitializer::DefineSubstance(GetId(), GetName(), diffusion, decay,
                                      sp->grid_resolution);
  }

  // Like DefineGrid but uses the coarser structural resolution when set.
  // Structural fields (D=0) store local densities that don't need fine
  // spatial resolution. Running them at a coarser grid saves 8x memory
  // at half resolution. Falls back to grid_resolution when structural = 0.
  void DefineStructuralGrid(Simulation* sim, real_t diffusion, real_t decay) {
    auto* sp = sim->GetParam()->Get<SimParam>();
    int res = (sp->grid_resolution_structural > 0)
                  ? sp->grid_resolution_structural
                  : sp->grid_resolution;
    ModelInitializer::DefineSubstance(GetId(), GetName(), diffusion, decay,
                                      res);
  }

  // For prescribed fields (diffusion=0, decay=0): skip the BDM diffusion
  // solver by setting a time step so large that IntegrateTimeAsynchronously
  // never triggers Step().  Biologically: these profiles are maintained by
  // active cellular processes (ion pumps, paracrine signaling), not by
  // physical diffusion.  Source terms still run via CompositeFieldOp.
  void MarkPrescribed(Simulation* sim) {
    Grid(sim)->SetTimeStep(1e30);
  }

 protected:
  mutable DiffusionGrid* grid_ = nullptr;
};

// Concrete base for agent-written fields (no profile, no source, no wound response).
// Subclass only needs to provide name/id/params via constructor.
struct SimplePDE : public PDE {
  SimplePDE(const char* name, int id, real_t diffusion, real_t decay)
      : name_(name), id_(id), diffusion_(diffusion), decay_(decay) {}

  const char* GetName() const override { return name_; }
  int GetId() const override { return id_; }

  void Init(Simulation* sim) override {
    DefineGrid(sim, diffusion_, decay_);
    if (diffusion_ == 0 && decay_ == 0) {
      MarkPrescribed(sim);  // no physics to solve, skip diffusion step
    }
  }

 private:
  const char* name_;
  int id_;
  real_t diffusion_;
  real_t decay_;
};

// Same as SimplePDE but registers on the structural (coarser) grid.
// Use for non-diffusing local density fields: collagen, fibronectin, elastin,
// fibrin, scar, dermis, biofilm, hyaluronan, pH, tumor.
struct SimpleStructuralPDE : public PDE {
  SimpleStructuralPDE(const char* name, int id, real_t diffusion, real_t decay)
      : name_(name), id_(id), diffusion_(diffusion), decay_(decay) {}

  const char* GetName() const override { return name_; }
  int GetId() const override { return id_; }

  void Init(Simulation* sim) override {
    DefineStructuralGrid(sim, diffusion_, decay_);
    if (diffusion_ == 0 && decay_ == 0) {
      MarkPrescribed(sim);
    }
  }

 private:
  const char* name_;
  int id_;
  real_t diffusion_;
  real_t decay_;
};

// ---------------------------------------------------------------------------
// PerfTimer -- lightweight wall-clock timer for debug_perf profiling.
//
// Usage:
//   PerfTimer timer(sp->debug_perf);
//   // ... do work ...
//   timer.Print("fused_source");  // prints elapsed ms if debug_perf=true
// ---------------------------------------------------------------------------
struct PerfTimer {
  bool enabled;
  std::chrono::high_resolution_clock::time_point t0;

  explicit PerfTimer(bool on) : enabled(on) {
    if (enabled) t0 = std::chrono::high_resolution_clock::now();
  }
  void Print(const char* label) const {
    if (!enabled) return;
    auto t1 = std::chrono::high_resolution_clock::now();
    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    std::cout << "[perf] " << label << " " << ms << "ms" << std::endl;
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // FIELDS_PDE_H_
