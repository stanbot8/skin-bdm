#ifndef DIFFERENTIATION_H_
#define DIFFERENTIATION_H_

#include <cmath>

#include "tissue/keratinocyte.h"
#include "core/field_names.h"
#include "infra/sim_param.h"
#include "infra/util.h"
#include "infra/volume.h"
#include "scar/scar_op.h"

namespace bdm {
namespace skibidy {

// ---------------------------------------------------------------------------
// Behavior: assign epidermal stratum from calcium concentration and height.
// Calcium is the primary differentiation signal (Grabe & Neuber 2005);
// z-height acts as reinforcement.
// Note: BioDynaMo models cells as spheres, so morphological flattening
// (squamous shape) is not represented geometrically.  Stratum identity is
// tracked via the stratum_ data member for visualization and analysis.
// ---------------------------------------------------------------------------
struct Differentiation : public Behavior {
  BDM_BEHAVIOR_HEADER(Differentiation, Behavior, 1);

  Differentiation() { AlwaysCopyToNew(); }
  virtual ~Differentiation() {}

  void Run(Agent* agent) override {
    auto* cell = dynamic_cast<Keratinocyte*>(agent);
    if (!cell) return;

    auto pos = cell->GetPosition();

    // Basement membrane: clamp cells above z=0
    if (pos[2] < 0) {
      cell->SetPosition({pos[0], pos[1], 0});
      pos = cell->GetPosition();
    }

    // Stem cell anchoring: basal stem cells express integrins (alpha6-beta4)
    // that maintain adhesion to the basement membrane, preventing mechanical
    // displacement from the basal niche.
    if (cell->IsStem()) {
      real_t max_z = cell->GetDiameter() * 0.5;
      if (pos[2] > max_z) {
        cell->SetPosition({pos[0], pos[1], max_z});
        pos = cell->GetPosition();
      }
    }

    auto* sim = Simulation::GetActive();
    auto* sp = sim->GetParam()->Get<SimParam>();

    // Wound wall: soft repulsion at wound boundary.
    // Cells can poke slightly past the edge but get pushed back gently.
    if (sp->wound_enabled) {
      real_t dx = pos[0] - sp->wound_center_x;
      real_t dy = pos[1] - sp->wound_center_y;
      real_t dist = std::sqrt(dx * dx + dy * dy);
      real_t r = sp->wound_radius;
      real_t overshoot = dist - r;
      if (overshoot > 0 && dist > 1e-6) {
        // Push back 30% of the overshoot per step (soft spring)
        real_t pushback = 0.3 * overshoot;
        pos[0] -= dx / dist * pushback;
        pos[1] -= dy / dist * pushback;
        cell->SetPosition({pos[0], pos[1], pos[2]});
      }
    }

    // Read calcium concentration at cell position.
    auto* rm = sim->GetResourceManager();
    auto* ca_grid = rm->GetDiffusionGrid(fields::kCalcium);
    Real3 qpos = ClampToBounds(cell->GetPosition(), sim->GetParam());
    real_t ca = ca_grid->GetValue(qpos);
    real_t z = cell->GetPosition()[2];

    // Stratum assignment: volume (z-based) gives coarse stratum,
    // calcium provides biological nuance near boundaries.
    auto* vm = VolumeManager::Get();
    Stratum vol_stratum = vm->GetStratumAt(cell->GetPosition());

    // Calcium can override volume near boundaries:
    // - Low calcium keeps cells basal even if z is above the volume boundary
    // - High calcium promotes differentiation even if z is near the boundary
    Stratum new_stratum;
    if (ca < sp->ca_spinous_threshold && z < sp->spinous_threshold) {
      new_stratum = kBasal;
    } else {
      new_stratum = vol_stratum;
    }

    // UWYN: continuum-only strata don't track individual agents.
    // The cell's upward transit is represented by the continuum fields.
    if (!vm->AreAgentsEnabled(new_stratum)) {
      cell->RemoveFromSimulation();
      return;
    }

    cell->SetStratum(new_stratum);
    cell->IncrementStratumSteps();

    // Stratum field recovery: write agent stratum back into the
    // continuum Stratum grid to visually close the wound breach.
    if (sp->wound_enabled) {
      real_t dx_w = pos[0] - sp->wound_center_x;
      real_t dy_w = pos[1] - sp->wound_center_y;
      if (dx_w * dx_w + dy_w * dy_w <=
          sp->wound_radius * sp->wound_radius) {
        auto* stratum_grid = rm->GetDiffusionGrid(fields::kStratum);
        // Use raw voxel value (not interpolated) to prevent accumulation
        size_t vidx = stratum_grid->GetBoxIndex(qpos);
        real_t current_val = stratum_grid->GetConcentration(vidx);
        real_t new_val = static_cast<real_t>(new_stratum);
        if (new_val > current_val) {
          stratum_grid->ChangeConcentrationBy(vidx, new_val - current_val);
        }
      }
    }

    // Per-cell UWYN handoff: only cornified cells dissolve back to
    // continuum individually. They're dead keratinized cells at the
    // top of the stack -- nothing needs to build above them. All
    // lower strata (basale, spinosum, granulosum) persist as
    // mechanical scaffolding for vertical stratification.
    if (sp->handoff_delay > 0 &&
        !cell->IsStem() &&
        new_stratum == kCornified &&
        cell->GetStratumSteps() > sp->handoff_delay) {
      WriteScarValue(cell);
      cell->RemoveFromSimulation();
      return;
    }
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // DIFFERENTIATION_H_
