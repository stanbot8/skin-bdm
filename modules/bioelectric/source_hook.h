#ifndef BIOELECTRIC_SOURCE_HOOK_H_
#define BIOELECTRIC_SOURCE_HOOK_H_

#include "core/grid_registry.h"
#include "core/voxel_snapshot.h"
#include "core/signal_board.h"

namespace bdm {
namespace skibidy {

// Source hook for bioelectric module.
// Na+/K+ ATPase maintains transepithelial potential (TEP) in intact
// epithelium. Wound edge (stratum = 0 or -1) has collapsed TEP,
// creating lateral voltage gradient that drives galvanotaxis.
// Barker et al. 1982 (doi:10.1152/ajpregu.1982.242.3.R358)
struct BioelectricSourceHook {
  DiffusionGrid* volt_grid = nullptr;
  const SimParam* sp_ = nullptr;
  bool active = false;

  inline void Init(const GridRegistry& reg, SignalBoard& sig) {
    sp_ = reg.Params();
    active = sp_->bioelectric.enabled;
    if (!active) return;
    volt_grid = reg.Get(fields::kVoltageId);
    if (!volt_grid) active = false;
  }

  // Epidermal wound: TEP source from intact epithelium. Fully self-contained.
  inline void ApplyEpiWound(const VoxelSnapshot& snap, SignalBoard& sig) {
    real_t epi_gate = std::max(static_cast<real_t>(0),
                               std::min(static_cast<real_t>(1), snap.stratum));
    real_t source_rate = sp_->voltage_epithelial_source * epi_gate;
    if (sp_->diabetic.mode) {
      source_rate *= sp_->diabetic.voltage_factor;
    }
    if (source_rate > 1e-10) {
      volt_grid->ChangeConcentrationBy(snap.idx, source_rate);
    }
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // BIOELECTRIC_SOURCE_HOOK_H_
