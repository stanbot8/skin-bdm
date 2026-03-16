#ifndef VOLUME_H_
#define VOLUME_H_

#include <cmath>
#include <limits>

#include "tissue/keratinocyte.h"

namespace bdm {
namespace skibidy {

// Properties attached to each epidermal volume (stratum)
struct VolumeProperties {
  Stratum stratum;
  real_t z_min;
  real_t z_max;
  bool agents_enabled = true;  // LOD: false = continuum only
};

// Singleton: maps z-positions to strata (epidermal and dermal sub-layers).
class VolumeManager {
 public:
  static VolumeManager* Get() {
    static VolumeManager instance;
    return &instance;
  }

  void Init(real_t z_spinous, real_t z_granular, real_t z_cornified,
            real_t z_papillary = -2.0, real_t z_reticular = -8.0) {
    volumes_[kBasal] = {kBasal, 0.0, z_spinous, true};
    volumes_[kSpinous] = {kSpinous, z_spinous, z_granular, false};
    volumes_[kGranular] = {kGranular, z_granular, z_cornified, false};
    volumes_[kCornified] = {kCornified, z_cornified,
                            std::numeric_limits<real_t>::infinity(), false};
    volumes_[kDermis] = {kDermis,
                         -std::numeric_limits<real_t>::infinity(), 0.0, false};
    volumes_[kPapillary] = {kPapillary, z_papillary, 0.0, false};
    volumes_[kReticular] = {kReticular, z_reticular, z_papillary, false};
    volumes_[kHypodermis] = {kHypodermis,
                             -std::numeric_limits<real_t>::infinity(),
                             z_reticular, false};
    initialized_ = true;
  }

  Stratum GetStratumAt(const Real3& pos) const {
    real_t z = pos[2];
    if (z >= 0) {
      if (z < volumes_[kBasal].z_max) return kBasal;
      if (z < volumes_[kSpinous].z_max) return kSpinous;
      if (z < volumes_[kGranular].z_max) return kGranular;
      return kCornified;
    }
    if (z >= volumes_[kPapillary].z_min) return kPapillary;
    if (z >= volumes_[kReticular].z_min) return kReticular;
    return kHypodermis;
  }

  const VolumeProperties& GetVolumeAt(const Real3& pos) const {
    return volumes_[static_cast<int>(GetStratumAt(pos))];
  }

  void SetAgentsEnabled(Stratum s, bool enabled) {
    volumes_[static_cast<int>(s)].agents_enabled = enabled;
  }

  bool AreAgentsEnabled(Stratum s) const {
    return volumes_[static_cast<int>(s)].agents_enabled;
  }

  const VolumeProperties& GetVolume(Stratum s) const {
    return volumes_[static_cast<int>(s)];
  }

  bool IsInitialized() const { return initialized_; }

 private:
  VolumeManager() = default;
  VolumeProperties volumes_[13] = {};
  bool initialized_ = false;
};

}  // namespace skibidy
}  // namespace bdm

#endif  // VOLUME_H_
