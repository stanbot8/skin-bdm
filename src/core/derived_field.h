#ifndef DERIVED_FIELD_H_
#define DERIVED_FIELD_H_

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "biodynamo.h"

namespace bdm {
namespace skibidy {

// Lightweight read-only derived field. Not a DiffusionGrid: no diffusion,
// no decay, no solver. A flat voxel array recomputed from source grids once
// per step, queryable by agents via Real3 position.
//
// Mirrors DiffusionGrid's spatial query interface (GetValue, GetConcentration,
// GetBoxIndex) so agents can treat derived fields like real PDE fields.
class DerivedField {
 public:
  DerivedField(const std::string& name, size_t resolution,
               real_t lo, real_t box_len)
      : name_(name), resolution_(resolution), lo_(lo), box_len_(box_len) {
    data_.resize(resolution * resolution * resolution, 0.0);
  }

  const std::string& GetName() const { return name_; }
  size_t GetResolution() const { return resolution_; }
  size_t GetNumBoxes() const { return data_.size(); }

  real_t GetConcentration(size_t idx) const { return data_[idx]; }
  void SetConcentration(size_t idx, real_t val) { data_[idx] = val; }

  // Coordinate query (mirrors DiffusionGrid::GetValue).
  // Nearest-voxel lookup matches BDM's DiffusionGrid behavior.
  real_t GetValue(const Real3& pos) const {
    return data_[GetBoxIndex(pos)];
  }

  size_t GetBoxIndex(const Real3& pos) const {
    auto clamp = [&](real_t v) -> size_t {
      int i = static_cast<int>((v - lo_) / box_len_);
      if (i < 0) return 0;
      if (static_cast<size_t>(i) >= resolution_) return resolution_ - 1;
      return static_cast<size_t>(i);
    };
    size_t ix = clamp(pos[0]);
    size_t iy = clamp(pos[1]);
    size_t iz = clamp(pos[2]);
    return ix + iy * resolution_ + iz * resolution_ * resolution_;
  }

  void Zero() {
    std::fill(data_.begin(), data_.end(), 0.0);
  }

 private:
  std::string name_;
  size_t resolution_;
  real_t lo_;
  real_t box_len_;
  std::vector<real_t> data_;
};

}  // namespace skibidy
}  // namespace bdm

#endif  // DERIVED_FIELD_H_
