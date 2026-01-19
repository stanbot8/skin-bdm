#ifndef FIBROBLAST_H_
#define FIBROBLAST_H_

#include "biodynamo.h"

namespace bdm {
namespace skibidy {

// Fibroblast differentiation states
enum FibroblastState {
  kFibroQuiescent = 0,
  kFibroActivated = 1,
  kMyofibroblast = 2,
};

// ---------------------------------------------------------------------------
// Fibroblast -- dermal cell recruited to wound site that deposits collagen
// (scar tissue). Quiescent fibroblasts activate in response to TGF-beta,
// then differentiate into myofibroblasts that produce TGF-beta (positive
// feedback) and deposit collagen proportional to local TGF-beta.
// ---------------------------------------------------------------------------
class Fibroblast : public Cell {
  BDM_AGENT_HEADER(Fibroblast, Cell, 1);

 public:
  Fibroblast() {}
  explicit Fibroblast(const Real3& position) : Base(position) {}
  virtual ~Fibroblast() {}

  void SetFibroblastState(FibroblastState s) {
    if (static_cast<int>(s) != state_) {
      state_age_ = 0;
    }
    state_ = static_cast<int>(s);
  }
  FibroblastState GetFibroblastState() const {
    return static_cast<FibroblastState>(state_);
  }

  void SetAge(int a) { age_ = a; }
  int GetAge() const { return age_; }
  void IncrementAge() { age_++; }

  int GetStateAge() const { return state_age_; }
  void IncrementStateAge() { state_age_++; }

 private:
  int state_ = kFibroQuiescent;
  int age_ = 0;
  int state_age_ = 0;
};

}  // namespace skibidy
}  // namespace bdm

#endif  // FIBROBLAST_H_
