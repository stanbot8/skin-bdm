#ifndef IMMUNE_CELL_H_
#define IMMUNE_CELL_H_

#include "biodynamo.h"

namespace bdm {
namespace skibidy {

// Immune cell types
enum ImmuneCellType { kNeutrophil = 0, kMacrophage = 1 };

// Immune cell state (lifecycle phase)
enum ImmuneState { kM1Active = 0, kM2Resolving = 1 };

// ---------------------------------------------------------------------------
// ImmuneCell — innate immune cell (neutrophil or macrophage) that produces
// or resolves inflammation in the wound bed.
// ---------------------------------------------------------------------------
class ImmuneCell : public Cell {
  BDM_AGENT_HEADER(ImmuneCell, Cell, 1);

 public:
  ImmuneCell() {}
  explicit ImmuneCell(const Real3& position) : Base(position) {}
  virtual ~ImmuneCell() {}

  // Type
  void SetImmuneCellType(ImmuneCellType t) {
    immune_type_ = static_cast<int>(t);
  }
  ImmuneCellType GetImmuneCellType() const {
    return static_cast<ImmuneCellType>(immune_type_);
  }

  // Age (steps since spawn)
  void SetAge(int a) { age_ = a; }
  int GetAge() const { return age_; }
  void IncrementAge() { age_++; }

  // State (M1/active vs M2/resolving)
  void SetState(ImmuneState s) {
    if (static_cast<int>(s) != state_) {
      state_age_ = 0;
    }
    state_ = static_cast<int>(s);
  }
  ImmuneState GetState() const { return static_cast<ImmuneState>(state_); }

  // Time spent in current state (resets on transition)
  int GetStateAge() const { return state_age_; }
  void IncrementStateAge() { state_age_++; }

 private:
  int immune_type_ = kNeutrophil;
  int age_ = 0;
  int state_ = kM1Active;
  int state_age_ = 0;
};

}  // namespace skibidy
}  // namespace bdm

#endif  // IMMUNE_CELL_H_
