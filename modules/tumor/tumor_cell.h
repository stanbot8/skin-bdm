#ifndef TUMOR_CELL_H_
#define TUMOR_CELL_H_

#include "biodynamo.h"
#include "tissue/keratinocyte.h"  // CellCyclePhase enum

namespace bdm {
namespace skibidy {

// ---------------------------------------------------------------------------
// TumorCell -- neoplastic cell with unlimited proliferation, faster cell
// cycle, and reduced contact inhibition.  Self-contained: composes with
// wound/immune events through shared fields (O2, calcium) but has no
// direct dependency on them.
// ---------------------------------------------------------------------------
class TumorCell : public Cell {
  BDM_AGENT_HEADER(TumorCell, Cell, 1);

 public:
  TumorCell() {}
  explicit TumorCell(const Real3& position) : Base(position) {}
  virtual ~TumorCell() {}

  void Initialize(const NewAgentEvent& event) override {
    Base::Initialize(event);
    if (auto* mother = dynamic_cast<TumorCell*>(event.existing_agent)) {
      // Daughter starts fresh cycle, inherits type
      cycle_phase_ = kG1;
      phase_elapsed_ = 0;
      cell_age_ = 0;
      g0_steps_ = 0;
      tumor_type_ = mother->tumor_type_;
    }
  }

  // Cell cycle phase (reuses CellCyclePhase enum from keratinocyte.h)
  void SetCyclePhase(CellCyclePhase p) { cycle_phase_ = static_cast<int>(p); }
  CellCyclePhase GetCyclePhase() const {
    return static_cast<CellCyclePhase>(cycle_phase_);
  }

  // Time elapsed in current cycle phase (hours)
  void SetPhaseElapsed(real_t t) { phase_elapsed_ = t; }
  real_t GetPhaseElapsed() const { return phase_elapsed_; }

  // Total cell age (hours)
  void SetCellAge(real_t a) { cell_age_ = a; }
  real_t GetCellAge() const { return cell_age_; }

  // Steps spent in G0 (quiescent); used for UWYN handoff
  void SetG0Steps(int s) { g0_steps_ = s; }
  int GetG0Steps() const { return g0_steps_; }
  void IncrementG0Steps() { g0_steps_++; }

  // Visualization attribute (exported to ParaView)
  void SetTumorType(int t) { tumor_type_ = t; }
  int GetTumorType() const { return tumor_type_; }

 private:
  int cycle_phase_ = kG1;
  real_t phase_elapsed_ = 0;
  real_t cell_age_ = 0;
  int g0_steps_ = 0;
  int tumor_type_ = 0;
};

}  // namespace skibidy
}  // namespace bdm

#endif  // TUMOR_CELL_H_
