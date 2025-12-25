#ifndef CELL_CELL_FORCE_H_
#define CELL_CELL_FORCE_H_

#include "biodynamo.h"
#include "core/interaction_force.h"
#include "infra/sim_param.h"

namespace bdm {
namespace skibidy {

// Custom cell-cell force with configurable adhesion (attraction) and
// repulsion coefficients.  Keratinocytes are held together by desmosomes
// and E-cadherin; the attraction term models this intercellular adhesion.
class CellCellForce : public InteractionForce {
 public:
  CellCellForce() {}
  ~CellCellForce() override {}

  Real4 Calculate(const Agent* lhs, const Agent* rhs) const override {
    auto* sp = Simulation::GetActive()->GetParam()->Get<SimParam>();
    auto* c1 = dynamic_cast<const Cell*>(lhs);
    auto* c2 = dynamic_cast<const Cell*>(rhs);
    if (!c1 || !c2) return {0, 0, 0, 0};

    Real3 diff = c2->GetPosition() - c1->GetPosition();
    real_t dist = std::sqrt(diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2]);
    if (dist < 1e-10) return {0, 0, 0, 0};

    real_t r1 = c1->GetDiameter() / 2.0;
    real_t r2 = c2->GetDiameter() / 2.0;
    real_t overlap = r1 + r2 - dist;

    Real3 dir = {diff[0]/dist, diff[1]/dist, diff[2]/dist};
    real_t force_mag = 0;

    if (overlap > 0) {
      // Hertz-like repulsion (proportional to overlap^1.5)
      force_mag = -sp->repulsion_coeff * std::pow(overlap, 1.5);
    } else {
      // Adhesion (desmosomes/E-cadherin) -- short-range attraction
      real_t gap = -overlap;
      real_t max_adhesion_range = 2.0;
      if (gap < max_adhesion_range) {
        force_mag = sp->attraction_coeff * (1.0 - gap / max_adhesion_range);
      }
    }

    return {force_mag * dir[0], force_mag * dir[1], force_mag * dir[2], 0};
  }

  InteractionForce* NewCopy() const override {
    return new CellCellForce(*this);
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // CELL_CELL_FORCE_H_
