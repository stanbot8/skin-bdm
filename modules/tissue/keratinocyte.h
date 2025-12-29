#ifndef KERATINOCYTE_H_
#define KERATINOCYTE_H_

#include "biodynamo.h"
#include "core/field_names.h"
#include "infra/sim_param.h"
#include "infra/util.h"

namespace bdm {
namespace skibidy {

// Epidermal stratum (mapped to color in visualization)
enum Stratum {
  kBasal = 0, kSpinous = 1, kGranular = 2, kCornified = 3,
  kDermis = 4,        // generic dermis (backward compat)
  // Scar encoding uses stratum+5 (values 5-8), so dermal sub-layers
  // start at 10 to avoid collision.
  kPapillary = 10,    // upper dermis: capillary loops, thin fibers
  kReticular = 11,    // lower dermis: thick collagen/elastin, appendages
  kHypodermis = 12    // subcutaneous fat, large vessels
};

// Cell cycle phases (G0 = quiescent, post-mitotic)
enum CellCyclePhase { kG1 = 0, kS = 1, kG2 = 2, kM = 3, kG0 = 4 };

// ---------------------------------------------------------------------------
// Keratinocyte — skin cell with stratum, cell cycle, and division state
// ---------------------------------------------------------------------------
class Keratinocyte : public Cell {
  BDM_AGENT_HEADER(Keratinocyte, Cell, 1);

 public:
  Keratinocyte() {}
  explicit Keratinocyte(const Real3& position) : Base(position) {}
  virtual ~Keratinocyte() {}

  void Initialize(const NewAgentEvent& event) override {
    Base::Initialize(event);
    if (auto* mother = dynamic_cast<Keratinocyte*>(event.existing_agent)) {
      stratum_ = mother->stratum_;
      ta_divisions_max_ = mother->ta_divisions_max_;

      // Asymmetric division: daughters of stem cells become TA cells;
      // daughters of TA cells inherit decremented count
      if (mother->divisions_left_ < 0) {
        divisions_left_ = mother->ta_divisions_max_;
      } else {
        divisions_left_ = mother->divisions_left_ - 1;
      }

      // Daughter starts fresh cell cycle
      cycle_phase_ = kG1;
      phase_elapsed_ = 0;
      cell_age_ = 0;
      stratum_steps_ = 0;
    }
  }

  // Stratum
  void SetStratum(Stratum s) {
    if (static_cast<int>(s) != stratum_) {
      stratum_steps_ = 0;
    }
    stratum_ = static_cast<int>(s);
  }
  Stratum GetStratum() const { return static_cast<Stratum>(stratum_); }

  // Stem / TA division tracking
  // -1 = stem cell (divides forever), >=0 = TA cell (divisions remaining)
  void SetDivisionsLeft(int d) { divisions_left_ = d; }
  int GetDivisionsLeft() const { return divisions_left_; }
  bool IsStem() const { return divisions_left_ < 0; }
  void SetTADivisionsMax(int d) { ta_divisions_max_ = d; }

  // Cell cycle phase
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

  // Steps spent in current stratum (for residence time tracking)
  void SetStratumSteps(int s) { stratum_steps_ = s; }
  int GetStratumSteps() const { return stratum_steps_; }
  void IncrementStratumSteps() { stratum_steps_++; }

  // Field-driven initialization: read local continuum state and derive
  // stem identity, TA division budget, and cell cycle readiness.
  // Called at spawn time so agents bootstrap from the continuum ground truth.
  void InitFromFields(Simulation* sim, const SimParam* sp) {
    auto* rm = sim->GetResourceManager();
    Real3 pos = ClampToBounds(GetPosition(), sim->GetParam());

    // Read local field values
    auto* ca_grid = rm->GetDiffusionGrid(fields::kCalcium);
    auto* kgf_grid = rm->GetDiffusionGrid(fields::kKGF);
    auto* o2_grid = rm->GetDiffusionGrid(fields::kOxygen);
    real_t ca = ca_grid->GetConcentration(ca_grid->GetBoxIndex(pos));
    real_t kgf = kgf_grid->GetConcentration(kgf_grid->GetBoxIndex(pos));
    real_t o2 = o2_grid->GetConcentration(o2_grid->GetBoxIndex(pos));

    // 1. Stem/TA decision: low calcium favors stem identity
    real_t ca_norm = ca / sp->ca_spinous_threshold;
    real_t stem_boost = std::max(0.0, 1.0 - ca_norm);
    real_t stem_prob = sp->stem_fraction * (1.0 + stem_boost);
    stem_prob = std::min(stem_prob, 1.0);

    auto* random = sim->GetRandom();
    if (random->Uniform(0, 1) < stem_prob) {
      SetDivisionsLeft(-1);
    } else {
      // 2. TA division budget: hypoxia reduces max divisions
      real_t o2_factor = std::min(1.0, o2 / sp->oxygen_prolif_threshold);
      int ta_divs = std::max(1,
          static_cast<int>(sp->max_ta_divisions * o2_factor));
      SetDivisionsLeft(ta_divs);
    }

    // 3. Cell cycle pre-advance: KGF readiness (Michaelis-Menten)
    real_t kgf_response = kgf / (sp->kgf_half_maximal + kgf);
    real_t g1_advance = kgf_response * sp->g1_duration * 0.5;
    SetPhaseElapsed(g1_advance);
  }

 private:
  int stratum_ = kBasal;
  int divisions_left_ = -1;   // -1 = stem cell
  int ta_divisions_max_ = 4;
  int cycle_phase_ = kG1;
  real_t phase_elapsed_ = 0;  // hours in current phase
  real_t cell_age_ = 0;       // total age in hours
  int stratum_steps_ = 0;     // steps in current stratum
};

}  // namespace skibidy
}  // namespace bdm

#endif  // KERATINOCYTE_H_
