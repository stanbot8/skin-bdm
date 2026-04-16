#ifndef BIOFILM_POST_HOOK_H_
#define BIOFILM_POST_HOOK_H_

#include "core/grid_registry.h"
#include "core/voxel_snapshot.h"
#include "core/signal_board.h"

namespace bdm {
namespace skibidy {

// Post hook for biofilm module.
// Biofilm (polymicrobial community) seeds at a configurable delay post-wound
// and grows logistically on open wound surface. Produces pro-inflammatory
// cytokines (IL-1, TNF-alpha) that sustain M1 macrophage polarization.
// Growth modulated by pH, temperature (Q10), glucose, and NO.
// James et al. 2008 (doi:10.1111/j.1524-475X.2007.00321.x)
struct BiofilmPostHook {
  DiffusionGrid* biofilm_grid = nullptr;
  DiffusionGrid* infl_grid = nullptr;
  DiffusionGrid* temp_grid = nullptr;
  DiffusionGrid* glucose_grid = nullptr;
  DiffusionGrid* no_grid = nullptr;
  const SimParam* sp_ = nullptr;
  bool active = false;

  inline void Init(const GridRegistry& reg, SignalBoard& sig) {
    sp_ = reg.Params();
    active = sp_->biofilm.enabled;
    if (!active) return;
    biofilm_grid = reg.Get(fields::kBiofilmId);
    if (!biofilm_grid) { active = false; return; }
    infl_grid = reg.InflammationGrid();
    if (sp_->temperature.enabled)
      temp_grid = reg.Get(fields::kTemperatureId);
    if (sp_->glucose_mod.enabled)
      glucose_grid = reg.Get(fields::kGlucoseId);
    if (sp_->nitric_oxide.enabled)
      no_grid = reg.Get(fields::kNitricOxideId);
  }

  // Called per epidermal wound voxel.
  inline void ApplyEpiWound(const VoxelSnapshot& snap, SignalBoard& sig) {
    if (snap.stratum >= 1.0) return;

    size_t bf_idx = snap.coarse_si;
    real_t bf = biofilm_grid->GetConcentration(bf_idx);

    // Seeding: first colonization event
    bool seeding_step =
        (snap.wound_age == static_cast<uint64_t>(sp_->biofilm.seed_delay));
    if (seeding_step && bf < 1e-10) {
      biofilm_grid->ChangeConcentrationBy(bf_idx,
          sp_->biofilm.seed_amount * snap.coarse_w);
      bf += sp_->biofilm.seed_amount * snap.coarse_w;
    }

    if (bf <= 1e-10) return;

    // Logistic growth with environmental modulation
    real_t eff_growth = sp_->biofilm.growth_rate *
        (1.0 + sp_->ph.biofilm_boost * snap.ph);

    // Temperature Q10: warmer wound promotes bacterial growth
    if (temp_grid) {
      eff_growth *= Q10Factor(temp_grid->GetConcentration(snap.idx),
                                sp_->temperature.q10_biofilm);
    }
    // Glucose fuels bacterial metabolism
    if (glucose_grid) {
      real_t gluc = glucose_grid->GetConcentration(snap.idx);
      eff_growth *= (1.0 + sp_->glucose_mod.bacterial_consumption * gluc);
    }
    // NO antimicrobial: reactive nitrogen species suppress growth
    if (no_grid) {
      real_t no_val = no_grid->GetConcentration(snap.idx);
      eff_growth *= std::max(0.0,
          1.0 - sp_->no_antimicrobial_factor * no_val);
    }

    real_t delta = eff_growth * bf *
        (1.0 - bf / sp_->biofilm.carrying_capacity);
    if (delta > 1e-10) {
      biofilm_grid->ChangeConcentrationBy(bf_idx, delta * snap.coarse_w);
    }
    // Inflammation is fine grid: no correction
    real_t infl_delta = sp_->biofilm.inflammation_rate * bf;
    if (infl_delta > 1e-10 && infl_grid) {
      infl_grid->ChangeConcentrationBy(snap.idx, infl_delta);
    }
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // BIOFILM_POST_HOOK_H_
