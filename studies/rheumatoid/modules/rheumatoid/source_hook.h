#ifndef RHEUMATOID_SOURCE_HOOK_H_
#define RHEUMATOID_SOURCE_HOOK_H_

#include "core/hook_api.h"

namespace bdm {
namespace skibidy {

// Source hook for rheumatoid arthritis module.
// Five TNF-alpha/IL-6 production pathways:
//   1. Autoimmune: chronic FLS and Th1/Th17 T cell secretion (flare dynamics)
//   2. Immune cell: M1 macrophage and neutrophil production (inflammation proxy)
//   3. TNF-induced IL-6: TNF-alpha -> NF-kB -> IL-6 transcription
//   4. Pannus amplification: synovial hyperplasia boosts local cytokine output
//   5. T cell density modulation: local T cell infiltration scales autoimmune production
// Plus T cell recruitment/proliferation, pannus growth, drug PK ramp,
// and treatment clearance (anti-TNF, anti-IL-6R).
// Feldmann & Maini 2003 (doi:10.1038/nm939)
// Kishimoto 2005 (doi:10.1146/annurev.immunol.23.021704.115806)
// Firestein 2003 (doi:10.1038/nature01661)
// McInnes & Schett 2011 (doi:10.1056/NEJMra1004965)
struct RASourceHook {
  DiffusionGrid* tnf_grid = nullptr;
  DiffusionGrid* il6_grid = nullptr;
  DiffusionGrid* infl_grid = nullptr;
  DiffusionGrid* synovial_grid = nullptr;
  DiffusionGrid* tcell_grid = nullptr;
  const SimParam* sp_ = nullptr;
  bool active = false;
  real_t flare_factor = 0;  // sigmoid flare envelope [0..1]
  real_t drug_factor = 0;   // drug PK ramp [0..1] (0=no drug, 1=full effect)

  inline void Init(const GridRegistry& reg, SignalBoard& sig) {
    sp_ = reg.Params();
    active = sp_->ra.enabled;
    if (!active) return;

    tnf_grid = reg.Get(fields::kTNFAlphaId);
    il6_grid = reg.Get(fields::kIL6Id);
    if (!tnf_grid || !il6_grid) { active = false; return; }

    infl_grid = reg.InflammationGrid();
    synovial_grid = reg.Get(fields::kSynovialFluidId);
    tcell_grid = reg.Get(fields::kTCellDensityId);

    // Autoimmune flare dynamics: sigmoid activation after wounding.
    uint64_t step = reg.Step();
    uint64_t wound_step = static_cast<uint64_t>(sp_->wound.trigger_step);
    if (step > wound_step) {
      real_t wound_age = static_cast<real_t>(step - wound_step);
      real_t delay = sp_->ra.flare_delay;
      real_t k = sp_->ra.flare_steepness;
      flare_factor = 1.0 / (1.0 + std::exp(-k * (wound_age - delay)));
    }

    // Drug PK ramp: linear onset from treatment_delay to treatment_delay + drug_onset.
    bool has_treatment = sp_->ra.anti_tnf_clearance > 0 ||
                         sp_->ra.anti_il6r_clearance > 0;
    if (has_treatment && step > wound_step) {
      real_t wound_age = static_cast<real_t>(step - wound_step);
      real_t start = static_cast<real_t>(sp_->ra.treatment_delay);
      real_t onset = static_cast<real_t>(sp_->ra.drug_onset);
      if (wound_age >= start) {
        drug_factor = (onset > 0 && wound_age < start + onset)
                      ? (wound_age - start) / onset : 1.0;
      }
    } else if (has_treatment && wound_step == 0) {
      drug_factor = 1.0;
    }
  }

  // Dermal voxels: autoimmune + immune cell cytokine production
  // + pannus growth + T cell dynamics + drug clearance
  inline void ApplyDermal(const VoxelSnapshot& snap, SignalBoard& sig) {
    if (!snap.post_wound) return;

    // T cell density modulation: local T cell infiltration scales autoimmune output
    real_t tcell_mod = 1.0;
    if (tcell_grid) {
      real_t tc_val = tcell_grid->GetConcentration(snap.idx);
      tcell_mod = 1.0 + sp_->ra.tcell_cytokine_weight * tc_val;

      // T cell recruitment: inflammation-driven extravasation
      real_t infl_for_tc = 0;
      if (infl_grid) infl_for_tc = infl_grid->GetConcentration(snap.idx);
      real_t cap = sp_->ra.tcell_carrying_capacity;
      if (tc_val < cap && infl_for_tc > 1e-10) {
        real_t recruit = sp_->ra.tcell_recruitment * infl_for_tc
                         * (1.0 - tc_val / cap);
        tcell_grid->ChangeConcentrationBy(snap.idx, recruit);
      }

      // T cell proliferation: antigen-driven clonal expansion (flare-gated)
      if (tc_val > 1e-10 && tc_val < cap) {
        real_t prolif = sp_->ra.tcell_proliferation * tc_val
                        * (1.0 - tc_val / cap) * flare_factor;
        tcell_grid->ChangeConcentrationBy(snap.idx, prolif);
      }
    }

    // 1. Autoimmune source (FLS, T cells) with flare envelope and T cell modulation
    real_t tnf_auto = sp_->ra.autoimmune_source * flare_factor * tcell_mod;
    real_t il6_auto = sp_->ra.il6_autoimmune_source * flare_factor * tcell_mod;

    // 2. Immune cell production: inflammation as proxy for M1/neutrophil density
    real_t infl_val = 0;
    if (infl_grid) infl_val = infl_grid->GetConcentration(snap.idx);

    real_t tnf_immune = (sp_->ra.tnf_m1_rate + sp_->ra.tnf_neutrophil_rate) * infl_val;
    real_t il6_immune = sp_->ra.il6_m1_rate * infl_val;

    // 3. TNF-induced IL-6 transcription (NF-kB cascade)
    real_t tnf_cur = tnf_grid->GetConcentration(snap.idx);
    real_t il6_from_tnf = sp_->ra.tnf_il6_induction * tnf_cur;

    // 4. Pannus amplification: synovial hyperplasia boosts cytokine output
    real_t pannus_boost = 1.0;
    if (synovial_grid) {
      size_t syn_si = snap.coarse_si;
      real_t syn_val = synovial_grid->GetConcentration(syn_si);
      pannus_boost = 1.0 + sp_->ra.synovial_cytokine_boost * syn_val;

      // Pannus growth: logistic, driven by local TNF+IL-6, flare-gated
      real_t il6_cur = il6_grid->GetConcentration(snap.idx);
      real_t cytokine_drive = tnf_cur + il6_cur;
      real_t cap = sp_->ra.synovial_carrying_capacity;
      if (cytokine_drive > 1e-10 && syn_val < cap) {
        real_t growth = sp_->ra.synovial_growth_rate * cytokine_drive
                        * (1.0 - syn_val / cap) * flare_factor;
        synovial_grid->ChangeConcentrationBy(syn_si, growth * snap.coarse_w);
      }
    }

    // Apply TNF-alpha sources (scaled by pannus density)
    real_t total_tnf = (tnf_auto + tnf_immune) * pannus_boost;
    if (total_tnf > 1e-10) {
      tnf_grid->ChangeConcentrationBy(snap.idx, total_tnf);
    }

    // Apply IL-6 sources (scaled by pannus density)
    real_t total_il6 = (il6_auto + il6_immune + il6_from_tnf) * pannus_boost;
    if (total_il6 > 1e-10) {
      il6_grid->ChangeConcentrationBy(snap.idx, total_il6);
    }

    // Anti-TNF biologic treatment: neutralize TNF-alpha (PK-ramped)
    if (sp_->ra.anti_tnf_clearance > 0 && tnf_cur > 1e-10 && drug_factor > 0) {
      real_t clear = std::min(tnf_cur,
          sp_->ra.anti_tnf_clearance * tnf_cur * drug_factor);
      tnf_grid->ChangeConcentrationBy(snap.idx, -clear);
    }

    // Anti-IL-6R treatment (tocilizumab): block IL-6 signaling (PK-ramped)
    if (sp_->ra.anti_il6r_clearance > 0 && drug_factor > 0) {
      real_t il6_cur = il6_grid->GetConcentration(snap.idx);
      if (il6_cur > 1e-10) {
        real_t clear = std::min(il6_cur,
            sp_->ra.anti_il6r_clearance * il6_cur * drug_factor);
        il6_grid->ChangeConcentrationBy(snap.idx, -clear);
      }
    }
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // RHEUMATOID_SOURCE_HOOK_H_
