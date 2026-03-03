#ifndef RA_PARAMS_H_
#define RA_PARAMS_H_

#include "biodynamo.h"

namespace bdm {
namespace skibidy {

struct RaParams {

  // Rheumatoid arthritis module (synovial inflammation and cartilage erosion)
  // Smolen et al. 2018 (doi:10.1038/nrdp.2018.1)
  bool enabled = false;
  // Autoimmune cytokine production (FLS, Th1/Th17 T cells)
  real_t autoimmune_source = 0.001;  // chronic TNF-alpha production
  real_t il6_autoimmune_source = 0.0015;  // chronic IL-6 production
  // Immune cell cytokine production (inflammation-proportional)
  real_t tnf_m1_rate = 0.001;  // TNF-alpha from M1 macrophages (scaled by inflammation)
  real_t tnf_neutrophil_rate = 0.0005;  // TNF-alpha from neutrophils (scaled by inflammation)
  real_t il6_m1_rate = 0.001;  // IL-6 from M1 macrophages (scaled by inflammation)
  // Autoimmune flare dynamics (sigmoid activation after wounding)
  real_t flare_delay = 200;  // steps before autoimmune flare onset (~0.8 days)
  real_t flare_steepness = 0.03;  // sigmoid steepness (higher = sharper onset)
  // TNF-alpha downstream effects
  real_t tnf_inflammation_coupling = 0.003;  // TNF -> NF-kB -> generic inflammation amplification
  real_t tnf_mmp_boost = 0.003;  // TNF -> MMP-1/3/13 transcription (Mengshol 2000)
  real_t tnf_vegf_boost = 0.002;  // TNF -> VEGF (pannus neovascularization)
  real_t tnf_il6_induction = 0.001;  // TNF -> NF-kB -> IL-6 transcription
  // IL-6 downstream effects
  real_t il6_inflammation_coupling = 0.002;  // IL-6 -> JAK/STAT3 -> inflammation amplification
  real_t il6_cartilage_boost = 0.00015;  // IL-6 -> RANKL -> osteoclast cartilage erosion
  // Cartilage integrity (structural tissue target)
  real_t cartilage_basal = 1.0;  // normalized healthy cartilage integrity
  real_t cartilage_mmp_degradation = 0.0005;  // aggrecanase/collagenase erosion
  real_t cartilage_tnf_degradation = 0.0002;  // direct TNF osteoclast activation
  // Synovial pathology
  real_t pannus_fibroblast_boost = 2.0;  // FLS hyperplasia multiplier on fibroblast spawn
  real_t m1_prolongation = 2.0;  // M1 duration multiplier (sustained Th1 polarization)
  // Synovial fluid / pannus tissue field
  real_t synovial_basal = 0.1;  // baseline synovial lining density (thin in healthy)
  real_t synovial_growth_rate = 0.001;  // pannus growth per unit TNF+IL-6
  real_t synovial_carrying_capacity = 1.0;  // max pannus density
  real_t synovial_cytokine_boost = 0.3;  // pannus amplifies local cytokine production
  real_t synovial_erosion_boost = 0.2;  // pannus proximity amplifies cartilage erosion
  // T cell density field (adaptive immunity proxy)
  // McInnes & Schett 2011 (doi:10.1056/NEJMra1004965): T cell infiltration
  // drives autoimmune cytokine production; Th1/Th17 polarization sustains
  // TNF-alpha and IL-17 in the synovium.
  real_t tcell_diffusion = 0.01;  // slow T cell migration in tissue
  real_t tcell_decay = 0.008;  // apoptosis and lymphatic egress
  real_t tcell_recruitment = 0.001;  // inflammation-driven extravasation
  real_t tcell_proliferation = 0.001;  // antigen-driven clonal expansion (flare-gated)
  real_t tcell_carrying_capacity = 1.0;  // local density cap (normalized)
  real_t tcell_cytokine_weight = 0.3;  // T cell contribution to autoimmune cytokine modulation
  // Subchondral bone integrity field
  // McInnes & Schett 2011: bone erosion via RANKL/osteoclast activation
  // is distinct from cartilage proteolysis and driven primarily by
  // TNF-alpha and IL-6 induced RANKL expression.
  real_t bone_basal = 1.0;  // normalized healthy bone integrity
  real_t bone_rankl_erosion = 0.00008;  // IL-6 -> RANKL -> osteoclast resorption
  real_t bone_tnf_erosion = 0.00004;  // TNF direct osteoclast activation
  real_t bone_depth_z = -6.0;  // z threshold for bone tissue (below reticular dermis)
  // Drug pharmacokinetics (treatment onset ramp)
  // Models delayed onset of biologic/DMARD action: linear ramp from
  // treatment start to full efficacy over drug_onset steps.
  int treatment_delay = 0;  // steps after wound before treatment starts (TOML: hours)
  int drug_onset = 0;  // steps for drug to reach peak effect (TOML: hours)
  // Treatment: anti-TNF biologic (infliximab/adalimumab/etanercept)
  real_t anti_tnf_clearance = 0.0;  // 0 = no treatment; >0 = TNF neutralization rate
  // Treatment: anti-IL-6R (tocilizumab)
  real_t anti_il6r_clearance = 0.0;  // 0 = no treatment; >0 = IL-6 neutralization rate
};

}  // namespace skibidy
}  // namespace bdm

#endif  // RA_PARAMS_H_
