#ifndef SIM_PARAM_H_
#define SIM_PARAM_H_

#include "biodynamo.h"
#include "infra/toml_compat.h"

// Module param structs
#include "angiogenesis/params.h"
#include "bioelectric/params.h"
#include "biofilm/params.h"
#include "blood/params.h"
#include "burn/params.h"
#include "dermis/params.h"
#include "diabetic/params.h"
#include "elastin/params.h"
#include "fibroblast/params.h"
#include "fibronectin/params.h"
#include "glucose/params.h"
#include "hemostasis/params.h"
#include "hyaluronan/params.h"
#include "immune/params.h"
#include "inflammation/params.h"
#include "lactate/params.h"
#include "lymphatic/params.h"
#include "mechanotransduction/params.h"
#include "mmp/params.h"
#include "neuropathy/params.h"
#include "nitric_oxide/params.h"
#include "perfusion/params.h"
#include "ph/params.h"
#include "photon/params.h"
#include "pressure/params.h"
#include "rheumatoid/params.h"
#include "ros/params.h"
#include "scab/params.h"
#include "scar/params.h"
#include "senescence/params.h"
#include "temperature/params.h"
#include "tumor/params.h"
#include "wound/params.h"

namespace bdm {
namespace skibidy {

struct SimParam : public ParamGroup {
  BDM_PARAM_GROUP_HEADER(SimParam, 1);

  // Cell cycle durations (hours)
  real_t g1_duration = 7.0;
  real_t s_duration = 6.0;
  real_t g2_duration = 3.0;
  real_t m_duration = 1.0;

  // Division mechanics
  real_t growth_rate = 5;
  int max_neighbors = 14;
  real_t lateral_scatter = 0.3;
  int max_ta_divisions = 4;
  real_t division_diameter = 5.0;
  real_t p_asymmetric = 0.7;     // probability of asymmetric stem division
  real_t stem_fraction = 0.50;   // fraction of seeded basal cells that are stem

  // Calcium gradient — prescribed sigmoid profile (Grabe & Neuber 2005)
  // The extracellular Ca2+ gradient is maintained by cellular transport
  // mechanisms; here it is modeled as a quasi-static environmental signal.
  real_t calcium_diffusion = 0.0;        // static prescribed field (no diffusion)
  real_t calcium_decay = 0.0;           // no decay: gradient persists
  real_t calcium_basal = 0.05;         // mM, concentration at basement membrane
  real_t calcium_peak = 1.5;           // mM, concentration at tissue surface
  real_t calcium_midpoint_z = 18.0;    // z-height where sigmoid reaches 50%
  real_t calcium_steepness = 2.0;      // sigmoid steepness (smaller = sharper)

  // Differentiation thresholds (dual: calcium + height)
  real_t ca_spinous_threshold = 0.1;   // mM, Ca2+ triggers spinous transition
  real_t ca_granular_threshold = 0.5;  // mM, Ca2+ triggers granular transition
  real_t ca_cornified_threshold = 1.0; // mM, Ca2+ triggers cornified transition
  real_t spinous_threshold = 6.0;      // z-height ceiling for calcium-driven basal identity

  // KGF (dermal growth factor) — models dermis as continuum source
  real_t kgf_diffusion = 0.0;      // static prescribed field (no diffusion)
  real_t kgf_decay = 0.0;          // no decay (persistent gradient)
  real_t kgf_basal_conc = 2.0;     // nM, concentration at basement membrane (z=0)
  real_t kgf_half_maximal = 0.5;   // nM, Michaelis-Menten Km (50% response)
  real_t kgf_max_boost = 1.0;      // max fold-increase in G1->S transition rate
  real_t kgf_decay_length = 5.0;   // z-scale for exponential decay of initial profile

  // Layer thicknesses (config values, relative to basement membrane at z=0)
  real_t basal_thickness = 3.0;       // basale layer height
  real_t spinous_thickness = 15.0;    // spinosum layer height
  real_t granular_thickness = 7.0;    // granulosum layer height
  real_t papillary_thickness = 2.0;   // papillary dermis depth
  real_t reticular_thickness = 6.0;   // reticular dermis depth

  // Computed z-boundaries (derived from thicknesses, not in config)
  real_t volume_z_spinous = 3.0;    // = basal_thickness
  real_t volume_z_granular = 18.0;  // = basal + spinous
  real_t volume_z_cornified = 25.0; // = basal + spinous + granular
  real_t dermal_z_papillary = -2.0; // = -papillary
  real_t dermal_z_reticular = -8.0; // = -(papillary + reticular)

  // Oxygen — dermal vasculature supplies O2 at basement membrane
  real_t oxygen_diffusion = 0.1;        // diffusion coefficient
  real_t oxygen_decay = 0.01;          // consumption/decay rate
  real_t oxygen_basal_conc = 1.0;      // normalized, 1.0 = arterial pO2 at z=0
  real_t oxygen_decay_length = 8.0;    // z-scale for initial exponential profile

  // Water / tissue moisture — hydration gradient from dermis to surface
  real_t water_diffusion = 1e-4;          // lateral diffusion (faster than O2)
  real_t water_decay = 0.002;            // background TEWL
  real_t water_basal_conc = 1.0;         // normalized, 1.0 = fully hydrated
  real_t water_decay_length = 12.0;      // z-scale for initial exponential profile
  real_t water_recovery_rate = 0.02;     // serum hydration rate per step
  real_t water_surface_loss_rate = 0.03; // evaporation at exposed wound surface per step
  real_t water_migration_threshold = 0.3;  // min moisture for full migration speed
  real_t water_prolif_threshold = 0.4;     // min moisture for full proliferation rate

  // Cell-cell mechanics (Hertz contact + adhesion)
  real_t repulsion_coeff = 5.0;     // steric repulsion strength
  real_t attraction_coeff = 0.0;    // adhesion disabled (desmosomes / E-cadherin)

  // Shedding & apoptosis (TOML: hours -> steps internally)
  int shedding_delay = 99999;   // 9999.9h = disabled
  int apoptosis_delay = 99999;  // 9999.9h = disabled

  // ---- Geometry (derived from [geometry] TOML section) ----
  // Physical dimensions in micrometers, converted to simulation units.
  real_t geo_patch_um = 2000.0;    // tissue lateral extent (um)
  real_t geo_depth_um = 1500.0;    // domain depth below basement membrane (um)
  real_t geo_height_um = 1500.0;   // domain height above basement membrane (um)
  real_t geo_margin_um = 500.0;    // empty margin beyond tissue (um)
  real_t geo_um_per_unit = 50.0;   // physical scale: 1 sim unit = N um
  real_t geo_voxel_um = 250.0;     // target voxel edge length (um)

  // Derived spatial bounds (computed from geometry, not directly in TOML)
  real_t tissue_min = -20;         // = -patch_um / (2 * um_per_unit)
  real_t tissue_max = 20;          // =  patch_um / (2 * um_per_unit)
  real_t domain_min = -30;         // = tissue_min - margin_um/um_per_unit
  real_t domain_max = 30;          // = tissue_max + margin_um/um_per_unit
  real_t domain_z_min = -30;       // = -depth_um / um_per_unit
  real_t domain_z_max = 30;        // =  height_um / um_per_unit

  // Per-axis simulation bounds for position clamping. BioDynaMo's Param
  // only has scalar bounds; these provide axis-aware clamping.
  Real3 bounds_min = {-20, -20, -30};
  Real3 bounds_max = {20, 20, 30};

  // Simulation (TOML: duration_days -> steps internally)
  int num_steps = 7200;         // 30 days = 7200 steps
  int grid_resolution = 12;     // DiffusionGrid boxes per axis (derived from geometry)
  real_t ref_box_length = 5.0;  // box_length (derived: voxel_um / um_per_unit)
  // wound_vascular_damage removed: superseded by VascularPDE
  int homeostatic_fold = 0;           // TOML: homeostatic_fold_h; steps in same stratum before agent folds to continuum (0=disabled)
  real_t dissolution_closure_pct = 90; // dissolve all agents once wound closure >= this % (0=disabled)

  // Per-layer angiogenesis rate multipliers
  real_t angio_papillary_factor = 1.5;   // capillaries regrow fastest
  real_t angio_reticular_factor = 1.0;   // baseline rate
  real_t angio_hypodermis_factor = 0.5;  // deep vessels recover slowly

  // Dynamic field coupling
  real_t calcium_recovery_rate = 0.002;  // fraction of target recovered per step
  real_t oxygen_prolif_threshold = 0.3;  // O2 below this suppresses proliferation
  bool oxygen_recovery_enabled = true;   // vasculature regenerates as wound heals

  // Cell migration -- multi-cue mechanistic chemotaxis/haptotaxis
  bool migration_enabled = true;
  real_t migration_speed = 3.0;               // intrinsic crawl force (calibrated: ~30 um/h)
  real_t migration_gradient_scale = 15.0;     // converts O2 gradient magnitude to speed [0.3,1]
  real_t cil_suppression = 1.0;               // speed multiplier when forward voxel is covered (CIL)
                                              // 1.0 = disabled; set < 1 when sheet-migration is added

  // Per-cell UWYN handoff (TOML: hours -> steps internally)
  int handoff_delay = 500;       // 50h before cell dissolves (0 = disabled)
  real_t m1_transition_threshold = 0.3;      // inflammation below this triggers M1->M2
  int m1_transition_min_age = 120;           // 12h min in M1 before cytokine-driven transition

  // Immune response -- shared
  real_t immune_cell_diameter = 3.0;         // neutrophil/macrophage diameter

  // Mechanistic replacements (test mode: false = parametric, true = mechanistic)
  bool mech_immune_recruitment = false;     // gradient-driven recruitment replaces spawn_rate/threshold/taper
  real_t mech_recruit_gradient_scale = 0.12; // gradient magnitude -> recruitment probability scaling
  real_t mech_recruit_saturation_k = 0.005;  // half-saturation for chemokine gradient (Michaelis-Menten)
  bool mech_m1_m2_transition = false;        // efferocytosis-count driven M1->M2 with timer ceiling
  int mech_efferocytosis_quota = 1;          // single engulfment triggers M2 program (PS receptor)
  real_t mech_m2_transition_rate = 0.25;     // probability per step of transition once quota met
  bool mech_collagen_deposition = false;     // TGF-beta receptor occupancy replaces flat rate
  real_t mech_collagen_tgfb_km = 0.035;      // TGF-beta half-max for collagen synthesis (Michaelis-Menten)
  real_t mech_collagen_vmax = 0.00025;         // TGF-b-responsive component Vmax
  real_t mech_collagen_basal = 0.00035;        // constitutive myofibroblast collagen rate (epigenetically locked)
  bool mech_vegf_production = false;         // HIF-1alpha stabilization replaces flat m2_vegf_rate
  real_t mech_hif_o2_threshold = 0.3;        // O2 below this stabilizes HIF-1alpha
  real_t mech_hif_vegf_rate = 0.003;         // max VEGF production at zero O2
  real_t anti_inflammation_diffusion = 0.05;   // anti-inflammatory diffusion coeff
  real_t anti_inflammation_decay = 0.002;      // anti-inflammatory decay rate
  real_t anti_inflammation_weight = 1.0;       // weight for net = pro - weight*anti

  // Efferocytosis-triggered M1-to-M2 (Extension 3)
  bool efferocytosis_enabled = false;          // proximity-based M1→M2
  real_t efferocytosis_radius = 5.0;           // search radius for dying neutrophils
  real_t efferocytosis_age_fraction = 0.8;     // neutrophil age fraction to be engulfable

  // Chemotaxis-driven immune migration (Extension 2)
  bool chemotaxis_enabled = false;             // follow inflammation gradient
  real_t chemotaxis_speed_scale = 1.0;         // gradient speed scaling factor
  real_t angio_vegf_rate = 0.01;               // perfusion recovery per unit VEGF per step
  real_t m2_vegf_rate = 0.001;                 // VEGF produced by M2 macrophages per step
  real_t dermal_fibroblast_depth = -3.0;    // z-coordinate for seeding (in dermis)
  real_t dermal_fibroblast_margin = 3.0;    // extra radius beyond wound for seeding ring
  real_t m2_tgfb_rate = 0.003;
  real_t decorin_sequestration_rate = 0.0;   // collagen-bound decorin neutralizes TGF-beta (Yamaguchi 1990)
  real_t matrikine_mmp_boost = 0.05;          // ECM fragment positive feedback on MMP production
  real_t fn_tgfb_max_boost = 2.0;            // max fold-increase in FN rate from TGF-beta (Ignotz & Massague 1986)
  real_t fn_tgfb_half_maximal = 0.01;        // TGF-beta Km for FN production (Michaelis-Menten half-max)
  real_t myofibroblast_fn_fraction = 0.3;    // EDA-FN production rate as fraction of activated rate (Serini et al. 1998)
  real_t age_decay = 0.0001;                   // AGE turnover (very slow; cross-linked proteins persist)
  real_t age_rage_inflammation = 0.001;        // RAGE-mediated NF-kB inflammation per AGE unit
  real_t age_collagen_crosslink = 0.15;        // fraction of MMP collagen degradation blocked by AGEs (0-1)
  real_t age_m1_prolongation = 0.05;           // AGE-RAGE extends M1 duration (multiplier on AGE level)
  real_t age_fibroblast_migration_impair = 0.1; // AGE-mediated fibroblast migration impairment
  real_t no_diffusion = 0.02;                  // fast diffusion (small molecule)
  real_t no_decay = 0.05;                      // rapid autoxidation
  real_t no_m1_production = 0.01;              // iNOS from M1 macrophages
  real_t no_neutrophil_production = 0.005;     // iNOS from neutrophils
  real_t no_vasodilation_factor = 0.05;         // perfusion boost from NO
  real_t no_antimicrobial_factor = 0.3;        // biofilm suppression factor
  real_t no_collagen_suppression = 0.03;       // anti-fibrotic effect
  real_t sasp_inflammation_rate = 0.001;          // SASP pro-inflammatory output
  real_t sasp_mmp_rate = 0.0005;                  // SASP MMP-3/9 output (as pro-MMP)
  real_t sasp_tgfb_rate = 0.0002;                 // SASP TGF-beta1 output
  real_t senolytic_clearance_rate = 0.0;          // treatment: senolytic drug clearance
  real_t nerve_regeneration_rate = 0.0002;        // regrowth rate
  real_t nerve_vegf_boost = 0.5;                  // VEGF promotes nerve sprouting
  real_t nerve_tgfb_inhibit = 0.3;               // TGF-beta inhibits neurite extension
  real_t substance_p_inflammation = 0.001;        // neurogenic inflammation
  real_t substance_p_proliferation = 0.15;        // keratinocyte proliferation boost
  real_t cgrp_vasodilation = 0.1;                 // CGRP perfusion boost
  real_t stiffness_yap_threshold = 0.3;        // stiffness above this promotes myofibroblast via YAP/TAZ
  real_t stiffness_contraction_rate = 0.0003;  // wound contraction from myofibroblast traction force
  real_t stiffness_scar_factor = 0.5;          // scar amplification from mechanical tension
  real_t edema_leak_rate = 0.003;              // vascular leak from inflammation
  real_t edema_drainage_rate = 0.01;           // lymphatic drainage of interstitial fluid
  real_t edema_o2_impairment = 0.3;            // edema reduces O2 delivery
  real_t edema_migration_impairment = 0.2;     // edema slows cell crawling
  real_t voltage_diffusion = 0.1;              // fast lateral spread (ionic current)
  real_t voltage_decay = 0.01;                 // charge leakage
  real_t voltage_epithelial_source = 0.05;     // Na+/K+ ATPase maintains TEP
  real_t galvanotaxis_strength = 0.4;          // voltage gradient biases fibroblast migration
  // TNF-alpha cytokine field
  real_t tnf_alpha_diffusion = 0.04;          // TNF-alpha diffusion (17 kDa homotrimer)
  real_t tnf_alpha_decay = 0.02;              // TNFR1/R2 receptor endocytosis clearance
  // IL-6 cytokine field (second inflammatory axis)
  // Kishimoto 2005 (doi:10.1146/annurev.immunol.23.021704.115806)
  real_t il6_diffusion = 0.05;                // IL-6 diffusion (21 kDa, slightly faster than TNF)
  real_t il6_decay = 0.015;                   // IL-6R/gp130 receptor endocytosis clearance

  // Derived composite fields (precomputed once per step from source grids)
  real_t ecm_weight_collagen = 0.4;       // collagen contribution to ECM quality
  real_t ecm_weight_fibronectin = 0.3;    // fibronectin contribution
  real_t ecm_weight_elastin = 0.2;        // elastin contribution
  real_t ecm_weight_fibrin = 0.1;         // fibrin provisional matrix contribution
  real_t ecm_migration_boost = 1.5;       // composite ECM migration speed boost

  // Multi-resolution: structural fields (D=0) at coarser grid
  // 0 = same as grid_resolution; nonzero = coarse resolution
  int grid_resolution_structural = 6;

  // Scar formation: dissolved agents stamp stratum+5 into Stratum field
  bool scar_enabled = true;

  // Metrics export (TOML: hours -> steps internally)
  int metrics_interval = 100;    // 10h between CSV rows (0 = disabled)
  bool metrics_autoopen = true;  // open CSV in default app after simulation
  bool headless = false;         // suppress all xdg-open/GUI calls (set true for batch/AI runs)
  real_t m1_cytokine_taper = 0.002;          // M1 macrophage cytokine exp(-k*state_age)
  real_t m2_tgfb_taper = 0.003;              // M2 macrophage TGF-beta exp(-k*state_age)

  // PDE sub-cycling: solve slow fields every N simulation steps instead of
  // every step. Reduces FTCS solver invocations by (N-1)/N. CFL-safe for
  // any D*subcycle*dt/dx^2 < 0.5 (verified at default parameters).
  // 1 = every step (no sub-cycling), 5 = every 5 steps, etc.
  int subcycle_slow = 1;           // water (D=1e-4), hyaluronan (D=0.01)
  int subcycle_medium = 1;         // inflammation (D=0.05), TGF-beta (D=0.03), MMP (D=0.02)

  // Agent behavior sub-cycling: reduce per-agent computation for cells in
  // steady-state homeostasis (far from wound) or with slowly-changing cues.
  // 1 = every step (no sub-cycling). Higher values skip expensive ops.
  int derived_fields_subcycle = 5;  // recompute ECM/viability composites every N steps
  int migration_subcycle = 3;       // recompute O2/VEGF/FN gradients every N steps
  int homeostasis_subcycle = 3;     // skip full behavior for distant cells every N steps

  // Basal density continuum: logistic PDE representing homeostatic basal layer
  // outside the wound zone. Agents that drift past demotion_radius_factor x
  // wound_radius fold back into this field -- reducing agent count to only the
  // active wound front. Cells are never cheaper than continuum.
  bool basal_density_enabled = true;
  real_t basal_density_max = 1.0;          // normalized carrying capacity
  real_t basal_density_div_rate = 0.01;    // per-step growth rate at full resources
  real_t basal_density_diff_rate = 0.005;  // per-step differentiation outflux
  int basal_density_subcycle = 5;          // run density PDE every N steps
  real_t demotion_radius_factor = 2.8;     // agents beyond factor*r demote to density

  // Hot-reload: check bdm.toml for parameter changes every metrics_interval
  // steps. Enables interactive tuning without restarting the simulation.
  bool hot_reload = false;

  // Debug profile: per-module diagnostic output (TOML: [skin.debug])
  bool debug_immune = false;       // immune cell writes, spawning, migration
  bool debug_fibroblast = false;   // collagen, TGF-beta, state transitions
  bool debug_wound = false;        // wound creation, resolution, coverage
  bool debug_scaled_grid = false;  // ScaledGrid init (agent_factor)
  bool debug_perf = false;         // wall-clock timing of major operations

  // Load all skin parameters from a toml++ table. Called after Simulation
  // construction and by hot-reload. Bypasses BDM's config chain so skibidy
  // works with any BDM version (v1.04 cpptoml or v1.05+ toml++).
  // Module parameter structs
  AngiogenesisParams angiogenesis;
  BioelectricParams bioelectric;
  BiofilmParams biofilm;
  BloodParams blood;
  BurnParams burn;
  DermisParams dermis;
  DiabeticParams diabetic;
  ElastinParams elastin;
  FibroblastParams fibroblast;
  FibronectinParams fibronectin;
  GlucoseParams glucose_mod;
  HemostasisParams hemostasis;
  HyaluronanParams hyaluronan;
  ImmuneParams immune;
  InflammationParams inflammation;
  LactateParams lactate;
  LymphaticParams lymphatic;
  MechanotransductionParams mechanotransduction;
  MmpParams mmp;
  NeuropathyParams neuropathy;
  NitricOxideParams nitric_oxide;
  PerfusionParams perfusion;
  PhParams ph;
  PhotonParams photon;
  PressureParams pressure;
  RaParams ra;
  RosParams ros;
  ScabParams scab;
  ScarParams scar;
  SenescenceParams senescence;
  TemperatureParams temperature;
  TumorParams tumor;
  WoundParams wound;

  void LoadConfig(const skibidy::TomlConfig& config);

  // Post-load sanity checks: abort early on obviously broken configs.
  void ValidateConfig() const;
};

}  // namespace skibidy
}  // namespace bdm

#endif  // SIM_PARAM_H_
