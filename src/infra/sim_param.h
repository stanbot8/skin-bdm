#ifndef SIM_PARAM_H_
#define SIM_PARAM_H_

#include "biodynamo.h"
#include "infra/toml_compat.h"

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

  // Volume boundaries (z-heights for epidermal strata)
  real_t volume_z_spinous = 3.0;    // basal/spinous boundary
  real_t volume_z_granular = 18.0;  // spinous/granular boundary
  real_t volume_z_cornified = 25.0; // granular/cornified boundary

  // Dermal sub-layer boundaries (z < 0)
  real_t dermal_z_papillary = -2.0;   // papillary/reticular boundary
  real_t dermal_z_reticular = -8.0;   // reticular/hypodermis boundary

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

  // Tissue extent (x/y range for boundary removal)
  real_t tissue_min = 0;
  real_t tissue_max = 30;

  // Per-axis simulation bounds (override Param::min_bound/max_bound for
  // position clamping). BioDynaMo's Param only has scalar bounds; these
  // provide axis-aware clamping until upstream supports per-axis bounds.
  // Defaults match skibidy.h set_param: x/y=[0,40], z=[-10,40].
  Real3 bounds_min = {0, 0, -10};
  Real3 bounds_max = {40, 40, 40};

  // Simulation (TOML: duration_days -> steps internally)
  int num_steps = 7200;         // 30 days = 7200 steps
  int grid_resolution = 10;     // DiffusionGrid boxes per axis
  real_t ref_box_length = 5.0;  // box_length at calibration resolution (res=10, domain=50)

  // Wound event
  real_t wound_center_x = 15.0;    // x-center of punch biopsy
  real_t wound_center_y = 15.0;    // y-center of punch biopsy
  real_t wound_radius = 5.0;       // radius of circular wound
  int wound_trigger_step = 0;      // TOML: trigger_h; 0h = start wounded
  bool wound_enabled = true;       // master switch
  real_t wound_inward_bias = 0.3;    // 0=pure random, 1=pure center-directed
  real_t wound_vascular_damage = 0.5;  // (superseded by VascularPDE; kept for compat)

  // Vascular perfusion field
  real_t perfusion_diffusion = 0.03;    // capillary sprouting / lateral spread
  real_t perfusion_decay = 0.0;         // no decay (stable vessels)
  real_t perfusion_basal = 1.0;         // healthy dermis = fully perfused
  real_t perfusion_angio_rate = 0.005;  // recovery speed per step
  int perfusion_angio_delay = 240;      // TOML: 24h before sprouting begins

  // Per-layer perfusion baselines (fraction of perfusion_basal)
  // Braverman 2000 (doi:10.1046/j.1087-0024.2000.00010.x)
  real_t perfusion_papillary_fraction = 1.0;    // dense capillary loops
  real_t perfusion_reticular_fraction = 0.7;    // arterioles and venules
  real_t perfusion_hypodermis_fraction = 0.3;   // sparse large vessels

  // Per-layer angiogenesis rate multipliers
  real_t angio_papillary_factor = 1.5;   // capillaries regrow fastest
  real_t angio_reticular_factor = 1.0;   // baseline rate
  real_t angio_hypodermis_factor = 0.5;  // deep vessels recover slowly

  // Dynamic field coupling
  real_t calcium_recovery_rate = 0.002;  // fraction of target recovered per step
  real_t oxygen_prolif_threshold = 0.3;  // O2 below this suppresses proliferation
  bool oxygen_recovery_enabled = true;   // vasculature regenerates as wound heals

  // Cell migration
  bool migration_enabled = true;       // active crawling toward wound center
  real_t migration_speed = 2.0;        // tractor force magnitude

  // Per-cell UWYN handoff (TOML: hours -> steps internally)
  int handoff_delay = 500;       // 50h before cell dissolves (0 = disabled)

  // Inflammation field
  real_t inflammation_diffusion = 0.05;          // lateral cytokine spread
  real_t inflammation_decay = 0.003;             // background clearance
  real_t inflammation_migration_threshold = 1.0; // above this, migration suppressed
  real_t inflammation_prolif_threshold = 1.5;    // above this, proliferation suppressed
  real_t wound_inflammation_source_rate = 0.003; // DAMPs/IL-1a from open wound (per step, scaled by 1-stratum)
  real_t wound_inflammation_source_taper = 0.0;  // exp decay rate for DAMP clearance (per step; 0=no taper)

  // Immune response -- neutrophils (TOML: hours -> steps internally)
  int neutrophil_spawn_delay = 20;           // 2h (Kim et al. 2008)
  int neutrophil_spawn_waves = 3;            // recruitment waves (gradual rise)
  int neutrophil_spawn_window = 240;         // 24h between first and last wave
  int neutrophil_lifespan = 480;             // 48h hard ceiling (Wilgus et al. 2013)
  int neutrophil_min_survival = 120;         // 12h before apoptosis possible
  real_t neutrophil_apoptosis_rate = 0.0029; // half-life ~24h (Wilgus et al. 2013, doi:10.1089/wound.2012.0383)

  // Immune response -- macrophages (TOML: hours -> steps internally)
  int macrophage_spawn_delay = 240;          // 24h (Rodero & Khosrotehrani 2010)
  real_t macrophage_spawn_threshold = 0.1;   // min inflammation for monocyte extravasation
  real_t macrophage_spawn_rate = 0.1;        // per-step spawn probability per unit excess inflammation
  real_t macrophage_spawn_taper = 0.0;       // per-hour decay of recruitment probability (chemokine gradient decline)
  int macrophage_m1_duration = 480;          // 48h in M1 (Krzyszczyk et al. 2018)
  real_t m1_transition_threshold = 0.3;      // inflammation below this triggers M1->M2
  int m1_transition_min_age = 120;           // 12h min in M1 before cytokine-driven transition
  int macrophage_lifespan = 6720;            // 672h (28d) hard ceiling (Lucas et al. 2010)
  int macrophage_min_survival = 240;         // 24h before apoptosis possible
  real_t macrophage_apoptosis_rate = 0.00072; // half-life ~4d (Lucas et al. 2010, doi:10.4049/jimmunol.0903356)
  real_t macrophage_emigration_rate = 0.02;   // M2 exit rate once stratum re-epithelialized

  // Immune response -- shared
  real_t immune_cell_diameter = 3.0;         // neutrophil/macrophage diameter
  real_t immune_cytokine_rate = 0.002;       // inflammation added per cell per step
  real_t immune_resolution_rate = 0.005;     // inflammation consumed by M2 per step
  real_t immune_migration_speed = 1.5;       // tractor force magnitude

  // Split pro/anti-inflammatory fields (Extension 4)
  bool split_inflammation_enabled = false;     // use separate pro/anti fields
  real_t anti_inflammation_diffusion = 0.05;   // anti-inflammatory diffusion coeff
  real_t anti_inflammation_decay = 0.002;      // anti-inflammatory decay rate
  real_t anti_inflammation_weight = 1.0;       // weight for net = pro - weight*anti

  // Diabetic / chronic wound model (Extension 5)
  bool diabetic_mode = false;                  // impair M1→M2 transition
  real_t diabetic_m1_duration_factor = 3.0;    // M1 duration multiplier
  real_t diabetic_resolution_factor = 0.3;     // M2 resolution rate multiplier
  real_t diabetic_efferocytosis_factor = 0.5;  // efferocytosis probability multiplier
  real_t diabetic_prolif_factor = 0.5;               // keratinocyte G1->S multiplier
  real_t diabetic_migration_factor = 0.45;            // keratinocyte migration speed multiplier
  real_t diabetic_fibroblast_activation_factor = 2.0; // activation delay multiplier (>1 = slower)
  real_t diabetic_collagen_factor = 0.4;              // collagen deposition multiplier
  real_t diabetic_fibroblast_lifespan_factor = 0.6;   // fibroblast lifespan multiplier
  real_t diabetic_baseline_inflammation = 0.0005;     // per-step inflammation added to wound
  real_t diabetic_neutrophil_factor = 1.8;            // neutrophil spawn count multiplier
  real_t diabetic_neutrophil_lifespan_factor = 3.0;   // neutrophil lifespan multiplier
  real_t diabetic_neutrophil_apoptosis_factor = 0.3; // impaired neutrophil apoptosis
  real_t diabetic_neutrophil_waves_factor = 2.0;     // recruitment wave multiplier
  real_t diabetic_neutrophil_window_factor = 3.0;    // spawn window multiplier
  real_t diabetic_macrophage_taper_factor = 0.3;     // recruitment taper multiplier
  real_t diabetic_macrophage_emigration_factor = 0.5; // M2 emigration rate multiplier (impaired lymphatic clearance)
  real_t diabetic_inflammation_sensitivity = 2.0;     // keratinocyte inflammation sensitivity multiplier
  real_t diabetic_migration_infl_K = 0.08;            // Hill half-max for migration inflammation gating (vs 1.0 normal)
  real_t diabetic_prolif_infl_K = 0.12;               // Hill half-max for proliferation inflammation gating (vs 1.5 normal)

  // Efferocytosis-triggered M1-to-M2 (Extension 3)
  bool efferocytosis_enabled = false;          // proximity-based M1→M2
  real_t efferocytosis_radius = 5.0;           // search radius for dying neutrophils
  real_t efferocytosis_age_fraction = 0.8;     // neutrophil age fraction to be engulfable

  // Chemotaxis-driven immune migration (Extension 2)
  bool chemotaxis_enabled = false;             // follow inflammation gradient
  real_t chemotaxis_speed_scale = 1.0;         // gradient speed scaling factor

  // Inflammation-proportional scarring (Extension 1)
  bool scar_proportional_enabled = false;      // cumulative inflammation integral
  real_t scar_accumulation_rate = 0.001;       // scar += inflammation * rate per step
  real_t scar_collagen_threshold = 0.003;      // local collagen above this -> scar, below -> normal

  // Biofilm dynamics (Extension 7)
  bool biofilm_enabled = false;
  real_t biofilm_growth_rate = 0.005;          // per-step logistic growth
  real_t biofilm_carrying_capacity = 1.0;      // max biofilm density per voxel
  int biofilm_seed_delay = 480;                // TOML: 48h post-wound
  real_t biofilm_seed_amount = 0.01;           // initial inoculum density
  real_t biofilm_neutrophil_clearance = 0.003; // clearance per neutrophil per step
  real_t biofilm_macrophage_clearance = 0.001; // clearance per macrophage per step
  real_t biofilm_m1_block_threshold = 0.3;     // biofilm level blocking M1->M2
  real_t biofilm_inflammation_rate = 0.001;    // PAMP-driven inflammation per biofilm
  real_t diabetic_biofilm_clearance_factor = 0.5; // impaired neutrophil oxidative burst

  // VEGF-driven angiogenesis (Extension 8)
  bool angiogenesis_enabled = false;
  real_t vegf_diffusion = 0.08;                // VEGF diffusion coefficient
  real_t vegf_decay = 0.01;                    // VEGF natural decay rate
  real_t vegf_production_rate = 0.002;         // hypoxia-driven VEGF production
  real_t vegf_hypoxia_threshold = 0.5;         // O2 level below which VEGF is produced
  real_t vegf_consumption_rate = 0.1;          // VEGF consumed per unit perfusion recovery
  real_t angio_vegf_rate = 0.01;               // perfusion recovery per unit VEGF per step
  real_t m2_vegf_rate = 0.001;                 // VEGF produced by M2 macrophages per step
  real_t vegf_production_taper = 0.001;         // PHD2/3 negative feedback: exp(-k*wound_age) on VEGF production
  real_t vegf_receptor_clearance = 0.005;      // VEGFR-2 internalization rate per vascular density (Ferrara 2003)
  real_t diabetic_vegf_factor = 0.4;           // HIF-1alpha impairment
  real_t diabetic_ph_recovery_factor = 0.5;   // impaired acid mantle restoration

  // Fibroblast / TGF-beta / Collagen module (TOML: hours -> steps internally)
  bool fibroblast_enabled = false;
  int fibroblast_spawn_delay = 10;           // 1h post-wound
  int fibroblast_spawn_waves = 4;            // recruitment pulses
  int fibroblast_spawn_window = 960;         // 96h window
  real_t fibroblast_diameter = 4.0;
  real_t fibroblast_density_factor = 1.0;
  int fibroblast_activation_delay = 100;     // 10h before auto-activation
  real_t fibroblast_activation_threshold = 0.1;
  int fibroblast_myofibroblast_delay = 960;  // 96h min in activated
  real_t fibroblast_myofibroblast_threshold = 0.3;
  real_t fibroblast_apoptosis_threshold = 0.003;
  int fibroblast_apoptosis_onset = 1920;     // 192h as myofibroblast
  real_t fibroblast_apoptosis_rate = 0.0008; // per-step probability
  int fibroblast_min_lifespan = 1680;        // 168h (7d)
  int fibroblast_lifespan = 6720;            // 672h (28d)
  real_t fibroblast_migration_speed = 1.0;
  real_t dermal_fibroblast_depth = -3.0;    // z-coordinate for seeding (in dermis)
  real_t dermal_fibroblast_margin = 3.0;    // extra radius beyond wound for seeding ring
  real_t tgfb_diffusion = 0.03;
  real_t tgfb_decay = 0.005;
  real_t tgfb_wound_seed = 0.15;                // platelet alpha-granule TGF-beta bolus (Shah et al. 1995)
  real_t fibroblast_tgfb_rate = 0.005;
  real_t m2_tgfb_rate = 0.003;
  real_t collagen_deposition_rate = 0.002;
  real_t collagen_decay = 0.0;

  // Tumor module
  bool tumor_enabled = false;
  int tumor_seed_time = 0;            // TOML: hours; 0h = frame 0
  real_t tumor_seed_x = 15.0;        // cluster center X
  real_t tumor_seed_y = 15.0;        // cluster center Y
  real_t tumor_seed_z = 2.0;         // cluster center Z (basal layer)
  int tumor_seed_count = 5;          // initial cluster size
  real_t tumor_diameter = 4.0;       // cell diameter
  real_t tumor_cycle_factor = 3.3;   // phase duration multiplier; BCC cycle ~56h (Td=148d, Ki67=27%, phi=0.96)
  int tumor_max_neighbors = 12;      // contact inhibition threshold (sphere packing max ~12-14)
  real_t tumor_ci_steepness = 2.0;   // soft CI exponent (0 = hard threshold; 2 = quadratic YAP/TAZ)
  real_t tumor_growth_rate = 5;      // volume growth rate (same as normal)
  int tumor_max_cells = 5000;        // carrying capacity (0 = unlimited)
  int tumor_handoff_delay = 200;     // TOML: 20h in G0 before field conversion (0 = disabled)
  real_t tumor_stratum_value = 10.0; // value written to Stratum field on handoff
  real_t tumor_apoptosis_rate = 0.00047; // per-step probability (BCC cell loss factor)

  // MMP dynamics (matrix metalloproteinase -- ECM remodeling)
  // Lobmann et al. 2002 (doi:10.1007/s00125-002-0868-8): MMP-1 65x, MMP-9 14x in diabetic wounds
  bool mmp_enabled = false;
  real_t mmp_diffusion = 0.02;                // MMP diffusion coefficient
  real_t mmp_decay = 0.01;                    // TIMP-mediated clearance
  real_t mmp_m1_rate = 0.003;                 // MMP-9 from M1 macrophages per step
  real_t mmp_neutrophil_rate = 0.002;         // MMP-8 from neutrophils per step (Nwomeh 1998)
  real_t mmp_fibroblast_rate = 0.001;         // MMP-1/3 from activated fibroblasts per step
  real_t mmp_keratinocyte_rate = 0.0005;      // MMP-1 from wound-edge keratinocytes per step
  real_t mmp_collagen_degradation = 0.005;    // collagen loss per MMP unit per step
  real_t mmp_fibronectin_degradation = 0.002; // fibronectin loss per MMP unit per step
  real_t diabetic_mmp_factor = 3.0;           // MMP overexpression multiplier
  real_t diabetic_timp_factor = 0.5;          // reduced TIMP (increased effective MMP)

  // Fibronectin (provisional wound matrix -- migration scaffold)
  // Grinnell 1984 (doi:10.1002/jcb.240260206): fibronectin essential for wound migration
  bool fibronectin_enabled = false;
  real_t fibronectin_decay = 0.005;           // ECM turnover (half-life ~5.8 days)
  real_t fibronectin_deposition_rate = 0.003; // deposited by activated fibroblasts
  real_t fibronectin_serum_rate = 0.001;      // plasma fibronectin from wound serum
  real_t fibronectin_migration_boost = 0.5;   // max migration speed boost on fibronectin
  real_t fibronectin_wound_seed = 0.3;        // plasma fibronectin co-deposited with fibrin clot (Clark 1996)
  real_t fibronectin_carrying_capacity = 1.0; // local saturation limit for fibroblast deposition
  real_t collagen_fn_transition = 0.15;      // collagen level at which FN deposition ceases (ECM maturation)

  // Elastin (elastic fiber network)
  bool elastin_enabled = false;
  real_t elastin_diffusion = 0.0;
  real_t elastin_decay = 0.0002;
  real_t elastin_basal_density = 0.5;
  real_t elastin_papillary_density = 0.3;
  real_t elastin_production_rate = 0.0005;
  real_t elastin_mmp_degradation = 0.003;

  // Hyaluronan (HA / GAG ground substance)
  bool hyaluronan_enabled = false;
  real_t hyaluronan_diffusion = 0.01;
  real_t hyaluronan_decay = 0.005;
  real_t hyaluronan_basal_density = 0.8;
  real_t hyaluronan_reticular_density = 0.4;
  real_t hyaluronan_production_rate = 0.002;
  real_t hyaluronan_water_retention_factor = 0.5;
  real_t hyaluronan_migration_scaffold_factor = 0.3;

  // Dermis (dermal tissue integrity -- ECM structural state)
  bool dermis_enabled = false;
  real_t dermis_diffusion = 0.001;
  real_t dermis_decay = 0.0;
  real_t dermis_papillary_density = 1.0;
  real_t dermis_reticular_density = 0.8;
  real_t dermis_hypodermis_density = 0.5;
  real_t dermis_collagen_threshold = 0.05;
  real_t dermis_collagen_recovery_rate = 0.002;
  real_t dermis_mmp_degradation = 0.004;
  real_t dermis_papillary_rate_factor = 1.5;
  real_t dermis_reticular_rate_factor = 1.0;
  real_t dermis_hypodermis_rate_factor = 0.5;

  // pH (wound alkalinity -- acid mantle disruption)
  // Schneider et al. 2007 (doi:10.1111/j.1524-475X.2007.00230.x): acidic wound
  // bed (pH 5.5-6.5) promotes keratinocyte migration and re-epithelialization.
  // Field stores alkalinity excess: 0 = normal skin, 1.0 = fresh wound.
  real_t ph_recovery_rate = 0.003;             // acidification rate per step (perfusion-scaled)
  real_t ph_migration_suppression = 0.5;       // max speed reduction at full alkalinity
  real_t ph_mmp_boost = 0.5;                   // MMP activity amplification at full alkalinity
  real_t ph_biofilm_boost = 0.4;               // biofilm growth boost at full alkalinity
  real_t ph_bohr_factor = 0.3;                 // O2 delivery reduction at full alkalinity (Bohr effect)

  // Hemostasis / fibrin clot (provisional wound matrix)
  bool hemostasis_enabled = false;
  real_t hemostasis_decay = 0.0003;                // slow turnover
  real_t hemostasis_wound_seed = 1.0;              // initial clot density
  real_t hemostasis_tgfb_coupling = 0.0005;        // fibrin -> TGF-beta rate
  real_t hemostasis_fibronectin_coupling = 0.001;  // fibrin -> fibronectin rate
  real_t hemostasis_mmp_degradation = 0.002;       // MMP fibrinolysis rate
  real_t hemostasis_migration_boost = 0.3;         // keratinocyte migration boost

  // Derived composite fields (precomputed once per step from source grids)
  real_t ecm_weight_collagen = 0.4;       // collagen contribution to ECM quality
  real_t ecm_weight_fibronectin = 0.3;    // fibronectin contribution
  real_t ecm_weight_elastin = 0.2;        // elastin contribution
  real_t ecm_weight_fibrin = 0.1;         // fibrin provisional matrix contribution
  real_t ecm_migration_boost = 1.5;       // composite ECM migration speed boost

  // Multi-resolution: structural fields (D=0) at coarser grid
  // 0 = same as grid_resolution; nonzero = coarse resolution
  int grid_resolution_structural = 5;

  // Scar formation: dissolved agents stamp stratum+5 into Stratum field
  bool scar_enabled = true;

  // Metrics export (TOML: hours -> steps internally)
  int metrics_interval = 100;    // 10h between CSV rows (0 = disabled)
  bool metrics_autoopen = true;  // open CSV in default app after simulation

  // Cytokine taper rates (configurable exponential decay constants)
  real_t neutrophil_cytokine_taper = 0.002;  // neutrophil cytokine exp(-k*age)
  real_t m1_cytokine_taper = 0.002;          // M1 macrophage cytokine exp(-k*state_age)
  real_t m2_tgfb_taper = 0.003;              // M2 macrophage TGF-beta exp(-k*state_age)
  real_t fibroblast_tgfb_taper = 0.0012;     // myofibroblast TGF-beta exp(-k*state_age)

  // PDE sub-cycling: solve slow fields every N simulation steps instead of
  // every step. Reduces FTCS solver invocations by (N-1)/N. CFL-safe for
  // any D*subcycle*dt/dx^2 < 0.5 (verified at default parameters).
  // 1 = every step (no sub-cycling), 5 = every 5 steps, etc.
  int subcycle_slow = 1;           // water (D=1e-4), hyaluronan (D=0.01)
  int subcycle_medium = 1;         // inflammation (D=0.05), TGF-beta (D=0.03), MMP (D=0.02)

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
  void LoadConfig(const skibidy::TomlConfig& config);

  // Post-load sanity checks: abort early on obviously broken configs.
  void ValidateConfig() const;
};

}  // namespace skibidy
}  // namespace bdm

#endif  // SIM_PARAM_H_
