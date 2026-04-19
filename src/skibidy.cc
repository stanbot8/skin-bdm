// -----------------------------------------------------------------------------
//
// Copyright (C) 2021 CERN & University of Surrey for the benefit of the
// BioDynaMo collaboration. All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
//
// See the LICENSE file distributed with this work for details.
// See the NOTICE file distributed with this work for additional information
// regarding copyright ownership.
//
// -----------------------------------------------------------------------------
#include "infra/toml_compat.h"  // BDM_ASSIGN_CONFIG_VALUE, BDM_ASSIGN_HOURS/DAYS
#include "infra/sim_param.h"
#include "skibidy.h"

const bdm::ParamGroupUid bdm::skibidy::SimParam::kUid =
    bdm::ParamGroupUidGenerator::Get()->NewUid();

void bdm::skibidy::SimParam::LoadConfig(const skibidy::TomlConfig& config) {
  // ---- Geometry: derive all spatial bounds from physical dimensions ----
  BDM_ASSIGN_CONFIG_VALUE(geo_patch_um, "geometry.patch_um");
  BDM_ASSIGN_CONFIG_VALUE(geo_depth_um, "geometry.depth_um");
  BDM_ASSIGN_CONFIG_VALUE(geo_height_um, "geometry.height_um");
  BDM_ASSIGN_CONFIG_VALUE(geo_margin_um, "geometry.margin_um");
  BDM_ASSIGN_CONFIG_VALUE(geo_um_per_unit, "geometry.um_per_unit");
  BDM_ASSIGN_CONFIG_VALUE(geo_voxel_um, "geometry.voxel_um");

  // Derive spatial bounds from geometry
  real_t tissue_half = geo_patch_um / (2.0 * geo_um_per_unit);
  real_t margin = geo_margin_um / geo_um_per_unit;
  tissue_min = -tissue_half;
  tissue_max = tissue_half;
  domain_min = tissue_min - margin;
  domain_max = tissue_max + margin;
  domain_z_min = -(geo_depth_um / geo_um_per_unit);
  domain_z_max = geo_height_um / geo_um_per_unit;

  // Derive grid resolution and box length from voxel size
  ref_box_length = geo_voxel_um / geo_um_per_unit;
  real_t domain_size = domain_max - domain_min;
  grid_resolution = static_cast<int>(std::round(domain_size / ref_box_length));
  if (grid_resolution < 1) grid_resolution = 1;

  // Update per-axis bounds (x/y = tissue range, z = full domain)
  bounds_min = {tissue_min, tissue_min, domain_z_min};
  bounds_max = {tissue_max, tissue_max, domain_z_max};

  // Cell cycle
  BDM_ASSIGN_CONFIG_VALUE(g1_duration, "skin.g1_duration");
  BDM_ASSIGN_CONFIG_VALUE(s_duration, "skin.s_duration");
  BDM_ASSIGN_CONFIG_VALUE(g2_duration, "skin.g2_duration");
  BDM_ASSIGN_CONFIG_VALUE(m_duration, "skin.m_duration");

  // Division mechanics
  BDM_ASSIGN_CONFIG_VALUE(growth_rate, "skin.growth_rate");
  BDM_ASSIGN_CONFIG_VALUE(max_neighbors, "skin.max_neighbors");
  BDM_ASSIGN_CONFIG_VALUE(lateral_scatter, "skin.lateral_scatter");
  BDM_ASSIGN_CONFIG_VALUE(max_ta_divisions, "skin.max_ta_divisions");
  BDM_ASSIGN_CONFIG_VALUE(division_diameter, "skin.division_diameter");
  BDM_ASSIGN_CONFIG_VALUE(p_asymmetric, "skin.p_asymmetric");
  BDM_ASSIGN_CONFIG_VALUE(stem_fraction, "skin.stem_fraction");

  // Calcium gradient
  BDM_ASSIGN_CONFIG_VALUE(calcium_diffusion, "skin.calcium_diffusion");
  BDM_ASSIGN_CONFIG_VALUE(calcium_decay, "skin.calcium_decay");
  BDM_ASSIGN_CONFIG_VALUE(calcium_basal, "skin.calcium_basal");
  BDM_ASSIGN_CONFIG_VALUE(calcium_peak, "skin.calcium_peak");
  BDM_ASSIGN_CONFIG_VALUE(calcium_midpoint_z, "skin.calcium_midpoint_z");
  BDM_ASSIGN_CONFIG_VALUE(calcium_steepness, "skin.calcium_steepness");

  // Differentiation thresholds
  BDM_ASSIGN_CONFIG_VALUE(ca_spinous_threshold, "skin.ca_spinous_threshold");
  BDM_ASSIGN_CONFIG_VALUE(ca_granular_threshold, "skin.ca_granular_threshold");
  BDM_ASSIGN_CONFIG_VALUE(ca_cornified_threshold, "skin.ca_cornified_threshold");
  BDM_ASSIGN_CONFIG_VALUE(spinous_threshold, "skin.spinous_threshold");

  // Layer thicknesses (relative to basement membrane at z=0)
  BDM_ASSIGN_CONFIG_VALUE(basal_thickness, "skin.basal_thickness");
  BDM_ASSIGN_CONFIG_VALUE(spinous_thickness, "skin.spinous_thickness");
  BDM_ASSIGN_CONFIG_VALUE(granular_thickness, "skin.granular_thickness");
  BDM_ASSIGN_CONFIG_VALUE(papillary_thickness, "skin.papillary_thickness");
  BDM_ASSIGN_CONFIG_VALUE(reticular_thickness, "skin.reticular_thickness");

  // Compute z-boundaries from thicknesses
  volume_z_spinous = basal_thickness;
  volume_z_granular = basal_thickness + spinous_thickness;
  volume_z_cornified = basal_thickness + spinous_thickness + granular_thickness;
  dermal_z_papillary = -papillary_thickness;
  dermal_z_reticular = -(papillary_thickness + reticular_thickness);

  // Cell-cell mechanics
  BDM_ASSIGN_CONFIG_VALUE(repulsion_coeff, "skin.repulsion_coeff");
  BDM_ASSIGN_CONFIG_VALUE(attraction_coeff, "skin.attraction_coeff");

  // KGF (dermal growth factor)
  BDM_ASSIGN_CONFIG_VALUE(kgf_diffusion, "skin.kgf_diffusion");
  BDM_ASSIGN_CONFIG_VALUE(kgf_decay, "skin.kgf_decay");
  BDM_ASSIGN_CONFIG_VALUE(kgf_basal_conc, "skin.kgf_basal_conc");
  BDM_ASSIGN_CONFIG_VALUE(kgf_half_maximal, "skin.kgf_half_maximal");
  BDM_ASSIGN_CONFIG_VALUE(kgf_max_boost, "skin.kgf_max_boost");
  BDM_ASSIGN_CONFIG_VALUE(kgf_decay_length, "skin.kgf_decay_length");

  // Oxygen
  BDM_ASSIGN_CONFIG_VALUE(oxygen_diffusion, "skin.oxygen_diffusion");
  BDM_ASSIGN_CONFIG_VALUE(oxygen_decay, "skin.oxygen_decay");
  BDM_ASSIGN_CONFIG_VALUE(oxygen_basal_conc, "skin.oxygen_basal_conc");
  BDM_ASSIGN_CONFIG_VALUE(oxygen_decay_length, "skin.oxygen_decay_length");

  // Water / tissue moisture
  BDM_ASSIGN_CONFIG_VALUE(water_diffusion, "skin.water_diffusion");
  BDM_ASSIGN_CONFIG_VALUE(water_decay, "skin.water_decay");
  BDM_ASSIGN_CONFIG_VALUE(water_basal_conc, "skin.water_basal_conc");
  BDM_ASSIGN_CONFIG_VALUE(water_decay_length, "skin.water_decay_length");
  BDM_ASSIGN_CONFIG_VALUE(water_recovery_rate, "skin.water_recovery_rate");
  BDM_ASSIGN_CONFIG_VALUE(water_surface_loss_rate, "skin.water_surface_loss_rate");
  BDM_ASSIGN_CONFIG_VALUE(water_migration_threshold, "skin.water_migration_threshold");
  BDM_ASSIGN_CONFIG_VALUE(water_prolif_threshold, "skin.water_prolif_threshold");

  // Shedding & apoptosis (TOML: hours)
  BDM_ASSIGN_HOURS(shedding_delay, "skin.shedding_delay_h");
  BDM_ASSIGN_HOURS(apoptosis_delay, "skin.apoptosis_delay_h");

  // Tissue extent: now derived from [geometry]; legacy keys ignored.
  // (tissue_min, tissue_max, grid_resolution, ref_box_length set above)

  // Simulation (TOML: days)
  BDM_ASSIGN_DAYS(num_steps, "skin.duration_days");

  // Wound event (modules/wound.toml -> [skin.wound])
  BDM_ASSIGN_CONFIG_VALUE(wound.enabled, "skin.wound.enabled");
  BDM_ASSIGN_CONFIG_VALUE(wound.center_x, "skin.wound.center_x");
  BDM_ASSIGN_CONFIG_VALUE(wound.center_y, "skin.wound.center_y");
  BDM_ASSIGN_CONFIG_VALUE(wound.radius, "skin.wound.radius");
  BDM_ASSIGN_HOURS(wound.trigger_step, "skin.wound.trigger_h");
  BDM_ASSIGN_CONFIG_VALUE(wound.inward_bias, "skin.wound.inward_bias");
  BDM_ASSIGN_HOURS(homeostatic_fold, "skin.wound.homeostatic_fold_h");
  BDM_ASSIGN_CONFIG_VALUE(dissolution_closure_pct, "skin.wound.dissolution_closure_pct");

  // Vascular perfusion (modules/perfusion.toml -> [skin.perfusion])
  BDM_ASSIGN_CONFIG_VALUE(perfusion.diffusion, "skin.perfusion.diffusion");
  BDM_ASSIGN_CONFIG_VALUE(perfusion.decay, "skin.perfusion.decay");
  BDM_ASSIGN_CONFIG_VALUE(perfusion.basal, "skin.perfusion.basal");
  BDM_ASSIGN_CONFIG_VALUE(perfusion.angio_rate, "skin.perfusion.angio_rate");
  BDM_ASSIGN_HOURS(perfusion.angio_delay, "skin.perfusion.angio_delay_h");
  BDM_ASSIGN_CONFIG_VALUE(perfusion.papillary_fraction, "skin.perfusion.papillary_fraction");
  BDM_ASSIGN_CONFIG_VALUE(perfusion.reticular_fraction, "skin.perfusion.reticular_fraction");
  BDM_ASSIGN_CONFIG_VALUE(perfusion.hypodermis_fraction, "skin.perfusion.hypodermis_fraction");
  BDM_ASSIGN_CONFIG_VALUE(angio_papillary_factor, "skin.perfusion.angio_papillary_factor");
  BDM_ASSIGN_CONFIG_VALUE(angio_reticular_factor, "skin.perfusion.angio_reticular_factor");
  BDM_ASSIGN_CONFIG_VALUE(angio_hypodermis_factor, "skin.perfusion.angio_hypodermis_factor");

  // Dynamic field coupling
  BDM_ASSIGN_CONFIG_VALUE(calcium_recovery_rate, "skin.calcium_recovery_rate");
  BDM_ASSIGN_CONFIG_VALUE(oxygen_prolif_threshold, "skin.oxygen_prolif_threshold");
  BDM_ASSIGN_CONFIG_VALUE(oxygen_recovery_enabled, "skin.oxygen_recovery_enabled");

  // Cell migration -- multi-cue mechanistic
  BDM_ASSIGN_CONFIG_VALUE(migration_enabled, "skin.migration_enabled");
  BDM_ASSIGN_CONFIG_VALUE(migration_speed, "skin.migration_speed");
  BDM_ASSIGN_CONFIG_VALUE(migration_gradient_scale, "skin.migration_gradient_scale");
  BDM_ASSIGN_CONFIG_VALUE(angiogenesis.vegf_migration_weight, "skin.vegf_migration_weight");
  BDM_ASSIGN_CONFIG_VALUE(fibronectin.haptotaxis_weight, "skin.fibronectin_haptotaxis_weight");
  BDM_ASSIGN_CONFIG_VALUE(cil_suppression, "skin.cil_suppression");

  // Per-cell UWYN handoff (TOML: hours)
  BDM_ASSIGN_HOURS(handoff_delay, "skin.handoff_delay_h");

  // Inflammation (modules/inflammation.toml -> [skin.inflammation])
  BDM_ASSIGN_CONFIG_VALUE(inflammation.diffusion, "skin.inflammation.diffusion");
  BDM_ASSIGN_CONFIG_VALUE(inflammation.decay, "skin.inflammation.decay");
  BDM_ASSIGN_CONFIG_VALUE(inflammation.migration_threshold, "skin.inflammation.migration_threshold");
  BDM_ASSIGN_CONFIG_VALUE(inflammation.prolif_threshold, "skin.inflammation.prolif_threshold");
  BDM_ASSIGN_CONFIG_VALUE(wound.inflammation_source_rate, "skin.inflammation.wound_source_rate");
  BDM_ASSIGN_CONFIG_VALUE(wound.inflammation_source_taper, "skin.inflammation.wound_source_taper");
  BDM_ASSIGN_CONFIG_VALUE(inflammation.split_inflammation_enabled, "skin.inflammation.split_enabled");
  BDM_ASSIGN_CONFIG_VALUE(anti_inflammation_diffusion, "skin.inflammation.anti_diffusion");
  BDM_ASSIGN_CONFIG_VALUE(anti_inflammation_decay, "skin.inflammation.anti_decay");
  BDM_ASSIGN_CONFIG_VALUE(anti_inflammation_weight, "skin.inflammation.anti_weight");

  // Immune response (modules/immune.toml -> [skin.immune]; TOML: hours)
  BDM_ASSIGN_HOURS(immune.neutrophil_spawn_delay, "skin.immune.neutrophil_spawn_delay_h");
  BDM_ASSIGN_CONFIG_VALUE(immune.neutrophil_spawn_waves, "skin.immune.neutrophil_spawn_waves");
  BDM_ASSIGN_HOURS(immune.neutrophil_spawn_window, "skin.immune.neutrophil_spawn_window_h");
  BDM_ASSIGN_HOURS(immune.neutrophil_lifespan, "skin.immune.neutrophil_lifespan_h");
  BDM_ASSIGN_HOURS(immune.neutrophil_min_survival, "skin.immune.neutrophil_min_survival_h");
  BDM_ASSIGN_CONFIG_VALUE(immune.neutrophil_apoptosis_rate, "skin.immune.neutrophil_apoptosis_rate");
  BDM_ASSIGN_HOURS(immune.macrophage_spawn_delay, "skin.immune.macrophage_spawn_delay_h");
  BDM_ASSIGN_CONFIG_VALUE(immune.macrophage_spawn_threshold, "skin.immune.macrophage_spawn_threshold");
  BDM_ASSIGN_CONFIG_VALUE(immune.macrophage_spawn_rate, "skin.immune.macrophage_spawn_rate");
  BDM_ASSIGN_CONFIG_VALUE(immune.macrophage_spawn_taper, "skin.immune.macrophage_spawn_taper");
  BDM_ASSIGN_HOURS(immune.macrophage_m1_duration, "skin.immune.macrophage_m1_duration_h");
  BDM_ASSIGN_CONFIG_VALUE(m1_transition_threshold, "skin.immune.m1_transition_threshold");
  BDM_ASSIGN_HOURS(m1_transition_min_age, "skin.immune.m1_transition_min_age_h");
  BDM_ASSIGN_HOURS(immune.macrophage_lifespan, "skin.immune.macrophage_lifespan_h");
  BDM_ASSIGN_HOURS(immune.macrophage_min_survival, "skin.immune.macrophage_min_survival_h");
  BDM_ASSIGN_CONFIG_VALUE(immune.macrophage_apoptosis_rate, "skin.immune.macrophage_apoptosis_rate");
  BDM_ASSIGN_CONFIG_VALUE(immune.macrophage_emigration_rate, "skin.immune.macrophage_emigration_rate");
  BDM_ASSIGN_CONFIG_VALUE(immune_cell_diameter, "skin.immune.cell_diameter");
  BDM_ASSIGN_CONFIG_VALUE(immune.cytokine_rate, "skin.immune.cytokine_rate");
  BDM_ASSIGN_CONFIG_VALUE(immune.resolution_rate, "skin.immune.resolution_rate");
  BDM_ASSIGN_CONFIG_VALUE(immune.migration_speed, "skin.immune.migration_speed");
  BDM_ASSIGN_CONFIG_VALUE(efferocytosis_enabled, "skin.immune.efferocytosis_enabled");
  BDM_ASSIGN_CONFIG_VALUE(efferocytosis_radius, "skin.immune.efferocytosis_radius");
  BDM_ASSIGN_CONFIG_VALUE(efferocytosis_age_fraction, "skin.immune.efferocytosis_age_fraction");
  BDM_ASSIGN_CONFIG_VALUE(chemotaxis_enabled, "skin.immune.chemotaxis_enabled");
  BDM_ASSIGN_CONFIG_VALUE(chemotaxis_speed_scale, "skin.immune.chemotaxis_speed_scale");

  // Mechanistic immune replacements
  BDM_ASSIGN_CONFIG_VALUE(mech_immune_recruitment, "skin.immune.mech_immune_recruitment");
  BDM_ASSIGN_CONFIG_VALUE(mech_recruit_gradient_scale, "skin.immune.mech_recruit_gradient_scale");
  BDM_ASSIGN_CONFIG_VALUE(mech_recruit_saturation_k, "skin.immune.mech_recruit_saturation_k");
  BDM_ASSIGN_CONFIG_VALUE(mech_m1_m2_transition, "skin.immune.mech_m1_m2_transition");
  BDM_ASSIGN_CONFIG_VALUE(mech_efferocytosis_quota, "skin.immune.mech_efferocytosis_quota");
  BDM_ASSIGN_CONFIG_VALUE(mech_m2_transition_rate, "skin.immune.mech_m2_transition_rate");

  // Biofilm dynamics (modules/biofilm.toml -> [skin.biofilm])
  BDM_ASSIGN_CONFIG_VALUE(biofilm.enabled, "skin.biofilm.enabled");
  BDM_ASSIGN_CONFIG_VALUE(biofilm.growth_rate, "skin.biofilm.growth_rate");
  BDM_ASSIGN_CONFIG_VALUE(biofilm.carrying_capacity, "skin.biofilm.carrying_capacity");
  BDM_ASSIGN_HOURS(biofilm.seed_delay, "skin.biofilm.seed_delay_h");
  BDM_ASSIGN_CONFIG_VALUE(biofilm.seed_amount, "skin.biofilm.seed_amount");
  BDM_ASSIGN_CONFIG_VALUE(biofilm.neutrophil_clearance, "skin.biofilm.neutrophil_clearance");
  BDM_ASSIGN_CONFIG_VALUE(biofilm.macrophage_clearance, "skin.biofilm.macrophage_clearance");
  BDM_ASSIGN_CONFIG_VALUE(biofilm.m1_block_threshold, "skin.biofilm.m1_block_threshold");
  BDM_ASSIGN_CONFIG_VALUE(biofilm.inflammation_rate, "skin.biofilm.inflammation_rate");

  // VEGF-driven angiogenesis (modules/angiogenesis.toml -> [skin.angiogenesis])
  BDM_ASSIGN_CONFIG_VALUE(angiogenesis.enabled, "skin.angiogenesis.enabled");
  BDM_ASSIGN_CONFIG_VALUE(angiogenesis.vegf_diffusion, "skin.angiogenesis.vegf_diffusion");
  BDM_ASSIGN_CONFIG_VALUE(angiogenesis.vegf_decay, "skin.angiogenesis.vegf_decay");
  BDM_ASSIGN_CONFIG_VALUE(angiogenesis.vegf_production_rate, "skin.angiogenesis.vegf_production_rate");
  BDM_ASSIGN_CONFIG_VALUE(angiogenesis.vegf_hypoxia_threshold, "skin.angiogenesis.vegf_hypoxia_threshold");
  BDM_ASSIGN_CONFIG_VALUE(angiogenesis.vegf_consumption_rate, "skin.angiogenesis.vegf_consumption_rate");
  BDM_ASSIGN_CONFIG_VALUE(angio_vegf_rate, "skin.angiogenesis.angio_vegf_rate");
  BDM_ASSIGN_CONFIG_VALUE(m2_vegf_rate, "skin.angiogenesis.m2_vegf_rate");
  BDM_ASSIGN_CONFIG_VALUE(angiogenesis.vegf_production_taper, "skin.angiogenesis.vegf_production_taper");
  BDM_ASSIGN_CONFIG_VALUE(angiogenesis.vegf_receptor_clearance, "skin.angiogenesis.vegf_receptor_clearance");

  // Mechanistic VEGF replacement
  BDM_ASSIGN_CONFIG_VALUE(mech_vegf_production, "skin.angiogenesis.mech_vegf_production");
  BDM_ASSIGN_CONFIG_VALUE(mech_hif_o2_threshold, "skin.angiogenesis.mech_hif_o2_threshold");
  BDM_ASSIGN_CONFIG_VALUE(mech_hif_vegf_rate, "skin.angiogenesis.mech_hif_vegf_rate");

  // Diabetic modifiers (modules/diabetic.toml -> [skin.diabetic])
  BDM_ASSIGN_CONFIG_VALUE(diabetic.mode, "skin.diabetic.mode");
  BDM_ASSIGN_CONFIG_VALUE(diabetic.m1_duration_factor, "skin.diabetic.m1_duration_factor");
  BDM_ASSIGN_CONFIG_VALUE(diabetic.resolution_factor, "skin.diabetic.resolution_factor");
  BDM_ASSIGN_CONFIG_VALUE(diabetic.efferocytosis_factor, "skin.diabetic.efferocytosis_factor");
  BDM_ASSIGN_CONFIG_VALUE(diabetic.prolif_factor, "skin.diabetic.prolif_factor");
  BDM_ASSIGN_CONFIG_VALUE(diabetic.migration_factor, "skin.diabetic.migration_factor");
  BDM_ASSIGN_CONFIG_VALUE(diabetic.fibroblast_activation_factor, "skin.diabetic.fibroblast_activation_factor");
  BDM_ASSIGN_CONFIG_VALUE(diabetic.collagen_factor, "skin.diabetic.collagen_factor");
  BDM_ASSIGN_CONFIG_VALUE(diabetic.fibroblast_lifespan_factor, "skin.diabetic.fibroblast_lifespan_factor");
  BDM_ASSIGN_CONFIG_VALUE(diabetic.baseline_inflammation, "skin.diabetic.baseline_inflammation");
  BDM_ASSIGN_CONFIG_VALUE(diabetic.neutrophil_factor, "skin.diabetic.neutrophil_factor");
  BDM_ASSIGN_CONFIG_VALUE(diabetic.neutrophil_lifespan_factor, "skin.diabetic.neutrophil_lifespan_factor");
  BDM_ASSIGN_CONFIG_VALUE(diabetic.neutrophil_apoptosis_factor, "skin.diabetic.neutrophil_apoptosis_factor");
  BDM_ASSIGN_CONFIG_VALUE(diabetic.neutrophil_waves_factor, "skin.diabetic.neutrophil_waves_factor");
  BDM_ASSIGN_CONFIG_VALUE(diabetic.neutrophil_window_factor, "skin.diabetic.neutrophil_window_factor");
  BDM_ASSIGN_CONFIG_VALUE(diabetic.macrophage_taper_factor, "skin.diabetic.macrophage_taper_factor");
  BDM_ASSIGN_CONFIG_VALUE(diabetic.macrophage_apoptosis_factor, "skin.diabetic.macrophage_apoptosis_factor");
  BDM_ASSIGN_CONFIG_VALUE(diabetic.macrophage_emigration_factor, "skin.diabetic.macrophage_emigration_factor");
  BDM_ASSIGN_CONFIG_VALUE(diabetic.biofilm_clearance_factor, "skin.diabetic.biofilm_clearance_factor");
  BDM_ASSIGN_CONFIG_VALUE(diabetic.vegf_factor, "skin.diabetic.vegf_factor");
  BDM_ASSIGN_CONFIG_VALUE(diabetic.ph_recovery_factor, "skin.diabetic.ph_recovery_factor");
  BDM_ASSIGN_CONFIG_VALUE(diabetic.inflammation_sensitivity, "skin.diabetic.inflammation_sensitivity");
  BDM_ASSIGN_CONFIG_VALUE(diabetic.migration_infl_K, "skin.diabetic.migration_infl_K");
  BDM_ASSIGN_CONFIG_VALUE(diabetic.prolif_infl_K, "skin.diabetic.prolif_infl_K");
  BDM_ASSIGN_CONFIG_VALUE(diabetic.mmp_factor, "skin.diabetic.mmp_factor");
  BDM_ASSIGN_CONFIG_VALUE(diabetic.timp_production_factor, "skin.diabetic.timp_production_factor");

  // Fibroblast / TGF-beta / Collagen (modules/fibroblast.toml -> [skin.fibroblast]; TOML: hours)
  BDM_ASSIGN_CONFIG_VALUE(fibroblast.enabled, "skin.fibroblast.enabled");
  BDM_ASSIGN_HOURS(fibroblast.spawn_delay, "skin.fibroblast.spawn_delay_h");
  BDM_ASSIGN_CONFIG_VALUE(fibroblast.spawn_waves, "skin.fibroblast.spawn_waves");
  BDM_ASSIGN_HOURS(fibroblast.spawn_window, "skin.fibroblast.spawn_window_h");
  BDM_ASSIGN_CONFIG_VALUE(fibroblast.diameter, "skin.fibroblast.diameter");
  BDM_ASSIGN_CONFIG_VALUE(fibroblast.density_factor, "skin.fibroblast.density_factor");
  BDM_ASSIGN_HOURS(fibroblast.activation_delay, "skin.fibroblast.activation_delay_h");
  BDM_ASSIGN_CONFIG_VALUE(fibroblast.activation_threshold, "skin.fibroblast.activation_threshold");
  BDM_ASSIGN_HOURS(fibroblast.myofibroblast_delay, "skin.fibroblast.myofibroblast_delay_h");
  BDM_ASSIGN_CONFIG_VALUE(fibroblast.myofibroblast_threshold, "skin.fibroblast.myofibroblast_threshold");
  BDM_ASSIGN_CONFIG_VALUE(fibroblast.apoptosis_threshold, "skin.fibroblast.apoptosis_threshold");
  BDM_ASSIGN_HOURS(fibroblast.apoptosis_onset, "skin.fibroblast.apoptosis_onset_h");
  BDM_ASSIGN_CONFIG_VALUE(fibroblast.apoptosis_rate, "skin.fibroblast.apoptosis_rate");
  BDM_ASSIGN_HOURS(fibroblast.min_lifespan, "skin.fibroblast.min_lifespan_h");
  BDM_ASSIGN_HOURS(fibroblast.lifespan, "skin.fibroblast.lifespan_h");
  BDM_ASSIGN_HOURS(fibroblast.activated_apoptosis_onset, "skin.fibroblast.activated_apoptosis_onset_h");
  BDM_ASSIGN_CONFIG_VALUE(fibroblast.activated_apoptosis_rate, "skin.fibroblast.activated_apoptosis_rate");
  BDM_ASSIGN_CONFIG_VALUE(fibroblast.migration_speed, "skin.fibroblast.migration_speed");
  BDM_ASSIGN_CONFIG_VALUE(dermal_fibroblast_depth, "skin.fibroblast.dermal_depth");
  BDM_ASSIGN_CONFIG_VALUE(dermal_fibroblast_margin, "skin.fibroblast.dermal_margin");
  BDM_ASSIGN_CONFIG_VALUE(fibroblast.tgfb_diffusion, "skin.fibroblast.tgfb_diffusion");
  BDM_ASSIGN_CONFIG_VALUE(fibroblast.tgfb_decay, "skin.fibroblast.tgfb_decay");
  BDM_ASSIGN_CONFIG_VALUE(fibroblast.tgfb_wound_seed, "skin.fibroblast.tgfb_wound_seed");
  BDM_ASSIGN_CONFIG_VALUE(fibroblast.tgfb_rate, "skin.fibroblast.tgfb_rate");
  BDM_ASSIGN_CONFIG_VALUE(m2_tgfb_rate, "skin.fibroblast.m2_tgfb_rate");
  BDM_ASSIGN_CONFIG_VALUE(fibroblast.collagen_deposition_rate, "skin.fibroblast.collagen_deposition_rate");
  BDM_ASSIGN_CONFIG_VALUE(fibroblast.collagen_decay, "skin.fibroblast.collagen_decay");
  BDM_ASSIGN_CONFIG_VALUE(decorin_sequestration_rate, "skin.fibroblast.decorin_sequestration_rate");
  BDM_ASSIGN_CONFIG_VALUE(fibroblast.tgfb_receptor_consumption, "skin.fibroblast.tgfb_receptor_consumption");
  BDM_ASSIGN_CONFIG_VALUE(fibroblast.tgfb_tissue_clearance, "skin.fibroblast.tgfb_tissue_clearance");
  BDM_ASSIGN_CONFIG_VALUE(perfusion.clearance_rate, "skin.fibroblast.perfusion_clearance_rate");
  BDM_ASSIGN_CONFIG_VALUE(fibroblast.collagen_o2_half_max, "skin.fibroblast.collagen_o2_half_max");

  // Mechanistic collagen replacement
  BDM_ASSIGN_CONFIG_VALUE(mech_collagen_deposition, "skin.fibroblast.mech_collagen_deposition");
  BDM_ASSIGN_CONFIG_VALUE(mech_collagen_tgfb_km, "skin.fibroblast.mech_collagen_tgfb_km");
  BDM_ASSIGN_CONFIG_VALUE(mech_collagen_vmax, "skin.fibroblast.mech_collagen_vmax");
  BDM_ASSIGN_CONFIG_VALUE(mech_collagen_basal, "skin.fibroblast.mech_collagen_basal");

  // Tumor module (modules/tumor.toml -> [skin.tumor])
  BDM_ASSIGN_CONFIG_VALUE(tumor.enabled, "skin.tumor.enabled");
  BDM_ASSIGN_HOURS(tumor.seed_time, "skin.tumor.seed_time_h");
  BDM_ASSIGN_CONFIG_VALUE(tumor.seed_x, "skin.tumor.seed_x");
  BDM_ASSIGN_CONFIG_VALUE(tumor.seed_y, "skin.tumor.seed_y");
  BDM_ASSIGN_CONFIG_VALUE(tumor.seed_z, "skin.tumor.seed_z");
  BDM_ASSIGN_CONFIG_VALUE(tumor.seed_count, "skin.tumor.seed_count");
  BDM_ASSIGN_CONFIG_VALUE(tumor.diameter, "skin.tumor.diameter");
  BDM_ASSIGN_CONFIG_VALUE(tumor.cycle_factor, "skin.tumor.cycle_factor");
  BDM_ASSIGN_CONFIG_VALUE(tumor.g1_factor, "skin.tumor.g1_factor");
  BDM_ASSIGN_CONFIG_VALUE(tumor.max_neighbors, "skin.tumor.max_neighbors");
  BDM_ASSIGN_CONFIG_VALUE(tumor.ci_steepness, "skin.tumor.ci_steepness");
  BDM_ASSIGN_CONFIG_VALUE(tumor.growth_rate, "skin.tumor.growth_rate");
  BDM_ASSIGN_CONFIG_VALUE(tumor.max_cells, "skin.tumor.max_cells");
  BDM_ASSIGN_HOURS(tumor.handoff_delay, "skin.tumor.handoff_delay_h");
  BDM_ASSIGN_CONFIG_VALUE(tumor.stratum_value, "skin.tumor.stratum_value");
  BDM_ASSIGN_CONFIG_VALUE(tumor.apoptosis_rate, "skin.tumor.apoptosis_rate");
  BDM_ASSIGN_CONFIG_VALUE(tumor.o2_threshold, "skin.tumor.o2_threshold");

  // MMP dynamics (modules/mmp.toml -> [skin.mmp])
  BDM_ASSIGN_CONFIG_VALUE(mmp.enabled, "skin.mmp.enabled");
  BDM_ASSIGN_CONFIG_VALUE(mmp.diffusion, "skin.mmp.diffusion");
  BDM_ASSIGN_CONFIG_VALUE(mmp.m1_rate, "skin.mmp.m1_rate");
  BDM_ASSIGN_CONFIG_VALUE(mmp.neutrophil_rate, "skin.mmp.neutrophil_rate");
  BDM_ASSIGN_CONFIG_VALUE(mmp.fibroblast_rate, "skin.mmp.fibroblast_rate");
  BDM_ASSIGN_CONFIG_VALUE(mmp.keratinocyte_rate, "skin.mmp.keratinocyte_rate");
  BDM_ASSIGN_CONFIG_VALUE(mmp.collagen_degradation, "skin.mmp.collagen_degradation");
  BDM_ASSIGN_CONFIG_VALUE(mmp.fibronectin_degradation, "skin.mmp.fibronectin_degradation");
  BDM_ASSIGN_CONFIG_VALUE(mmp.residual_decay, "skin.mmp.residual_decay");
  // TIMP dynamics
  BDM_ASSIGN_CONFIG_VALUE(mmp.timp_diffusion, "skin.mmp.timp_diffusion");
  BDM_ASSIGN_CONFIG_VALUE(mmp.timp_decay, "skin.mmp.timp_decay");
  BDM_ASSIGN_CONFIG_VALUE(mmp.timp_fibroblast_rate, "skin.mmp.timp_fibroblast_rate");
  BDM_ASSIGN_CONFIG_VALUE(mmp.timp_macrophage_rate, "skin.mmp.timp_macrophage_rate");
  BDM_ASSIGN_CONFIG_VALUE(mmp.timp_keratinocyte_rate, "skin.mmp.timp_keratinocyte_rate");
  BDM_ASSIGN_CONFIG_VALUE(mmp.timp_inhibition_rate, "skin.mmp.mmp_timp_inhibition_rate");
  BDM_ASSIGN_CONFIG_VALUE(matrikine_mmp_boost, "skin.mmp.matrikine_mmp_boost");
  // Pro-MMP zymogen cascade
  BDM_ASSIGN_CONFIG_VALUE(mmp.prommp_decay, "skin.mmp.prommp_decay");
  BDM_ASSIGN_CONFIG_VALUE(mmp.prommp_activation_rate, "skin.mmp.prommp_activation_rate");
  BDM_ASSIGN_CONFIG_VALUE(mmp.prommp_autocatalytic_rate, "skin.mmp.prommp_autocatalytic_rate");

  // Fibronectin (modules/fibronectin.toml -> [skin.fibronectin])
  BDM_ASSIGN_CONFIG_VALUE(fibronectin.enabled, "skin.fibronectin.enabled");
  BDM_ASSIGN_CONFIG_VALUE(fibronectin.decay, "skin.fibronectin.decay");
  BDM_ASSIGN_CONFIG_VALUE(fibronectin.deposition_rate, "skin.fibronectin.deposition_rate");
  BDM_ASSIGN_CONFIG_VALUE(fibronectin.serum_rate, "skin.fibronectin.serum_rate");
  BDM_ASSIGN_CONFIG_VALUE(fibronectin.migration_boost, "skin.fibronectin.migration_boost");
  BDM_ASSIGN_CONFIG_VALUE(fibronectin.wound_seed, "skin.fibronectin.wound_seed");
  BDM_ASSIGN_CONFIG_VALUE(fibronectin.carrying_capacity, "skin.fibronectin.carrying_capacity");
  BDM_ASSIGN_CONFIG_VALUE(fibroblast.collagen_fn_transition, "skin.fibronectin.collagen_fn_transition");
  BDM_ASSIGN_CONFIG_VALUE(fn_tgfb_max_boost, "skin.fibronectin.tgfb_max_boost");
  BDM_ASSIGN_CONFIG_VALUE(fn_tgfb_half_maximal, "skin.fibronectin.tgfb_half_maximal");
  BDM_ASSIGN_CONFIG_VALUE(myofibroblast_fn_fraction, "skin.fibronectin.myofibroblast_fn_fraction");

  // Elastin (modules/elastin/config.toml -> [skin.elastin])
  BDM_ASSIGN_CONFIG_VALUE(elastin.enabled, "skin.elastin.enabled");
  BDM_ASSIGN_CONFIG_VALUE(elastin.diffusion, "skin.elastin.diffusion");
  BDM_ASSIGN_CONFIG_VALUE(elastin.decay, "skin.elastin.decay");
  BDM_ASSIGN_CONFIG_VALUE(elastin.basal_density, "skin.elastin.basal_density");
  BDM_ASSIGN_CONFIG_VALUE(elastin.papillary_density, "skin.elastin.papillary_density");
  BDM_ASSIGN_CONFIG_VALUE(elastin.production_rate, "skin.elastin.production_rate");
  BDM_ASSIGN_CONFIG_VALUE(elastin.mmp_degradation, "skin.elastin.mmp_degradation");

  // Hyaluronan (modules/hyaluronan/config.toml -> [skin.hyaluronan])
  BDM_ASSIGN_CONFIG_VALUE(hyaluronan.enabled, "skin.hyaluronan.enabled");
  BDM_ASSIGN_CONFIG_VALUE(hyaluronan.diffusion, "skin.hyaluronan.diffusion");
  BDM_ASSIGN_CONFIG_VALUE(hyaluronan.decay, "skin.hyaluronan.decay");
  BDM_ASSIGN_CONFIG_VALUE(hyaluronan.basal_density, "skin.hyaluronan.basal_density");
  BDM_ASSIGN_CONFIG_VALUE(hyaluronan.reticular_density, "skin.hyaluronan.reticular_density");
  BDM_ASSIGN_CONFIG_VALUE(hyaluronan.production_rate, "skin.hyaluronan.production_rate");
  BDM_ASSIGN_CONFIG_VALUE(hyaluronan.water_retention_factor, "skin.hyaluronan.water_retention_factor");
  BDM_ASSIGN_CONFIG_VALUE(hyaluronan.migration_scaffold_factor, "skin.hyaluronan.migration_scaffold_factor");

  // Dermis module (modules/dermis/config.toml -> [skin.dermis])
  BDM_ASSIGN_CONFIG_VALUE(dermis.enabled, "skin.dermis.enabled");
  BDM_ASSIGN_CONFIG_VALUE(dermis.diffusion, "skin.dermis.diffusion");
  BDM_ASSIGN_CONFIG_VALUE(dermis.decay, "skin.dermis.decay");
  BDM_ASSIGN_CONFIG_VALUE(dermis.papillary_density, "skin.dermis.papillary_density");
  BDM_ASSIGN_CONFIG_VALUE(dermis.reticular_density, "skin.dermis.reticular_density");
  BDM_ASSIGN_CONFIG_VALUE(dermis.hypodermis_density, "skin.dermis.hypodermis_density");
  BDM_ASSIGN_CONFIG_VALUE(dermis.collagen_threshold, "skin.dermis.collagen_threshold");
  BDM_ASSIGN_CONFIG_VALUE(dermis.collagen_recovery_rate, "skin.dermis.collagen_recovery_rate");
  BDM_ASSIGN_CONFIG_VALUE(dermis.mmp_degradation, "skin.dermis.mmp_degradation");
  BDM_ASSIGN_CONFIG_VALUE(dermis.papillary_rate_factor, "skin.dermis.papillary_rate_factor");
  BDM_ASSIGN_CONFIG_VALUE(dermis.reticular_rate_factor, "skin.dermis.reticular_rate_factor");
  BDM_ASSIGN_CONFIG_VALUE(dermis.hypodermis_rate_factor, "skin.dermis.hypodermis_rate_factor");

  // pH module (modules/ph/config.toml -> [skin.ph])
  BDM_ASSIGN_CONFIG_VALUE(ph.recovery_rate, "skin.ph.recovery_rate");
  BDM_ASSIGN_CONFIG_VALUE(ph.migration_suppression, "skin.ph.migration_suppression");
  BDM_ASSIGN_CONFIG_VALUE(ph.mmp_boost, "skin.ph.mmp_boost");
  BDM_ASSIGN_CONFIG_VALUE(ph.biofilm_boost, "skin.ph.biofilm_boost");
  BDM_ASSIGN_CONFIG_VALUE(ph.bohr_factor, "skin.ph.bohr_factor");

  // Hemostasis / fibrin clot (modules/hemostasis/config.toml -> [skin.hemostasis])
  BDM_ASSIGN_CONFIG_VALUE(hemostasis.enabled, "skin.hemostasis.enabled");
  BDM_ASSIGN_CONFIG_VALUE(hemostasis.decay, "skin.hemostasis.decay");
  BDM_ASSIGN_CONFIG_VALUE(hemostasis.wound_seed, "skin.hemostasis.wound_seed");
  BDM_ASSIGN_CONFIG_VALUE(hemostasis.tgfb_coupling, "skin.hemostasis.tgfb_coupling");
  BDM_ASSIGN_CONFIG_VALUE(hemostasis.fibronectin_coupling, "skin.hemostasis.fibronectin_coupling");
  BDM_ASSIGN_CONFIG_VALUE(hemostasis.mmp_degradation, "skin.hemostasis.mmp_degradation");
  BDM_ASSIGN_CONFIG_VALUE(hemostasis.migration_boost, "skin.hemostasis.migration_boost");

  // Scab module (modules/scab/config.toml -> [skin.scab])
  BDM_ASSIGN_CONFIG_VALUE(scab.enabled, "skin.scab.enabled");
  BDM_ASSIGN_CONFIG_VALUE(scab.wound_seed, "skin.scab.wound_seed");
  BDM_ASSIGN_CONFIG_VALUE(scab.decay, "skin.scab.decay");
  BDM_ASSIGN_CONFIG_VALUE(scab.mmp_degradation, "skin.scab.mmp_degradation");
  BDM_ASSIGN_CONFIG_VALUE(scab.reepith_rate, "skin.scab.reepith_rate");
  BDM_ASSIGN_CONFIG_VALUE(scab.moisture_softening, "skin.scab.moisture_softening");
  BDM_ASSIGN_CONFIG_VALUE(scab.evaporation_shield, "skin.scab.evaporation_shield");
  BDM_ASSIGN_CONFIG_VALUE(scab.migration_penalty, "skin.scab.migration_penalty");

  // Blood module (modules/blood/config.toml -> [skin.blood])
  BDM_ASSIGN_CONFIG_VALUE(blood.enabled, "skin.blood.enabled");
  BDM_ASSIGN_CONFIG_VALUE(blood.bleed_rate, "skin.blood.bleed_rate");
  BDM_ASSIGN_CONFIG_VALUE(blood.depth_bleed_factor, "skin.blood.depth_bleed_factor");
  BDM_ASSIGN_CONFIG_VALUE(blood.vascularity_coupling, "skin.blood.vascularity_coupling");
  BDM_ASSIGN_CONFIG_VALUE(blood.intrinsic_rate, "skin.blood.intrinsic_rate");
  BDM_ASSIGN_CONFIG_VALUE(blood.extrinsic_rate, "skin.blood.extrinsic_rate");
  BDM_ASSIGN_CONFIG_VALUE(blood.thrombin_rate, "skin.blood.thrombin_rate");
  BDM_ASSIGN_CONFIG_VALUE(blood.thrombin_fibrin_coupling, "skin.blood.thrombin_fibrin_coupling");
  BDM_ASSIGN_CONFIG_VALUE(blood.platelet_aggregation, "skin.blood.platelet_aggregation_rate");
  BDM_ASSIGN_CONFIG_VALUE(blood.platelet_cytokine, "skin.blood.platelet_cytokine_coupling");
  BDM_ASSIGN_CONFIG_VALUE(blood.clot_maturation_h, "skin.blood.clot_maturation_h");
  BDM_ASSIGN_CONFIG_VALUE(blood.fibrinolysis_delay_h, "skin.blood.fibrinolysis_delay_h");
  BDM_ASSIGN_CONFIG_VALUE(blood.initial_volume, "skin.blood.initial_volume");
  BDM_ASSIGN_CONFIG_VALUE(blood.volume_loss_rate, "skin.blood.volume_loss_rate");
  BDM_ASSIGN_CONFIG_VALUE(blood.volume_recovery_rate, "skin.blood.volume_recovery_rate");
  BDM_ASSIGN_CONFIG_VALUE(blood.shock_threshold, "skin.blood.shock_threshold");
  BDM_ASSIGN_CONFIG_VALUE(blood.shock_perfusion_penalty, "skin.blood.shock_perfusion_penalty");
  BDM_ASSIGN_CONFIG_VALUE(blood.hemoglobin, "skin.blood.hemoglobin");
  BDM_ASSIGN_CONFIG_VALUE(blood.anemia_threshold, "skin.blood.anemia_threshold");
  BDM_ASSIGN_CONFIG_VALUE(blood.anemia_o2_penalty, "skin.blood.anemia_o2_penalty");
  BDM_ASSIGN_CONFIG_VALUE(blood.anticoagulant_level, "skin.blood.anticoagulant_level");
  BDM_ASSIGN_CONFIG_VALUE(blood.anticoag_coag_penalty, "skin.blood.anticoag_coag_penalty");
  BDM_ASSIGN_CONFIG_VALUE(blood.anticoag_bleed_boost, "skin.blood.anticoag_bleed_boost");
  BDM_ASSIGN_CONFIG_VALUE(blood.platelet_count, "skin.blood.platelet_count");
  BDM_ASSIGN_CONFIG_VALUE(blood.thrombocytopenia_threshold, "skin.blood.thrombocytopenia_threshold");

  // Burn module (modules/burn/config.toml -> [skin.burn])
  BDM_ASSIGN_CONFIG_VALUE(burn.enabled, "skin.burn.enabled");
  BDM_ASSIGN_CONFIG_VALUE(burn.coagulation_necrosis, "skin.burn.coagulation_necrosis");
  BDM_ASSIGN_CONFIG_VALUE(burn.stasis_viability, "skin.burn.stasis_initial_viability");
  BDM_ASSIGN_CONFIG_VALUE(burn.stasis_deterioration, "skin.burn.stasis_deterioration_rate");
  BDM_ASSIGN_CONFIG_VALUE(burn.hyperemia_boost, "skin.burn.hyperemia_inflammation_boost");
  BDM_ASSIGN_CONFIG_VALUE(burn.depth_fraction, "skin.burn.depth_fraction");
  BDM_ASSIGN_CONFIG_VALUE(burn.eschar_rate, "skin.burn.eschar_rate");
  BDM_ASSIGN_CONFIG_VALUE(burn.tewl_multiplier, "skin.burn.tewl_multiplier");
  BDM_ASSIGN_CONFIG_VALUE(burn.fluid_loss_rate, "skin.burn.fluid_loss_rate");
  BDM_ASSIGN_CONFIG_VALUE(burn.infection_susceptibility, "skin.burn.infection_susceptibility");
  BDM_ASSIGN_CONFIG_VALUE(burn.contracture_rate, "skin.burn.contracture_rate");

  // Pressure ulcer module (modules/pressure/config.toml -> [skin.pressure])
  BDM_ASSIGN_CONFIG_VALUE(pressure.enabled, "skin.pressure.enabled");
  BDM_ASSIGN_CONFIG_VALUE(pressure.compression_threshold, "skin.pressure.compression_threshold");
  BDM_ASSIGN_CONFIG_VALUE(pressure.ischemia_rate, "skin.pressure.ischemia_rate");
  BDM_ASSIGN_CONFIG_VALUE(pressure.reperfusion_ros_burst, "skin.pressure.reperfusion_ros_burst");
  BDM_ASSIGN_CONFIG_VALUE(pressure.shear_factor, "skin.pressure.shear_factor");
  BDM_ASSIGN_CONFIG_VALUE(pressure.tissue_damage_rate, "skin.pressure.tissue_damage_rate");
  BDM_ASSIGN_CONFIG_VALUE(pressure.reposition_interval_h, "skin.pressure.reposition_interval_h");
  BDM_ASSIGN_CONFIG_VALUE(pressure.moisture_damage_rate, "skin.pressure.moisture_damage_rate");
  BDM_ASSIGN_CONFIG_VALUE(pressure.stage_1, "skin.pressure.stage_1_threshold");
  BDM_ASSIGN_CONFIG_VALUE(pressure.stage_2, "skin.pressure.stage_2_threshold");
  BDM_ASSIGN_CONFIG_VALUE(pressure.stage_3, "skin.pressure.stage_3_threshold");
  BDM_ASSIGN_CONFIG_VALUE(pressure.stage_4, "skin.pressure.stage_4_threshold");

  // Photon transport module (modules/photon/config.toml -> [skin.photon])
  BDM_ASSIGN_CONFIG_VALUE(photon.enabled, "skin.photon.enabled");
  BDM_ASSIGN_CONFIG_VALUE(photon.absorption_coeff, "skin.photon.absorption_coeff");
  BDM_ASSIGN_CONFIG_VALUE(photon.scattering_coeff, "skin.photon.scattering_coeff");
  BDM_ASSIGN_CONFIG_VALUE(photon.anisotropy, "skin.photon.anisotropy");
  BDM_ASSIGN_CONFIG_VALUE(photon.irradiance, "skin.photon.irradiance");
  BDM_ASSIGN_CONFIG_VALUE(photon.beam_radius, "skin.photon.beam_radius");
  BDM_ASSIGN_CONFIG_VALUE(photon.beam_center_x, "skin.photon.beam_center_x");
  BDM_ASSIGN_CONFIG_VALUE(photon.beam_center_y, "skin.photon.beam_center_y");
  BDM_ASSIGN_CONFIG_VALUE(photon.diffusion, "skin.photon.diffusion");
  BDM_ASSIGN_CONFIG_VALUE(photon.decay, "skin.photon.decay");
  BDM_ASSIGN_CONFIG_VALUE(photon.opsin_activation_rate, "skin.photon.opsin_activation_rate");
  BDM_ASSIGN_CONFIG_VALUE(photon.opsin_deactivation_rate, "skin.photon.opsin_deactivation_rate");
  BDM_ASSIGN_CONFIG_VALUE(photon.opsin_saturation, "skin.photon.opsin_saturation");
  BDM_ASSIGN_CONFIG_VALUE(photon.phototoxicity_threshold, "skin.photon.phototoxicity_threshold");
  BDM_ASSIGN_CONFIG_VALUE(photon.phototoxicity_rate, "skin.photon.phototoxicity_rate");
  BDM_ASSIGN_CONFIG_VALUE(photon.thermal_coupling, "skin.photon.thermal_coupling");

  // Temperature module (modules/temperature/config.toml -> [skin.temperature])
  BDM_ASSIGN_CONFIG_VALUE(temperature.enabled, "skin.temperature.enabled");
  BDM_ASSIGN_CONFIG_VALUE(temperature.diffusion, "skin.temperature.diffusion");
  BDM_ASSIGN_CONFIG_VALUE(temperature.decay, "skin.temperature.decay");
  BDM_ASSIGN_CONFIG_VALUE(temperature.wound_surface, "skin.temperature.wound_surface");
  BDM_ASSIGN_CONFIG_VALUE(temperature.perfusion_warming, "skin.temperature.perfusion_warming_rate");
  BDM_ASSIGN_CONFIG_VALUE(temperature.surface_cooling, "skin.temperature.surface_cooling_rate");
  BDM_ASSIGN_CONFIG_VALUE(temperature.q10_migration, "skin.temperature.q10_migration");
  BDM_ASSIGN_CONFIG_VALUE(temperature.q10_proliferation, "skin.temperature.q10_proliferation");
  BDM_ASSIGN_CONFIG_VALUE(temperature.q10_mmp, "skin.temperature.q10_mmp");
  BDM_ASSIGN_CONFIG_VALUE(temperature.q10_biofilm, "skin.temperature.q10_biofilm");

  // Glucose module (modules/glucose/config.toml -> [skin.glucose])
  BDM_ASSIGN_CONFIG_VALUE(glucose_mod.enabled, "skin.glucose.enabled");
  BDM_ASSIGN_CONFIG_VALUE(glucose_mod.diffusion, "skin.glucose.diffusion");
  BDM_ASSIGN_CONFIG_VALUE(glucose_mod.decay, "skin.glucose.decay");
  BDM_ASSIGN_CONFIG_VALUE(glucose_mod.basal_conc, "skin.glucose.basal_conc");
  BDM_ASSIGN_CONFIG_VALUE(glucose_mod.perfusion_supply, "skin.glucose.perfusion_supply");
  BDM_ASSIGN_CONFIG_VALUE(glucose_mod.age_rate, "skin.glucose.age_rate");
  BDM_ASSIGN_CONFIG_VALUE(glucose_mod.age_inflammation, "skin.glucose.age_inflammation");
  BDM_ASSIGN_CONFIG_VALUE(glucose_mod.bacterial_consumption, "skin.glucose.bacterial_consumption");
  BDM_ASSIGN_CONFIG_VALUE(glucose_mod.prolif_threshold, "skin.glucose.prolif_threshold");
  // AGE (Advanced Glycation End-product) dynamics
  BDM_ASSIGN_CONFIG_VALUE(age_decay, "skin.glucose.age_decay");
  BDM_ASSIGN_CONFIG_VALUE(age_rage_inflammation, "skin.glucose.age_rage_inflammation");
  BDM_ASSIGN_CONFIG_VALUE(age_collagen_crosslink, "skin.glucose.age_collagen_crosslink");
  BDM_ASSIGN_CONFIG_VALUE(age_m1_prolongation, "skin.glucose.age_m1_prolongation");
  BDM_ASSIGN_CONFIG_VALUE(age_fibroblast_migration_impair, "skin.glucose.age_fibroblast_migration_impair");

  // Lactate module (modules/lactate/config.toml -> [skin.lactate])
  BDM_ASSIGN_CONFIG_VALUE(lactate.enabled, "skin.lactate.enabled");
  BDM_ASSIGN_CONFIG_VALUE(lactate.diffusion, "skin.lactate.diffusion");
  BDM_ASSIGN_CONFIG_VALUE(lactate.decay, "skin.lactate.decay");
  BDM_ASSIGN_CONFIG_VALUE(lactate.production_rate, "skin.lactate.production_rate");
  BDM_ASSIGN_CONFIG_VALUE(lactate.o2_threshold, "skin.lactate.o2_threshold");
  BDM_ASSIGN_CONFIG_VALUE(lactate.vegf_boost, "skin.lactate.vegf_boost");
  BDM_ASSIGN_CONFIG_VALUE(lactate.collagen_boost, "skin.lactate.collagen_boost");
  BDM_ASSIGN_CONFIG_VALUE(lactate.perfusion_clearance, "skin.lactate.perfusion_clearance");

  // Nitric oxide module (modules/nitric_oxide/config.toml -> [skin.nitric_oxide])
  BDM_ASSIGN_CONFIG_VALUE(nitric_oxide.enabled, "skin.nitric_oxide.enabled");
  BDM_ASSIGN_CONFIG_VALUE(no_diffusion, "skin.nitric_oxide.diffusion");
  BDM_ASSIGN_CONFIG_VALUE(no_decay, "skin.nitric_oxide.decay");
  BDM_ASSIGN_CONFIG_VALUE(no_m1_production, "skin.nitric_oxide.m1_production");
  BDM_ASSIGN_CONFIG_VALUE(no_neutrophil_production, "skin.nitric_oxide.neutrophil_production");
  BDM_ASSIGN_CONFIG_VALUE(no_vasodilation_factor, "skin.nitric_oxide.vasodilation_factor");
  BDM_ASSIGN_CONFIG_VALUE(no_antimicrobial_factor, "skin.nitric_oxide.antimicrobial_factor");
  BDM_ASSIGN_CONFIG_VALUE(no_collagen_suppression, "skin.nitric_oxide.collagen_suppression");
  BDM_ASSIGN_CONFIG_VALUE(diabetic.no_factor, "skin.diabetic.no_factor");

  // Senescence module (modules/senescence/config.toml -> [skin.senescence])
  BDM_ASSIGN_CONFIG_VALUE(senescence.enabled, "skin.senescence.enabled");
  BDM_ASSIGN_CONFIG_VALUE(senescence.diffusion, "skin.senescence.diffusion");
  BDM_ASSIGN_CONFIG_VALUE(senescence.decay, "skin.senescence.decay");
  BDM_ASSIGN_CONFIG_VALUE(senescence.wound_rate, "skin.senescence.wound_accumulation_rate");
  BDM_ASSIGN_CONFIG_VALUE(senescence.infl_rate, "skin.senescence.inflammation_accumulation");
  BDM_ASSIGN_CONFIG_VALUE(senescence.age_rate, "skin.senescence.age_accumulation");
  BDM_ASSIGN_CONFIG_VALUE(sasp_inflammation_rate, "skin.senescence.sasp_inflammation_rate");
  BDM_ASSIGN_CONFIG_VALUE(sasp_mmp_rate, "skin.senescence.sasp_mmp_rate");
  BDM_ASSIGN_CONFIG_VALUE(sasp_tgfb_rate, "skin.senescence.sasp_tgfb_rate");
  BDM_ASSIGN_CONFIG_VALUE(senolytic_clearance_rate, "skin.senescence.senolytic_clearance_rate");
  BDM_ASSIGN_CONFIG_VALUE(senescence.immune_clearance_rate, "skin.senescence.immune_clearance_rate");
  BDM_ASSIGN_CONFIG_VALUE(diabetic.senescence_factor, "skin.senescence.diabetic_accumulation_factor");

  // Neuropathy module (modules/neuropathy/config.toml -> [skin.neuropathy])
  BDM_ASSIGN_CONFIG_VALUE(neuropathy.enabled, "skin.neuropathy.enabled");
  BDM_ASSIGN_CONFIG_VALUE(neuropathy.diffusion, "skin.neuropathy.diffusion");
  BDM_ASSIGN_CONFIG_VALUE(neuropathy.decay, "skin.neuropathy.decay");
  BDM_ASSIGN_CONFIG_VALUE(neuropathy.basal_density, "skin.neuropathy.basal_nerve_density");
  BDM_ASSIGN_CONFIG_VALUE(nerve_regeneration_rate, "skin.neuropathy.regeneration_rate");
  BDM_ASSIGN_CONFIG_VALUE(nerve_vegf_boost, "skin.neuropathy.regeneration_vegf_boost");
  BDM_ASSIGN_CONFIG_VALUE(nerve_tgfb_inhibit, "skin.neuropathy.regeneration_tgfb_inhibit");
  BDM_ASSIGN_CONFIG_VALUE(substance_p_inflammation, "skin.neuropathy.substance_p_inflammation");
  BDM_ASSIGN_CONFIG_VALUE(substance_p_proliferation, "skin.neuropathy.substance_p_proliferation");
  BDM_ASSIGN_CONFIG_VALUE(cgrp_vasodilation, "skin.neuropathy.cgrp_vasodilation");
  BDM_ASSIGN_CONFIG_VALUE(diabetic.nerve_factor, "skin.neuropathy.diabetic_nerve_factor");
  BDM_ASSIGN_CONFIG_VALUE(diabetic.nerve_regen_factor, "skin.neuropathy.diabetic_regeneration_factor");

  // ROS module (modules/ros/config.toml -> [skin.ros])
  BDM_ASSIGN_CONFIG_VALUE(ros.enabled, "skin.ros.enabled");
  BDM_ASSIGN_CONFIG_VALUE(ros.diffusion, "skin.ros.diffusion");
  BDM_ASSIGN_CONFIG_VALUE(ros.decay, "skin.ros.decay");
  BDM_ASSIGN_CONFIG_VALUE(ros.neutrophil_burst, "skin.ros.neutrophil_burst");
  BDM_ASSIGN_CONFIG_VALUE(ros.m1_burst, "skin.ros.m1_burst");
  BDM_ASSIGN_CONFIG_VALUE(ros.mitochondrial_rate, "skin.ros.mitochondrial_rate");
  BDM_ASSIGN_CONFIG_VALUE(ros.hypoxia_threshold, "skin.ros.hypoxia_threshold");
  BDM_ASSIGN_CONFIG_VALUE(ros.hypoxia_amplification, "skin.ros.hypoxia_amplification");
  BDM_ASSIGN_CONFIG_VALUE(diabetic.mito_ros_factor, "skin.ros.diabetic_mito_factor");
  BDM_ASSIGN_CONFIG_VALUE(diabetic.antioxidant_factor, "skin.ros.diabetic_antioxidant_factor");
  BDM_ASSIGN_CONFIG_VALUE(ros.senescence_rate, "skin.ros.senescence_ros_rate");
  BDM_ASSIGN_CONFIG_VALUE(ros.mmp_activation, "skin.ros.mmp_activation_rate");
  BDM_ASSIGN_CONFIG_VALUE(ros.angiogenesis_impairment, "skin.ros.angiogenesis_impairment");
  BDM_ASSIGN_CONFIG_VALUE(ros.inflammation_amplification, "skin.ros.inflammation_amplification");
  BDM_ASSIGN_CONFIG_VALUE(ros.collagen_damage, "skin.ros.collagen_damage_rate");
  BDM_ASSIGN_CONFIG_VALUE(ros.tissue_antioxidant, "skin.ros.tissue_antioxidant_rate");
  BDM_ASSIGN_CONFIG_VALUE(ros.perfusion_clearance, "skin.ros.perfusion_clearance");

  // Mechanotransduction module (modules/mechanotransduction/config.toml -> [skin.mechanotransduction])
  BDM_ASSIGN_CONFIG_VALUE(mechanotransduction.enabled, "skin.mechanotransduction.enabled");
  BDM_ASSIGN_CONFIG_VALUE(stiffness_yap_threshold, "skin.mechanotransduction.stiffness_yap_threshold");
  BDM_ASSIGN_CONFIG_VALUE(stiffness_contraction_rate, "skin.mechanotransduction.stiffness_contraction_rate");
  BDM_ASSIGN_CONFIG_VALUE(stiffness_scar_factor, "skin.mechanotransduction.stiffness_scar_factor");

  // Lymphatic module (modules/lymphatic/config.toml -> [skin.lymphatic])
  BDM_ASSIGN_CONFIG_VALUE(lymphatic.enabled, "skin.lymphatic.enabled");
  BDM_ASSIGN_CONFIG_VALUE(lymphatic.diffusion, "skin.lymphatic.diffusion");
  BDM_ASSIGN_CONFIG_VALUE(lymphatic.basal_density, "skin.lymphatic.basal_density");
  BDM_ASSIGN_CONFIG_VALUE(lymphatic.regen_rate, "skin.lymphatic.regen_rate");
  BDM_ASSIGN_CONFIG_VALUE(lymphatic.vegf_boost, "skin.lymphatic.vegf_boost");
  BDM_ASSIGN_CONFIG_VALUE(edema_leak_rate, "skin.lymphatic.edema_leak_rate");
  BDM_ASSIGN_CONFIG_VALUE(edema_drainage_rate, "skin.lymphatic.edema_drainage_rate");
  BDM_ASSIGN_CONFIG_VALUE(edema_o2_impairment, "skin.lymphatic.edema_o2_impairment");
  BDM_ASSIGN_CONFIG_VALUE(edema_migration_impairment, "skin.lymphatic.edema_migration_impairment");
  BDM_ASSIGN_CONFIG_VALUE(diabetic.lymphatic_factor, "skin.lymphatic.diabetic_lymphatic_factor");

  // Bioelectric module (modules/bioelectric/config.toml -> [skin.bioelectric])
  BDM_ASSIGN_CONFIG_VALUE(bioelectric.enabled, "skin.bioelectric.enabled");
  BDM_ASSIGN_CONFIG_VALUE(voltage_diffusion, "skin.bioelectric.voltage_diffusion");
  BDM_ASSIGN_CONFIG_VALUE(voltage_decay, "skin.bioelectric.voltage_decay");
  BDM_ASSIGN_CONFIG_VALUE(voltage_epithelial_source, "skin.bioelectric.voltage_epithelial_source");
  BDM_ASSIGN_CONFIG_VALUE(galvanotaxis_strength, "skin.bioelectric.galvanotaxis_strength");
  BDM_ASSIGN_CONFIG_VALUE(diabetic.voltage_factor, "skin.bioelectric.diabetic_voltage_factor");

  // Rheumatoid arthritis module (studies/rheumatoid/modules/rheumatoid/config.toml -> [skin.ra])
  BDM_ASSIGN_CONFIG_VALUE(ra.enabled, "skin.ra.enabled");
  BDM_ASSIGN_CONFIG_VALUE(tnf_alpha_diffusion, "skin.ra.tnf_diffusion");
  BDM_ASSIGN_CONFIG_VALUE(tnf_alpha_decay, "skin.ra.tnf_decay");
  BDM_ASSIGN_CONFIG_VALUE(il6_diffusion, "skin.ra.il6_diffusion");
  BDM_ASSIGN_CONFIG_VALUE(il6_decay, "skin.ra.il6_decay");
  BDM_ASSIGN_CONFIG_VALUE(ra.autoimmune_source, "skin.ra.autoimmune_source");
  BDM_ASSIGN_CONFIG_VALUE(ra.il6_autoimmune_source, "skin.ra.il6_autoimmune_source");
  BDM_ASSIGN_CONFIG_VALUE(ra.tnf_m1_rate, "skin.ra.tnf_m1_rate");
  BDM_ASSIGN_CONFIG_VALUE(ra.tnf_neutrophil_rate, "skin.ra.tnf_neutrophil_rate");
  BDM_ASSIGN_CONFIG_VALUE(ra.il6_m1_rate, "skin.ra.il6_m1_rate");
  BDM_ASSIGN_CONFIG_VALUE(ra.flare_delay, "skin.ra.flare_delay");
  BDM_ASSIGN_CONFIG_VALUE(ra.flare_steepness, "skin.ra.flare_steepness");
  BDM_ASSIGN_CONFIG_VALUE(ra.tnf_inflammation_coupling, "skin.ra.tnf_inflammation_coupling");
  BDM_ASSIGN_CONFIG_VALUE(ra.tnf_mmp_boost, "skin.ra.tnf_mmp_boost");
  BDM_ASSIGN_CONFIG_VALUE(ra.tnf_vegf_boost, "skin.ra.tnf_vegf_boost");
  BDM_ASSIGN_CONFIG_VALUE(ra.tnf_il6_induction, "skin.ra.tnf_il6_induction");
  BDM_ASSIGN_CONFIG_VALUE(ra.il6_inflammation_coupling, "skin.ra.il6_inflammation_coupling");
  BDM_ASSIGN_CONFIG_VALUE(ra.il6_cartilage_boost, "skin.ra.il6_cartilage_boost");
  BDM_ASSIGN_CONFIG_VALUE(ra.cartilage_basal, "skin.ra.cartilage_basal");
  BDM_ASSIGN_CONFIG_VALUE(ra.cartilage_mmp_degradation, "skin.ra.cartilage_mmp_degradation");
  BDM_ASSIGN_CONFIG_VALUE(ra.cartilage_tnf_degradation, "skin.ra.cartilage_tnf_degradation");
  BDM_ASSIGN_CONFIG_VALUE(ra.pannus_fibroblast_boost, "skin.ra.pannus_fibroblast_boost");
  BDM_ASSIGN_CONFIG_VALUE(ra.m1_prolongation, "skin.ra.m1_prolongation");
  BDM_ASSIGN_CONFIG_VALUE(ra.synovial_basal, "skin.ra.synovial_basal");
  BDM_ASSIGN_CONFIG_VALUE(ra.synovial_growth_rate, "skin.ra.synovial_growth_rate");
  BDM_ASSIGN_CONFIG_VALUE(ra.synovial_carrying_capacity, "skin.ra.synovial_carrying_capacity");
  BDM_ASSIGN_CONFIG_VALUE(ra.synovial_cytokine_boost, "skin.ra.synovial_cytokine_boost");
  BDM_ASSIGN_CONFIG_VALUE(ra.synovial_erosion_boost, "skin.ra.synovial_erosion_boost");
  BDM_ASSIGN_CONFIG_VALUE(ra.tcell_diffusion, "skin.ra.tcell_diffusion");
  BDM_ASSIGN_CONFIG_VALUE(ra.tcell_decay, "skin.ra.tcell_decay");
  BDM_ASSIGN_CONFIG_VALUE(ra.tcell_recruitment, "skin.ra.tcell_recruitment");
  BDM_ASSIGN_CONFIG_VALUE(ra.tcell_proliferation, "skin.ra.tcell_proliferation");
  BDM_ASSIGN_CONFIG_VALUE(ra.tcell_carrying_capacity, "skin.ra.tcell_carrying_capacity");
  BDM_ASSIGN_CONFIG_VALUE(ra.tcell_cytokine_weight, "skin.ra.tcell_cytokine_weight");
  BDM_ASSIGN_CONFIG_VALUE(ra.bone_basal, "skin.ra.bone_basal");
  BDM_ASSIGN_CONFIG_VALUE(ra.bone_rankl_erosion, "skin.ra.bone_rankl_erosion");
  BDM_ASSIGN_CONFIG_VALUE(ra.bone_tnf_erosion, "skin.ra.bone_tnf_erosion");
  BDM_ASSIGN_CONFIG_VALUE(ra.bone_depth_z, "skin.ra.bone_depth_z");
  BDM_ASSIGN_HOURS(ra.treatment_delay, "skin.ra.treatment_delay_h");
  BDM_ASSIGN_HOURS(ra.drug_onset, "skin.ra.drug_onset_h");
  BDM_ASSIGN_CONFIG_VALUE(ra.anti_tnf_clearance, "skin.ra.anti_tnf_clearance");
  BDM_ASSIGN_CONFIG_VALUE(ra.anti_il6r_clearance, "skin.ra.anti_il6r_clearance");

  // Scar formation (modules/scar.toml -> [skin.scar])
  BDM_ASSIGN_CONFIG_VALUE(scar_enabled, "skin.scar.enabled");
  BDM_ASSIGN_CONFIG_VALUE(scar.proportional_enabled, "skin.scar.proportional_enabled");
  BDM_ASSIGN_CONFIG_VALUE(scar.accumulation_rate, "skin.scar.accumulation_rate");
  BDM_ASSIGN_CONFIG_VALUE(scar.collagen_threshold, "skin.scar.collagen_threshold");
  BDM_ASSIGN_CONFIG_VALUE(scar.maturity_enabled, "skin.scar.maturity_enabled");
  BDM_ASSIGN_CONFIG_VALUE(scar.maturation_rate, "skin.scar.maturation_rate");
  BDM_ASSIGN_CONFIG_VALUE(scar.maturation_mmp_boost, "skin.scar.maturation_mmp_boost");
  BDM_ASSIGN_CONFIG_VALUE(scar.maturation_infl_block, "skin.scar.maturation_infl_block");
  BDM_ASSIGN_CONFIG_VALUE(scar.maturation_myofib_block, "skin.scar.maturation_myofib_block");

  // Metrics export (TOML: hours)
  BDM_ASSIGN_HOURS(metrics_interval, "skin.metrics_interval_h");
  BDM_ASSIGN_CONFIG_VALUE(metrics_autoopen, "skin.metrics_autoopen");
  BDM_ASSIGN_CONFIG_VALUE(headless, "skin.headless");

  // Cytokine taper rates
  BDM_ASSIGN_CONFIG_VALUE(immune.neutrophil_cytokine_taper, "skin.immune.neutrophil_cytokine_taper");
  BDM_ASSIGN_CONFIG_VALUE(m1_cytokine_taper, "skin.immune.m1_cytokine_taper");
  BDM_ASSIGN_CONFIG_VALUE(m2_tgfb_taper, "skin.immune.m2_tgfb_taper");
  BDM_ASSIGN_CONFIG_VALUE(fibroblast.tgfb_taper, "skin.fibroblast.tgfb_taper_rate");

  // PDE sub-cycling
  BDM_ASSIGN_CONFIG_VALUE(subcycle_slow, "skin.subcycle_slow");
  BDM_ASSIGN_CONFIG_VALUE(subcycle_medium, "skin.subcycle_medium");

  // Agent behavior sub-cycling
  BDM_ASSIGN_CONFIG_VALUE(derived_fields_subcycle, "skin.derived_fields_subcycle");
  BDM_ASSIGN_CONFIG_VALUE(migration_subcycle, "skin.migration_subcycle");
  BDM_ASSIGN_CONFIG_VALUE(homeostasis_subcycle, "skin.homeostasis_subcycle");

  // Basal density continuum
  BDM_ASSIGN_CONFIG_VALUE(basal_density_enabled, "skin.basal_density_enabled");
  BDM_ASSIGN_CONFIG_VALUE(basal_density_max, "skin.basal_density_max");
  BDM_ASSIGN_CONFIG_VALUE(basal_density_div_rate, "skin.basal_density_div_rate");
  BDM_ASSIGN_CONFIG_VALUE(basal_density_diff_rate, "skin.basal_density_diff_rate");
  BDM_ASSIGN_CONFIG_VALUE(basal_density_subcycle, "skin.basal_density_subcycle");
  BDM_ASSIGN_CONFIG_VALUE(demotion_radius_factor, "skin.demotion_radius_factor");

  // Hot-reload
  BDM_ASSIGN_CONFIG_VALUE(hot_reload, "skin.hot_reload");

  // Derived composite fields
  BDM_ASSIGN_CONFIG_VALUE(ecm_weight_collagen, "skin.derived.ecm_weight_collagen");
  BDM_ASSIGN_CONFIG_VALUE(ecm_weight_fibronectin, "skin.derived.ecm_weight_fibronectin");
  BDM_ASSIGN_CONFIG_VALUE(ecm_weight_elastin, "skin.derived.ecm_weight_elastin");
  BDM_ASSIGN_CONFIG_VALUE(ecm_weight_fibrin, "skin.derived.ecm_weight_fibrin");
  BDM_ASSIGN_CONFIG_VALUE(ecm_migration_boost, "skin.derived.ecm_migration_boost");

  // Multi-resolution structural fields
  BDM_ASSIGN_CONFIG_VALUE(grid_resolution_structural, "skin.derived.grid_resolution_structural");

  // Debug profile
  BDM_ASSIGN_CONFIG_VALUE(debug_immune, "skin.debug.immune");
  BDM_ASSIGN_CONFIG_VALUE(debug_fibroblast, "skin.debug.fibroblast");
  BDM_ASSIGN_CONFIG_VALUE(debug_wound, "skin.debug.wound");
  BDM_ASSIGN_CONFIG_VALUE(debug_scaled_grid, "skin.debug.scaled_grid");
  BDM_ASSIGN_CONFIG_VALUE(debug_perf, "skin.debug.perf");
}

void bdm::skibidy::SimParam::ValidateConfig() const {
  auto check = [](bool condition, const char* msg) {
    if (!condition) {
      Log::Fatal("SimParam::ValidateConfig", msg);
    }
  };

  // Geometry
  check(geo_patch_um > 0, "geometry.patch_um must be positive");
  check(geo_depth_um > 0, "geometry.depth_um must be positive");
  check(geo_height_um > 0, "geometry.height_um must be positive");
  check(geo_margin_um >= 0, "geometry.margin_um must be non-negative");
  check(geo_um_per_unit > 0, "geometry.um_per_unit must be positive");
  check(geo_voxel_um > 0, "geometry.voxel_um must be positive");
  check(tissue_max > tissue_min, "derived tissue_max must exceed tissue_min");
  check(domain_max > domain_min, "derived domain_max must exceed domain_min");

  // Wound geometry
  check(wound.radius > 0, "wound.radius must be positive");
  check(wound.center_x >= 0, "wound.center_x must be non-negative");
  check(wound.center_y >= 0, "wound.center_y must be non-negative");

  // Simulation grid
  check(grid_resolution >= 1, "grid_resolution must be >= 1");
  check(num_steps >= 1, "num_steps must be >= 1");

  // Decay rates must be non-negative
  check(inflammation.decay >= 0, "inflammation.decay must be non-negative");
  check(oxygen_decay >= 0, "oxygen_decay must be non-negative");
  check(water_decay >= 0, "water_decay must be non-negative");
  check(fibroblast.tgfb_decay >= 0, "fibroblast.tgfb_decay must be non-negative");
  check(fibroblast.collagen_decay >= 0, "fibroblast.collagen_decay must be non-negative");
  check(mmp.residual_decay >= 0, "mmp.residual_decay must be non-negative");
  check(mmp.timp_decay >= 0, "mmp.timp_decay must be non-negative");

  // Diffusion coefficients must be non-negative
  check(inflammation.diffusion >= 0, "inflammation.diffusion must be non-negative");
  check(oxygen_diffusion >= 0, "oxygen_diffusion must be non-negative");
  check(mmp.diffusion >= 0, "mmp.diffusion must be non-negative");
  check(mmp.timp_diffusion >= 0, "mmp.timp_diffusion must be non-negative");
  check(fibroblast.tgfb_diffusion >= 0, "fibroblast.tgfb_diffusion must be non-negative");
  check(angiogenesis.vegf_diffusion >= 0, "angiogenesis.vegf_diffusion must be non-negative");

  // Diabetic factors must be positive (multipliers)
  if (diabetic.mode) {
    check(diabetic.m1_duration_factor > 0, "diabetic.m1_duration_factor must be positive");
    check(diabetic.resolution_factor > 0, "diabetic.resolution_factor must be positive");
    check(diabetic.prolif_factor > 0, "diabetic.prolif_factor must be positive");
    check(diabetic.migration_factor > 0, "diabetic.migration_factor must be positive");
    check(diabetic.migration_infl_K > 0, "diabetic.migration_infl_K must be positive (Hill half-max)");
    check(diabetic.prolif_infl_K > 0, "diabetic.prolif_infl_K must be positive (Hill half-max)");
  }

  // Probabilities in [0,1]
  check(immune.neutrophil_apoptosis_rate >= 0 && immune.neutrophil_apoptosis_rate <= 1,
        "immune.neutrophil_apoptosis_rate must be in [0,1]");
  check(immune.macrophage_apoptosis_rate >= 0 && immune.macrophage_apoptosis_rate <= 1,
        "immune.macrophage_apoptosis_rate must be in [0,1]");
  check(immune.macrophage_emigration_rate >= 0 && immune.macrophage_emigration_rate <= 1,
        "immune.macrophage_emigration_rate must be in [0,1]");

  // Cell sizes
  check(division_diameter > 0, "division_diameter must be positive");
  check(immune_cell_diameter > 0, "immune_cell_diameter must be positive");

  // Perfusion
  check(perfusion.basal >= 0 && perfusion.basal <= 2.0,
        "perfusion.basal should be in [0,2] (normalized)");
  check(perfusion.angio_rate >= 0, "perfusion.angio_rate must be non-negative");
}
