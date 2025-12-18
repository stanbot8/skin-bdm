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

  // Volume boundaries
  BDM_ASSIGN_CONFIG_VALUE(volume_z_spinous, "skin.volume_z_spinous");
  BDM_ASSIGN_CONFIG_VALUE(volume_z_granular, "skin.volume_z_granular");
  BDM_ASSIGN_CONFIG_VALUE(volume_z_cornified, "skin.volume_z_cornified");

  // Dermal sub-layer boundaries
  BDM_ASSIGN_CONFIG_VALUE(dermal_z_papillary, "skin.dermal_z_papillary");
  BDM_ASSIGN_CONFIG_VALUE(dermal_z_reticular, "skin.dermal_z_reticular");

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

  // Tissue extent
  BDM_ASSIGN_CONFIG_VALUE(tissue_min, "skin.tissue_min");
  BDM_ASSIGN_CONFIG_VALUE(tissue_max, "skin.tissue_max");

  // Simulation (TOML: days)
  BDM_ASSIGN_DAYS(num_steps, "skin.duration_days");
  BDM_ASSIGN_CONFIG_VALUE(grid_resolution, "skin.grid_resolution");
  BDM_ASSIGN_CONFIG_VALUE(ref_box_length, "skin.ref_box_length");

  // Wound event (modules/wound.toml -> [skin.wound])
  BDM_ASSIGN_CONFIG_VALUE(wound_enabled, "skin.wound.enabled");
  BDM_ASSIGN_CONFIG_VALUE(wound_center_x, "skin.wound.center_x");
  BDM_ASSIGN_CONFIG_VALUE(wound_center_y, "skin.wound.center_y");
  BDM_ASSIGN_CONFIG_VALUE(wound_radius, "skin.wound.radius");
  BDM_ASSIGN_HOURS(wound_trigger_step, "skin.wound.trigger_h");
  BDM_ASSIGN_CONFIG_VALUE(wound_inward_bias, "skin.wound.inward_bias");
  BDM_ASSIGN_CONFIG_VALUE(wound_vascular_damage, "skin.wound.vascular_damage");

  // Vascular perfusion (modules/perfusion.toml -> [skin.perfusion])
  BDM_ASSIGN_CONFIG_VALUE(perfusion_diffusion, "skin.perfusion.diffusion");
  BDM_ASSIGN_CONFIG_VALUE(perfusion_decay, "skin.perfusion.decay");
  BDM_ASSIGN_CONFIG_VALUE(perfusion_basal, "skin.perfusion.basal");
  BDM_ASSIGN_CONFIG_VALUE(perfusion_angio_rate, "skin.perfusion.angio_rate");
  BDM_ASSIGN_HOURS(perfusion_angio_delay, "skin.perfusion.angio_delay_h");
  BDM_ASSIGN_CONFIG_VALUE(perfusion_papillary_fraction, "skin.perfusion.papillary_fraction");
  BDM_ASSIGN_CONFIG_VALUE(perfusion_reticular_fraction, "skin.perfusion.reticular_fraction");
  BDM_ASSIGN_CONFIG_VALUE(perfusion_hypodermis_fraction, "skin.perfusion.hypodermis_fraction");
  BDM_ASSIGN_CONFIG_VALUE(angio_papillary_factor, "skin.perfusion.angio_papillary_factor");
  BDM_ASSIGN_CONFIG_VALUE(angio_reticular_factor, "skin.perfusion.angio_reticular_factor");
  BDM_ASSIGN_CONFIG_VALUE(angio_hypodermis_factor, "skin.perfusion.angio_hypodermis_factor");

  // Dynamic field coupling
  BDM_ASSIGN_CONFIG_VALUE(calcium_recovery_rate, "skin.calcium_recovery_rate");
  BDM_ASSIGN_CONFIG_VALUE(oxygen_prolif_threshold, "skin.oxygen_prolif_threshold");
  BDM_ASSIGN_CONFIG_VALUE(oxygen_recovery_enabled, "skin.oxygen_recovery_enabled");

  // Cell migration
  BDM_ASSIGN_CONFIG_VALUE(migration_enabled, "skin.migration_enabled");
  BDM_ASSIGN_CONFIG_VALUE(migration_speed, "skin.migration_speed");

  // Per-cell UWYN handoff (TOML: hours)
  BDM_ASSIGN_HOURS(handoff_delay, "skin.handoff_delay_h");

  // Inflammation (modules/inflammation.toml -> [skin.inflammation])
  BDM_ASSIGN_CONFIG_VALUE(inflammation_diffusion, "skin.inflammation.diffusion");
  BDM_ASSIGN_CONFIG_VALUE(inflammation_decay, "skin.inflammation.decay");
  BDM_ASSIGN_CONFIG_VALUE(inflammation_migration_threshold, "skin.inflammation.migration_threshold");
  BDM_ASSIGN_CONFIG_VALUE(inflammation_prolif_threshold, "skin.inflammation.prolif_threshold");
  BDM_ASSIGN_CONFIG_VALUE(wound_inflammation_source_rate, "skin.inflammation.wound_source_rate");
  BDM_ASSIGN_CONFIG_VALUE(wound_inflammation_source_taper, "skin.inflammation.wound_source_taper");
  BDM_ASSIGN_CONFIG_VALUE(split_inflammation_enabled, "skin.inflammation.split_enabled");
  BDM_ASSIGN_CONFIG_VALUE(anti_inflammation_diffusion, "skin.inflammation.anti_diffusion");
  BDM_ASSIGN_CONFIG_VALUE(anti_inflammation_decay, "skin.inflammation.anti_decay");
  BDM_ASSIGN_CONFIG_VALUE(anti_inflammation_weight, "skin.inflammation.anti_weight");

  // Immune response (modules/immune.toml -> [skin.immune]; TOML: hours)
  BDM_ASSIGN_HOURS(neutrophil_spawn_delay, "skin.immune.neutrophil_spawn_delay_h");
  BDM_ASSIGN_CONFIG_VALUE(neutrophil_spawn_waves, "skin.immune.neutrophil_spawn_waves");
  BDM_ASSIGN_HOURS(neutrophil_spawn_window, "skin.immune.neutrophil_spawn_window_h");
  BDM_ASSIGN_HOURS(neutrophil_lifespan, "skin.immune.neutrophil_lifespan_h");
  BDM_ASSIGN_HOURS(neutrophil_min_survival, "skin.immune.neutrophil_min_survival_h");
  BDM_ASSIGN_CONFIG_VALUE(neutrophil_apoptosis_rate, "skin.immune.neutrophil_apoptosis_rate");
  BDM_ASSIGN_HOURS(macrophage_spawn_delay, "skin.immune.macrophage_spawn_delay_h");
  BDM_ASSIGN_CONFIG_VALUE(macrophage_spawn_threshold, "skin.immune.macrophage_spawn_threshold");
  BDM_ASSIGN_CONFIG_VALUE(macrophage_spawn_rate, "skin.immune.macrophage_spawn_rate");
  BDM_ASSIGN_CONFIG_VALUE(macrophage_spawn_taper, "skin.immune.macrophage_spawn_taper");
  BDM_ASSIGN_HOURS(macrophage_m1_duration, "skin.immune.macrophage_m1_duration_h");
  BDM_ASSIGN_CONFIG_VALUE(m1_transition_threshold, "skin.immune.m1_transition_threshold");
  BDM_ASSIGN_HOURS(m1_transition_min_age, "skin.immune.m1_transition_min_age_h");
  BDM_ASSIGN_HOURS(macrophage_lifespan, "skin.immune.macrophage_lifespan_h");
  BDM_ASSIGN_HOURS(macrophage_min_survival, "skin.immune.macrophage_min_survival_h");
  BDM_ASSIGN_CONFIG_VALUE(macrophage_apoptosis_rate, "skin.immune.macrophage_apoptosis_rate");
  BDM_ASSIGN_CONFIG_VALUE(macrophage_emigration_rate, "skin.immune.macrophage_emigration_rate");
  BDM_ASSIGN_CONFIG_VALUE(immune_cell_diameter, "skin.immune.cell_diameter");
  BDM_ASSIGN_CONFIG_VALUE(immune_cytokine_rate, "skin.immune.cytokine_rate");
  BDM_ASSIGN_CONFIG_VALUE(immune_resolution_rate, "skin.immune.resolution_rate");
  BDM_ASSIGN_CONFIG_VALUE(immune_migration_speed, "skin.immune.migration_speed");
  BDM_ASSIGN_CONFIG_VALUE(efferocytosis_enabled, "skin.immune.efferocytosis_enabled");
  BDM_ASSIGN_CONFIG_VALUE(efferocytosis_radius, "skin.immune.efferocytosis_radius");
  BDM_ASSIGN_CONFIG_VALUE(efferocytosis_age_fraction, "skin.immune.efferocytosis_age_fraction");
  BDM_ASSIGN_CONFIG_VALUE(chemotaxis_enabled, "skin.immune.chemotaxis_enabled");
  BDM_ASSIGN_CONFIG_VALUE(chemotaxis_speed_scale, "skin.immune.chemotaxis_speed_scale");

  // Biofilm dynamics (modules/biofilm.toml -> [skin.biofilm])
  BDM_ASSIGN_CONFIG_VALUE(biofilm_enabled, "skin.biofilm.enabled");
  BDM_ASSIGN_CONFIG_VALUE(biofilm_growth_rate, "skin.biofilm.growth_rate");
  BDM_ASSIGN_CONFIG_VALUE(biofilm_carrying_capacity, "skin.biofilm.carrying_capacity");
  BDM_ASSIGN_HOURS(biofilm_seed_delay, "skin.biofilm.seed_delay_h");
  BDM_ASSIGN_CONFIG_VALUE(biofilm_seed_amount, "skin.biofilm.seed_amount");
  BDM_ASSIGN_CONFIG_VALUE(biofilm_neutrophil_clearance, "skin.biofilm.neutrophil_clearance");
  BDM_ASSIGN_CONFIG_VALUE(biofilm_macrophage_clearance, "skin.biofilm.macrophage_clearance");
  BDM_ASSIGN_CONFIG_VALUE(biofilm_m1_block_threshold, "skin.biofilm.m1_block_threshold");
  BDM_ASSIGN_CONFIG_VALUE(biofilm_inflammation_rate, "skin.biofilm.inflammation_rate");

  // VEGF-driven angiogenesis (modules/angiogenesis.toml -> [skin.angiogenesis])
  BDM_ASSIGN_CONFIG_VALUE(angiogenesis_enabled, "skin.angiogenesis.enabled");
  BDM_ASSIGN_CONFIG_VALUE(vegf_diffusion, "skin.angiogenesis.vegf_diffusion");
  BDM_ASSIGN_CONFIG_VALUE(vegf_decay, "skin.angiogenesis.vegf_decay");
  BDM_ASSIGN_CONFIG_VALUE(vegf_production_rate, "skin.angiogenesis.vegf_production_rate");
  BDM_ASSIGN_CONFIG_VALUE(vegf_hypoxia_threshold, "skin.angiogenesis.vegf_hypoxia_threshold");
  BDM_ASSIGN_CONFIG_VALUE(vegf_consumption_rate, "skin.angiogenesis.vegf_consumption_rate");
  BDM_ASSIGN_CONFIG_VALUE(angio_vegf_rate, "skin.angiogenesis.angio_vegf_rate");
  BDM_ASSIGN_CONFIG_VALUE(m2_vegf_rate, "skin.angiogenesis.m2_vegf_rate");
  BDM_ASSIGN_CONFIG_VALUE(vegf_production_taper, "skin.angiogenesis.vegf_production_taper");
  BDM_ASSIGN_CONFIG_VALUE(vegf_receptor_clearance, "skin.angiogenesis.vegf_receptor_clearance");

  // Diabetic modifiers (modules/diabetic.toml -> [skin.diabetic])
  BDM_ASSIGN_CONFIG_VALUE(diabetic_mode, "skin.diabetic.mode");
  BDM_ASSIGN_CONFIG_VALUE(diabetic_m1_duration_factor, "skin.diabetic.m1_duration_factor");
  BDM_ASSIGN_CONFIG_VALUE(diabetic_resolution_factor, "skin.diabetic.resolution_factor");
  BDM_ASSIGN_CONFIG_VALUE(diabetic_efferocytosis_factor, "skin.diabetic.efferocytosis_factor");
  BDM_ASSIGN_CONFIG_VALUE(diabetic_prolif_factor, "skin.diabetic.prolif_factor");
  BDM_ASSIGN_CONFIG_VALUE(diabetic_migration_factor, "skin.diabetic.migration_factor");
  BDM_ASSIGN_CONFIG_VALUE(diabetic_fibroblast_activation_factor, "skin.diabetic.fibroblast_activation_factor");
  BDM_ASSIGN_CONFIG_VALUE(diabetic_collagen_factor, "skin.diabetic.collagen_factor");
  BDM_ASSIGN_CONFIG_VALUE(diabetic_fibroblast_lifespan_factor, "skin.diabetic.fibroblast_lifespan_factor");
  BDM_ASSIGN_CONFIG_VALUE(diabetic_baseline_inflammation, "skin.diabetic.baseline_inflammation");
  BDM_ASSIGN_CONFIG_VALUE(diabetic_neutrophil_factor, "skin.diabetic.neutrophil_factor");
  BDM_ASSIGN_CONFIG_VALUE(diabetic_neutrophil_lifespan_factor, "skin.diabetic.neutrophil_lifespan_factor");
  BDM_ASSIGN_CONFIG_VALUE(diabetic_neutrophil_apoptosis_factor, "skin.diabetic.neutrophil_apoptosis_factor");
  BDM_ASSIGN_CONFIG_VALUE(diabetic_neutrophil_waves_factor, "skin.diabetic.neutrophil_waves_factor");
  BDM_ASSIGN_CONFIG_VALUE(diabetic_neutrophil_window_factor, "skin.diabetic.neutrophil_window_factor");
  BDM_ASSIGN_CONFIG_VALUE(diabetic_macrophage_taper_factor, "skin.diabetic.macrophage_taper_factor");
  BDM_ASSIGN_CONFIG_VALUE(diabetic_macrophage_emigration_factor, "skin.diabetic.macrophage_emigration_factor");
  BDM_ASSIGN_CONFIG_VALUE(diabetic_biofilm_clearance_factor, "skin.diabetic.biofilm_clearance_factor");
  BDM_ASSIGN_CONFIG_VALUE(diabetic_vegf_factor, "skin.diabetic.vegf_factor");
  BDM_ASSIGN_CONFIG_VALUE(diabetic_ph_recovery_factor, "skin.diabetic.ph_recovery_factor");
  BDM_ASSIGN_CONFIG_VALUE(diabetic_inflammation_sensitivity, "skin.diabetic.inflammation_sensitivity");
  BDM_ASSIGN_CONFIG_VALUE(diabetic_migration_infl_K, "skin.diabetic.migration_infl_K");
  BDM_ASSIGN_CONFIG_VALUE(diabetic_prolif_infl_K, "skin.diabetic.prolif_infl_K");
  BDM_ASSIGN_CONFIG_VALUE(diabetic_mmp_factor, "skin.diabetic.mmp_factor");
  BDM_ASSIGN_CONFIG_VALUE(diabetic_timp_factor, "skin.diabetic.timp_factor");

  // Fibroblast / TGF-beta / Collagen (modules/fibroblast.toml -> [skin.fibroblast]; TOML: hours)
  BDM_ASSIGN_CONFIG_VALUE(fibroblast_enabled, "skin.fibroblast.enabled");
  BDM_ASSIGN_HOURS(fibroblast_spawn_delay, "skin.fibroblast.spawn_delay_h");
  BDM_ASSIGN_CONFIG_VALUE(fibroblast_spawn_waves, "skin.fibroblast.spawn_waves");
  BDM_ASSIGN_HOURS(fibroblast_spawn_window, "skin.fibroblast.spawn_window_h");
  BDM_ASSIGN_CONFIG_VALUE(fibroblast_diameter, "skin.fibroblast.diameter");
  BDM_ASSIGN_CONFIG_VALUE(fibroblast_density_factor, "skin.fibroblast.density_factor");
  BDM_ASSIGN_HOURS(fibroblast_activation_delay, "skin.fibroblast.activation_delay_h");
  BDM_ASSIGN_CONFIG_VALUE(fibroblast_activation_threshold, "skin.fibroblast.activation_threshold");
  BDM_ASSIGN_HOURS(fibroblast_myofibroblast_delay, "skin.fibroblast.myofibroblast_delay_h");
  BDM_ASSIGN_CONFIG_VALUE(fibroblast_myofibroblast_threshold, "skin.fibroblast.myofibroblast_threshold");
  BDM_ASSIGN_CONFIG_VALUE(fibroblast_apoptosis_threshold, "skin.fibroblast.apoptosis_threshold");
  BDM_ASSIGN_HOURS(fibroblast_apoptosis_onset, "skin.fibroblast.apoptosis_onset_h");
  BDM_ASSIGN_CONFIG_VALUE(fibroblast_apoptosis_rate, "skin.fibroblast.apoptosis_rate");
  BDM_ASSIGN_HOURS(fibroblast_min_lifespan, "skin.fibroblast.min_lifespan_h");
  BDM_ASSIGN_HOURS(fibroblast_lifespan, "skin.fibroblast.lifespan_h");
  BDM_ASSIGN_CONFIG_VALUE(fibroblast_migration_speed, "skin.fibroblast.migration_speed");
  BDM_ASSIGN_CONFIG_VALUE(dermal_fibroblast_depth, "skin.fibroblast.dermal_depth");
  BDM_ASSIGN_CONFIG_VALUE(dermal_fibroblast_margin, "skin.fibroblast.dermal_margin");
  BDM_ASSIGN_CONFIG_VALUE(tgfb_diffusion, "skin.fibroblast.tgfb_diffusion");
  BDM_ASSIGN_CONFIG_VALUE(tgfb_decay, "skin.fibroblast.tgfb_decay");
  BDM_ASSIGN_CONFIG_VALUE(tgfb_wound_seed, "skin.fibroblast.tgfb_wound_seed");
  BDM_ASSIGN_CONFIG_VALUE(fibroblast_tgfb_rate, "skin.fibroblast.tgfb_rate");
  BDM_ASSIGN_CONFIG_VALUE(m2_tgfb_rate, "skin.fibroblast.m2_tgfb_rate");
  BDM_ASSIGN_CONFIG_VALUE(collagen_deposition_rate, "skin.fibroblast.collagen_deposition_rate");
  BDM_ASSIGN_CONFIG_VALUE(collagen_decay, "skin.fibroblast.collagen_decay");

  // Tumor module (modules/tumor.toml -> [skin.tumor])
  BDM_ASSIGN_CONFIG_VALUE(tumor_enabled, "skin.tumor.enabled");
  BDM_ASSIGN_HOURS(tumor_seed_time, "skin.tumor.seed_time_h");
  BDM_ASSIGN_CONFIG_VALUE(tumor_seed_x, "skin.tumor.seed_x");
  BDM_ASSIGN_CONFIG_VALUE(tumor_seed_y, "skin.tumor.seed_y");
  BDM_ASSIGN_CONFIG_VALUE(tumor_seed_z, "skin.tumor.seed_z");
  BDM_ASSIGN_CONFIG_VALUE(tumor_seed_count, "skin.tumor.seed_count");
  BDM_ASSIGN_CONFIG_VALUE(tumor_diameter, "skin.tumor.diameter");
  BDM_ASSIGN_CONFIG_VALUE(tumor_cycle_factor, "skin.tumor.cycle_factor");
  BDM_ASSIGN_CONFIG_VALUE(tumor_max_neighbors, "skin.tumor.max_neighbors");
  BDM_ASSIGN_CONFIG_VALUE(tumor_ci_steepness, "skin.tumor.ci_steepness");
  BDM_ASSIGN_CONFIG_VALUE(tumor_growth_rate, "skin.tumor.growth_rate");
  BDM_ASSIGN_CONFIG_VALUE(tumor_max_cells, "skin.tumor.max_cells");
  BDM_ASSIGN_HOURS(tumor_handoff_delay, "skin.tumor.handoff_delay_h");
  BDM_ASSIGN_CONFIG_VALUE(tumor_stratum_value, "skin.tumor.stratum_value");
  BDM_ASSIGN_CONFIG_VALUE(tumor_apoptosis_rate, "skin.tumor.apoptosis_rate");

  // MMP dynamics (modules/mmp.toml -> [skin.mmp])
  BDM_ASSIGN_CONFIG_VALUE(mmp_enabled, "skin.mmp.enabled");
  BDM_ASSIGN_CONFIG_VALUE(mmp_diffusion, "skin.mmp.diffusion");
  BDM_ASSIGN_CONFIG_VALUE(mmp_decay, "skin.mmp.decay");
  BDM_ASSIGN_CONFIG_VALUE(mmp_m1_rate, "skin.mmp.m1_rate");
  BDM_ASSIGN_CONFIG_VALUE(mmp_neutrophil_rate, "skin.mmp.neutrophil_rate");
  BDM_ASSIGN_CONFIG_VALUE(mmp_fibroblast_rate, "skin.mmp.fibroblast_rate");
  BDM_ASSIGN_CONFIG_VALUE(mmp_keratinocyte_rate, "skin.mmp.keratinocyte_rate");
  BDM_ASSIGN_CONFIG_VALUE(mmp_collagen_degradation, "skin.mmp.collagen_degradation");
  BDM_ASSIGN_CONFIG_VALUE(mmp_fibronectin_degradation, "skin.mmp.fibronectin_degradation");

  // Fibronectin (modules/fibronectin.toml -> [skin.fibronectin])
  BDM_ASSIGN_CONFIG_VALUE(fibronectin_enabled, "skin.fibronectin.enabled");
  BDM_ASSIGN_CONFIG_VALUE(fibronectin_decay, "skin.fibronectin.decay");
  BDM_ASSIGN_CONFIG_VALUE(fibronectin_deposition_rate, "skin.fibronectin.deposition_rate");
  BDM_ASSIGN_CONFIG_VALUE(fibronectin_serum_rate, "skin.fibronectin.serum_rate");
  BDM_ASSIGN_CONFIG_VALUE(fibronectin_migration_boost, "skin.fibronectin.migration_boost");
  BDM_ASSIGN_CONFIG_VALUE(fibronectin_wound_seed, "skin.fibronectin.wound_seed");
  BDM_ASSIGN_CONFIG_VALUE(fibronectin_carrying_capacity, "skin.fibronectin.carrying_capacity");
  BDM_ASSIGN_CONFIG_VALUE(collagen_fn_transition, "skin.fibronectin.collagen_fn_transition");

  // Elastin (modules/elastin/config.toml -> [skin.elastin])
  BDM_ASSIGN_CONFIG_VALUE(elastin_enabled, "skin.elastin.enabled");
  BDM_ASSIGN_CONFIG_VALUE(elastin_diffusion, "skin.elastin.diffusion");
  BDM_ASSIGN_CONFIG_VALUE(elastin_decay, "skin.elastin.decay");
  BDM_ASSIGN_CONFIG_VALUE(elastin_basal_density, "skin.elastin.basal_density");
  BDM_ASSIGN_CONFIG_VALUE(elastin_papillary_density, "skin.elastin.papillary_density");
  BDM_ASSIGN_CONFIG_VALUE(elastin_production_rate, "skin.elastin.production_rate");
  BDM_ASSIGN_CONFIG_VALUE(elastin_mmp_degradation, "skin.elastin.mmp_degradation");

  // Hyaluronan (modules/hyaluronan/config.toml -> [skin.hyaluronan])
  BDM_ASSIGN_CONFIG_VALUE(hyaluronan_enabled, "skin.hyaluronan.enabled");
  BDM_ASSIGN_CONFIG_VALUE(hyaluronan_diffusion, "skin.hyaluronan.diffusion");
  BDM_ASSIGN_CONFIG_VALUE(hyaluronan_decay, "skin.hyaluronan.decay");
  BDM_ASSIGN_CONFIG_VALUE(hyaluronan_basal_density, "skin.hyaluronan.basal_density");
  BDM_ASSIGN_CONFIG_VALUE(hyaluronan_reticular_density, "skin.hyaluronan.reticular_density");
  BDM_ASSIGN_CONFIG_VALUE(hyaluronan_production_rate, "skin.hyaluronan.production_rate");
  BDM_ASSIGN_CONFIG_VALUE(hyaluronan_water_retention_factor, "skin.hyaluronan.water_retention_factor");
  BDM_ASSIGN_CONFIG_VALUE(hyaluronan_migration_scaffold_factor, "skin.hyaluronan.migration_scaffold_factor");

  // Dermis module (modules/dermis/config.toml -> [skin.dermis])
  BDM_ASSIGN_CONFIG_VALUE(dermis_enabled, "skin.dermis.enabled");
  BDM_ASSIGN_CONFIG_VALUE(dermis_diffusion, "skin.dermis.diffusion");
  BDM_ASSIGN_CONFIG_VALUE(dermis_decay, "skin.dermis.decay");
  BDM_ASSIGN_CONFIG_VALUE(dermis_papillary_density, "skin.dermis.papillary_density");
  BDM_ASSIGN_CONFIG_VALUE(dermis_reticular_density, "skin.dermis.reticular_density");
  BDM_ASSIGN_CONFIG_VALUE(dermis_hypodermis_density, "skin.dermis.hypodermis_density");
  BDM_ASSIGN_CONFIG_VALUE(dermis_collagen_threshold, "skin.dermis.collagen_threshold");
  BDM_ASSIGN_CONFIG_VALUE(dermis_collagen_recovery_rate, "skin.dermis.collagen_recovery_rate");
  BDM_ASSIGN_CONFIG_VALUE(dermis_mmp_degradation, "skin.dermis.mmp_degradation");
  BDM_ASSIGN_CONFIG_VALUE(dermis_papillary_rate_factor, "skin.dermis.papillary_rate_factor");
  BDM_ASSIGN_CONFIG_VALUE(dermis_reticular_rate_factor, "skin.dermis.reticular_rate_factor");
  BDM_ASSIGN_CONFIG_VALUE(dermis_hypodermis_rate_factor, "skin.dermis.hypodermis_rate_factor");

  // pH module (modules/ph/config.toml -> [skin.ph])
  BDM_ASSIGN_CONFIG_VALUE(ph_recovery_rate, "skin.ph.recovery_rate");
  BDM_ASSIGN_CONFIG_VALUE(ph_migration_suppression, "skin.ph.migration_suppression");
  BDM_ASSIGN_CONFIG_VALUE(ph_mmp_boost, "skin.ph.mmp_boost");
  BDM_ASSIGN_CONFIG_VALUE(ph_biofilm_boost, "skin.ph.biofilm_boost");
  BDM_ASSIGN_CONFIG_VALUE(ph_bohr_factor, "skin.ph.bohr_factor");

  // Hemostasis / fibrin clot (modules/hemostasis/config.toml -> [skin.hemostasis])
  BDM_ASSIGN_CONFIG_VALUE(hemostasis_enabled, "skin.hemostasis.enabled");
  BDM_ASSIGN_CONFIG_VALUE(hemostasis_decay, "skin.hemostasis.decay");
  BDM_ASSIGN_CONFIG_VALUE(hemostasis_wound_seed, "skin.hemostasis.wound_seed");
  BDM_ASSIGN_CONFIG_VALUE(hemostasis_tgfb_coupling, "skin.hemostasis.tgfb_coupling");
  BDM_ASSIGN_CONFIG_VALUE(hemostasis_fibronectin_coupling, "skin.hemostasis.fibronectin_coupling");
  BDM_ASSIGN_CONFIG_VALUE(hemostasis_mmp_degradation, "skin.hemostasis.mmp_degradation");
  BDM_ASSIGN_CONFIG_VALUE(hemostasis_migration_boost, "skin.hemostasis.migration_boost");

  // Scar formation (modules/scar.toml -> [skin.scar])
  BDM_ASSIGN_CONFIG_VALUE(scar_enabled, "skin.scar.enabled");
  BDM_ASSIGN_CONFIG_VALUE(scar_proportional_enabled, "skin.scar.proportional_enabled");
  BDM_ASSIGN_CONFIG_VALUE(scar_accumulation_rate, "skin.scar.accumulation_rate");
  BDM_ASSIGN_CONFIG_VALUE(scar_collagen_threshold, "skin.scar.collagen_threshold");

  // Metrics export (TOML: hours)
  BDM_ASSIGN_HOURS(metrics_interval, "skin.metrics_interval_h");
  BDM_ASSIGN_CONFIG_VALUE(metrics_autoopen, "skin.metrics_autoopen");

  // Cytokine taper rates
  BDM_ASSIGN_CONFIG_VALUE(neutrophil_cytokine_taper, "skin.immune.neutrophil_cytokine_taper");
  BDM_ASSIGN_CONFIG_VALUE(m1_cytokine_taper, "skin.immune.m1_cytokine_taper");
  BDM_ASSIGN_CONFIG_VALUE(m2_tgfb_taper, "skin.immune.m2_tgfb_taper");
  BDM_ASSIGN_CONFIG_VALUE(fibroblast_tgfb_taper, "skin.fibroblast.tgfb_taper_rate");

  // PDE sub-cycling
  BDM_ASSIGN_CONFIG_VALUE(subcycle_slow, "skin.subcycle_slow");
  BDM_ASSIGN_CONFIG_VALUE(subcycle_medium, "skin.subcycle_medium");

  // Hot-reload
  BDM_ASSIGN_CONFIG_VALUE(hot_reload, "skin.hot_reload");

  // Derived composite fields
  BDM_ASSIGN_CONFIG_VALUE(ecm_weight_collagen, "skin.derived.ecm_weight_collagen");
  BDM_ASSIGN_CONFIG_VALUE(ecm_weight_fibronectin, "skin.derived.ecm_weight_fibronectin");
  BDM_ASSIGN_CONFIG_VALUE(ecm_weight_elastin, "skin.derived.ecm_weight_elastin");
  BDM_ASSIGN_CONFIG_VALUE(ecm_weight_fibrin, "skin.derived.ecm_weight_fibrin");
  BDM_ASSIGN_CONFIG_VALUE(ecm_migration_boost, "skin.derived.ecm_migration_boost");

  // Multi-resolution structural fields
  BDM_ASSIGN_CONFIG_VALUE(grid_resolution_structural, "skin.grid_resolution_structural");

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

  // Wound geometry
  check(wound_radius > 0, "wound_radius must be positive");
  check(wound_center_x >= 0, "wound_center_x must be non-negative");
  check(wound_center_y >= 0, "wound_center_y must be non-negative");

  // Simulation grid
  check(grid_resolution >= 1, "grid_resolution must be >= 1");
  check(num_steps >= 1, "num_steps must be >= 1");

  // Decay rates must be non-negative
  check(inflammation_decay >= 0, "inflammation_decay must be non-negative");
  check(oxygen_decay >= 0, "oxygen_decay must be non-negative");
  check(water_decay >= 0, "water_decay must be non-negative");
  check(tgfb_decay >= 0, "tgfb_decay must be non-negative");
  check(collagen_decay >= 0, "collagen_decay must be non-negative");
  check(mmp_decay >= 0, "mmp_decay must be non-negative");

  // Diffusion coefficients must be non-negative
  check(inflammation_diffusion >= 0, "inflammation_diffusion must be non-negative");
  check(oxygen_diffusion >= 0, "oxygen_diffusion must be non-negative");
  check(mmp_diffusion >= 0, "mmp_diffusion must be non-negative");
  check(tgfb_diffusion >= 0, "tgfb_diffusion must be non-negative");
  check(vegf_diffusion >= 0, "vegf_diffusion must be non-negative");

  // Diabetic factors must be positive (multipliers)
  if (diabetic_mode) {
    check(diabetic_m1_duration_factor > 0, "diabetic_m1_duration_factor must be positive");
    check(diabetic_resolution_factor > 0, "diabetic_resolution_factor must be positive");
    check(diabetic_prolif_factor > 0, "diabetic_prolif_factor must be positive");
    check(diabetic_migration_factor > 0, "diabetic_migration_factor must be positive");
    check(diabetic_migration_infl_K > 0, "diabetic_migration_infl_K must be positive (Hill half-max)");
    check(diabetic_prolif_infl_K > 0, "diabetic_prolif_infl_K must be positive (Hill half-max)");
  }

  // Probabilities in [0,1]
  check(neutrophil_apoptosis_rate >= 0 && neutrophil_apoptosis_rate <= 1,
        "neutrophil_apoptosis_rate must be in [0,1]");
  check(macrophage_apoptosis_rate >= 0 && macrophage_apoptosis_rate <= 1,
        "macrophage_apoptosis_rate must be in [0,1]");
  check(macrophage_emigration_rate >= 0 && macrophage_emigration_rate <= 1,
        "macrophage_emigration_rate must be in [0,1]");

  // Cell sizes
  check(division_diameter > 0, "division_diameter must be positive");
  check(immune_cell_diameter > 0, "immune_cell_diameter must be positive");

  // Perfusion
  check(perfusion_basal >= 0 && perfusion_basal <= 2.0,
        "perfusion_basal should be in [0,2] (normalized)");
  check(perfusion_angio_rate >= 0, "perfusion_angio_rate must be non-negative");
}
