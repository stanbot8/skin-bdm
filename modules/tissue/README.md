> [Home](../../README.md) / [Modules](../README.md) / Tissue

# Tissue

Core epidermal biology module representing the stratified squamous epithelium of human skin.

## Biology

Human epidermis is a self-renewing stratified epithelium organized into four principal strata. At the base, the stratum basale houses proliferative keratinocytes anchored to the basement membrane by integrin adhesion complexes. As cells divide, daughters are displaced upward into the stratum spinosum, where they enlarge and begin expressing differentiation markers. Continued ascent carries cells into the stratum granulosum, where they extrude lipid-filled lamellar bodies and assemble the cornified envelope. Finally, cells in the stratum corneum are terminally differentiated, anucleate corneocytes that form the primary permeability barrier before desquamating from the skin surface.

Extracellular calcium is the master regulator of this vertical program. The epidermis maintains a characteristic calcium gradient, low in the basale and rising sharply through the granulosum, that is sustained by tight junctions and ion pumps rather than by free diffusion. When extracellular calcium rises above successive thresholds, keratinocytes activate differentiation gene programs (involucrin, loricrin, filaggrin) and withdraw from the cell cycle. Below the epidermis, dermal fibroblasts secrete keratinocyte growth factor (KGF/FGF7), a paracrine mitogen that diffuses upward and binds FGFR2b on basal cells, accelerating their transit from G1 into S phase and sustaining the proliferative pool.

The basal layer itself is organized as a stem and transit-amplifying (TA) cell hierarchy. Stem cells divide asymmetrically, producing one daughter that retains stem identity and one that becomes a TA cell committed to a finite number of divisions before exiting the cycle. This architecture ensures long-term tissue renewal while buffering against clonal exhaustion. Oxygen and water, supplied by the dermal vasculature, impose metabolic constraints on proliferation and migration: hypoxia suppresses cell cycle entry, while dehydration from barrier loss impairs both crawling and mitosis.

The UWYN (Use When You Need) hybrid model leverages the observation that healthy epidermis exists in a steady state that can be fully described by continuum PDE fields (calcium, KGF, oxygen, water, stratum identity) without tracking individual cells. Discrete keratinocyte agents are spawned only when an event such as wounding disrupts the continuum equilibrium. After agents have restored local tissue integrity, cornified cells dissolve back into the continuum following a configurable handoff delay. This approach dramatically reduces computational cost in regions of the tissue that are quiescent while preserving single-cell resolution where it matters.

## Model

The tissue module maintains five continuum fields and a suite of agent behaviors that together reproduce epidermal homeostasis and wound re-epithelialization. The **VolumeManager** partitions the simulation domain into horizontal layer profiles using z-height boundaries: basale below 6 um, spinosum from 6 to 15 um, granulosum from 15 to 25 um, and corneum above 25 um. The **calcium** field is initialized as a sigmoid profile rising from 0.05 mM at the basement membrane to 1.5 mM at the tissue surface, with midpoint and steepness parameters controlling the shape of the curve. **KGF** is initialized as an exponential decay from 2.0 nM at z=0, reflecting its dermal paracrine origin. **Oxygen** and **water** are both initialized as exponential decays from dermal baseline values, with their dermal source voxels pinned proportionally to the local Vascular perfusion field so that vascular damage or recovery directly modulates nutrient supply. The **stratum** field is a static scalar encoding layer identity at each voxel, used for visualization and gating of calcium recovery.

Each keratinocyte agent carries an explicit four-phase cell cycle (G1, S, G2, M) with literature-derived phase durations. Transition probabilities accumulate over time and are modulated multiplicatively by local field values: KGF boosts G1 to S via Michaelis-Menten kinetics, oxygen gates proliferation linearly below a hypoxia threshold, water gates proliferation below a dehydration threshold, and immune pressure (from the inflammation module) suppresses division through a Hill-function dose response. During S phase the cell grows in volume until it reaches a minimum division diameter. At M phase completion, the cell attempts division subject to contact inhibition: if the neighbor count exceeds `max_neighbors`, the cell enters quiescence (G0) instead. Successful divisions orient with a lateral scatter randomness and a wound-inward bias when applicable, and crowding in the basal layer triggers vertical extrusion of daughter cells to initiate stratification.

The stem and TA hierarchy governs clonal dynamics. Stem cells (identified by `divisions_left = -1`) divide indefinitely and are anchored to the basement membrane by simulated integrin adhesion. On division, a stem cell produces a TA daughter (with probability `p_asymmetric`) that inherits a finite division budget of `max_ta_divisions`. TA daughters decrement their count at each division and enter G0 when exhausted. The `stem_fraction` parameter controls the initial ratio of stem to TA cells when agents are seeded at the wound margin, with low-calcium voxels biased toward stem identity to reflect the basal niche environment.

Migration is implemented as a tractor force directed toward the wound center, composing naturally with BioDynaMo's mechanical repulsion. Speed scales with radial distance (leader cells at the wound edge move fastest, with a 30% floor to prevent convergence stall) and is gated by water availability, fibronectin concentration (integrin-mediated crawling enhancement), and immune pressure (Hill-function suppression). The **differentiation** behavior reads local calcium and z-height to assign stratum identity, with calcium able to override volume boundaries near layer transitions. Cornified cells that exceed the `handoff_delay` residence time dissolve back into the continuum, completing the per-cell UWYN lifecycle. The **shedding** behavior removes cornified cells after a configurable desquamation delay and clears exhausted TA cells that linger in the basale past the apoptosis timeout.

## Parameters

Core tissue parameters live in `bdm.core.toml` under `[skin]`. The `modules/tissue/config.toml` file is minimal, serving as a placeholder for the merge script.

### Cell cycle

| Parameter | Default | Units | Description | Source |
|-----------|---------|-------|-------------|--------|
| `g1_duration` | 7.0 | hours | G1 phase duration | Grabe & Bhatt-Neuber 2005 |
| `s_duration` | 6.0 | hours | S phase duration | Grabe & Bhatt-Neuber 2005 |
| `g2_duration` | 3.0 | hours | G2 phase duration | Grabe & Bhatt-Neuber 2005 |
| `m_duration` | 1.0 | hours | M phase duration | Grabe & Bhatt-Neuber 2005 |

### Division mechanics

| Parameter | Default | Units | Description | Source |
|-----------|---------|-------|-------------|--------|
| `growth_rate` | 5 | - | Volume growth rate during S phase | Calibrated |
| `max_neighbors` | 14 | count | Contact inhibition threshold | Eisenhoffer et al. 2012 |
| `lateral_scatter` | 0.3 | - | Lateral randomness in division axis | Calibrated |
| `max_ta_divisions` | 4 | count | Max TA divisions before G0 | Jones et al. 1993 |
| `division_diameter` | 5.0 | um | Minimum diameter to divide | Calibrated |
| `p_asymmetric` | 0.7 | probability | Asymmetric stem division probability | Clayton et al. 2007 |
| `stem_fraction` | 0.50 | fraction | Fraction of seeded basal cells that are stem | Safferling et al. 2013 |

### Calcium gradient

| Parameter | Default | Units | Description | Source |
|-----------|---------|-------|-------------|--------|
| `calcium_diffusion` | 0 | - | Diffusion coefficient (0 = static prescribed) | Convention |
| `calcium_decay` | 0 | - | Decay rate (0 = persistent gradient) | Convention |
| `calcium_basal` | 0.05 | mM | Concentration at basement membrane | Menon et al. 1985 |
| `calcium_peak` | 1.5 | mM | Concentration at tissue surface | Menon et al. 1985 |
| `calcium_midpoint_z` | 15.0 | um | Z-height where sigmoid reaches 50% | Calibrated |
| `calcium_steepness` | 2.0 | - | Sigmoid steepness | Calibrated |

### Differentiation thresholds

| Parameter | Default | Units | Description | Source |
|-----------|---------|-------|-------------|--------|
| `ca_spinous_threshold` | 0.1 | mM | Ca2+ for spinous transition | Bikle et al. 2012 |
| `ca_granular_threshold` | 0.5 | mM | Ca2+ for granular transition | Bikle et al. 2012 |
| `ca_cornified_threshold` | 1.0 | mM | Ca2+ for cornified transition | Bikle et al. 2012 |
| `spinous_threshold` | 6.0 | um | Z-height reinforcement for basal identity | Calibrated |

### Volume boundaries

| Parameter | Default | Units | Description | Source |
|-----------|---------|-------|-------------|--------|
| `volume_z_spinous` | 6.0 | um | Basal/spinous boundary | Calibrated |
| `volume_z_granular` | 15.0 | um | Spinous/granular boundary | Calibrated |
| `volume_z_cornified` | 25.0 | um | Granular/cornified boundary | Calibrated |

### KGF

| Parameter | Default | Units | Description | Source |
|-----------|---------|-------|-------------|--------|
| `kgf_diffusion` | 0 | - | Diffusion (0 = static) | Convention |
| `kgf_decay` | 0 | - | Decay (0 = persistent) | Convention |
| `kgf_basal_conc` | 2.0 | nM | Concentration at basement membrane | Calibrated |
| `kgf_half_maximal` | 0.5 | nM | Michaelis-Menten Km | Calibrated |
| `kgf_max_boost` | 1.0 | fold | Max fold-increase in G1 to S rate | Calibrated |
| `kgf_decay_length` | 5.0 | um | Z-scale for exponential decay | Calibrated |

### Oxygen

| Parameter | Default | Units | Description | Source |
|-----------|---------|-------|-------------|--------|
| `oxygen_diffusion` | 0.1 | - | Diffusion coefficient | Calibrated |
| `oxygen_decay` | 0.01 | - | Consumption rate | Calibrated |
| `oxygen_basal_conc` | 1.0 | normalized | Arterial pO2 at z=0 | Convention |
| `oxygen_decay_length` | 8.0 | um | Z-scale for initial profile | Calibrated |
| `oxygen_prolif_threshold` | 0.3 | normalized | O2 below this suppresses proliferation | Calibrated |
| `oxygen_recovery_enabled` | true | bool | Vasculature regenerates with healing | Convention |

### Water

| Parameter | Default | Units | Description | Source |
|-----------|---------|-------|-------------|--------|
| `water_diffusion` | 1e-4 | - | Lateral diffusion | Calibrated |
| `water_decay` | 0.002 | - | Background TEWL | Calibrated |
| `water_basal_conc` | 1.0 | normalized | Hydration at dermis | Convention |
| `water_decay_length` | 12.0 | um | Z-scale for initial profile | Calibrated |
| `water_recovery_rate` | 0.02 | per step | Serum hydration rate | Sakai et al. 2005 |
| `water_surface_loss_rate` | 0.03 | per step | Evaporation at exposed surface | Calibrated |
| `water_migration_threshold` | 0.3 | normalized | Min moisture for full migration | Calibrated |
| `water_prolif_threshold` | 0.4 | normalized | Min moisture for full proliferation | Calibrated |

### Mechanics

| Parameter | Default | Units | Description | Source |
|-----------|---------|-------|-------------|--------|
| `repulsion_coeff` | 5.0 | - | Hertz-like steric repulsion | Van Liedekerke et al. 2015 |
| `attraction_coeff` | 0.0 | - | Desmosome/E-cadherin adhesion | Calibrated |

### Shedding

| Parameter | Default | Units | Description | Source |
|-----------|---------|-------|-------------|--------|
| `shedding_delay` | 99999 | steps | Steps before cornified cells shed | Convention |
| `apoptosis_delay` | 99999 | steps | Steps before exhausted TA cells die | Convention |

### Dynamic coupling

| Parameter | Default | Units | Description | Source |
|-----------|---------|-------|-------------|--------|
| `calcium_recovery_rate` | 0.002 | per step | Sigmoid recovery rate in healed voxels | Calibrated |
| `migration_enabled` | true | bool | Active crawling toward wound center | Convention |
| `migration_speed` | 2.0 | - | Tractor force magnitude | Calibrated |
| `handoff_delay` | 500 | steps | Steps in stable stratum before dissolving | Calibrated |

## Coupling

### Reads

| Field | Source module | How used |
|-------|-------------|----------|
| Vascular | perfusion | O2 and Water dermal pinning proportional to local perfusion |
| ImmunePressure | inflammation | Hill-function gating of migration speed and G1 to S probability |
| Fibronectin | fibronectin | Migration speed boost via integrin-FN binding |
| Hyaluronan | hyaluronan | Water retention modulation |

### Writes

| Field | Consumer modules | What is written |
|-------|-----------------|-----------------|
| Calcium | (self) | Sigmoid gradient, wound recovery |
| KGF | (self) | Exponential decay from dermis |
| O2 | (self), fibroblast | Perfusion-gated dermal source |
| Water | (self) | Perfusion-gated dermal source, TEWL evaporation |
| Stratum | wound, scar, tumor | Layer identity from differentiating agents |

## Validation

| Dataset | Observable | Sources | Notes |
|---------|-----------|---------|-------|
| `core_tissue_parameters` | Cell cycle, calcium, KGF, stem hierarchy | Grabe 2005, Menon 1985, Bikle 2012, Adra 2010, Jones 1993, Clayton 2007, Safferling 2013 | Parameter derivation |
| `normal_skin_profile` | Turnover time, stem fraction | Dover & Potten 1988, Potten 1988 | Healthy adult epidermis |

## Metrics

| Column | Units | Description |
|--------|-------|-------------|
| `n_agents` | count | Active keratinocyte agents |
| `n_stem` | count | Stem cell agents |
| `n_ta` | count | Transit-amplifying agents |
| `n_g0` | count | Quiescent (G0) agents |
| `n_basal` | count | Agents in stratum basale |
| `n_spinous` | count | Agents in stratum spinosum |
| `n_granular` | count | Agents in stratum granulosum |
| `n_cornified` | count | Agents in stratum corneum |
| `mean_o2_wound` | normalized | Mean O2 in wound cylinder |
| `mean_ca_wound` | mM | Mean Ca2+ in wound cylinder |

## Source files

| File | Purpose |
|------|---------|
| `keratinocyte.h` | Keratinocyte agent with stratum, cell cycle, stem/TA identity |
| `basal_division.h` | G1/S/G2/M cell cycle with O2, water, KGF, inflammation gating |
| `migration.h` | Tractor force toward wound center with inflammation, water, fibronectin gating |
| `differentiation.h` | Calcium-driven stratum transitions |
| `shedding.h` | Desquamation and UWYN continuum handoff |
| `stratum.h` | Stratum continuum field PDE |
| `calcium.h` | Calcium gradient field with wound recovery |
| `kgf.h` | KGF growth factor field |
| `oxygen.h` | Oxygen diffusion with vascular source |
| `water.h` | Water diffusion, evaporation, and dermal pinning |
| `config.toml` | Module configuration (placeholder; params in bdm.core.toml) |
