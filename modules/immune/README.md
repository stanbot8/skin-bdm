> [Home](../../README.md) / [Modules](../README.md) / Immune

# Immune

Innate immune response to wounding, modeling neutrophil infiltration and macrophage recruitment with M1/M2 polarization.

## Biology
The innate immune response is the first cellular reaction to wounding. Neutrophils arrive within hours as first responders, producing pro-inflammatory cytokines and clearing debris. Macrophages follow, initially in a pro-inflammatory M1 state that sustains inflammation. As the wound stabilizes, macrophages transition to an anti-inflammatory M2 state that resolves inflammation and produces growth factors (TGF-beta, VEGF) to drive tissue repair. This M1 to M2 transition is a critical checkpoint: delayed or failed transition leads to chronic inflammation and impaired healing.

Cytokine production uses age-dependent exponential tapers rather than flat rates. Neutrophils produce at `cytokine_rate * exp(-taper * age)`, so young neutrophils at the wound site are the most active producers while NF-kB signaling declines as they age toward apoptosis. M1 macrophages produce at `cytokine_rate * exp(-taper * state_age)`, so output declines with time in the M1 state. The inflammation peak-and-decay timecourse emerges from the interplay of these individual cell tapers with population dynamics.

Efferocytosis (phagocytosis of apoptotic neutrophils) is a key M1 to M2 transition trigger. M1 macrophages search for dying neutrophils via ForEachNeighbor; on detection, the neutrophil is engulfed and the macrophage transitions to M2.

## Model
Agent-based immune cells with two types (Neutrophil, Macrophage) and macrophage states (M1, M2).

**Neutrophil behavior:**
- Spawn in waves (default 3 over 24h) starting at `neutrophil_spawn_delay_h` post-wound
- Age-dependent cytokine production with exponential taper
- After `neutrophil_min_survival_h`, stochastic apoptosis at `neutrophil_apoptosis_rate`
- Hard ceiling at `neutrophil_lifespan_h`

**Macrophage behavior:**
- Continuous recruitment driven by inflammation: spawn probability = `spawn_rate * max(0, infl - threshold) * max(0, 1 - n_mac / capacity)`
- Carrying capacity models limited ICAM-1 adhesion sites at postcapillary venules
- M1 state: pro-inflammatory cytokine production with age taper
- M1 to M2 transition triggers: (1) inflammation below `m1_transition_threshold` with min age, (2) efferocytosis of dying neutrophil, (3) `macrophage_m1_duration_h` ceiling fallback
- M2 state: actively consumes inflammation from the field at `resolution_rate`
- M2 macrophages also produce TGF-beta and VEGF

**Chemotaxis:** When enabled, immune cells follow the inflammation gradient via BioDynaMo's DiffusionGrid::GetGradient(). Falls back to geometric wound-center migration when gradient is flat.

### Mechanistic toggles

Two optional mechanistic replacements can be enabled for testing. Both default to `false` (parametric behavior unchanged).

**Gradient-driven recruitment** (`mech_immune_recruitment`): Replaces the threshold + rate + taper recruitment model with chemokine gradient magnitude driving monocyte extravasation via Michaelis-Menten saturation: `prob = scale * G/(K+G) * saturation * taper`. The endothelial adhesion taper (ICAM-1 decline) is shared by both modes.

**Efferocytosis-count M1 to M2** (`mech_m1_m2_transition`): Replaces cytokine-threshold M1 to M2 with engulfment-count driven transition. M1 macrophages accumulate an engulfment counter via TryEfferocytosis; once the quota is met (default 1), each step has a `mech_m2_transition_rate` probability of transitioning to M2 (PS receptor signaling via MerTK/Tim-4 reprogramming PPARgamma/LXR). The timer ceiling (`macrophage_m1_duration_h`) remains as fallback for macrophages that never encounter neutrophils, representing alternative M2 signals (IL-4, IL-10, glucocorticoids).

## Parameters
From modules/immune/config.toml:
| Parameter | Default | Units | Description | Source |
|-----------|---------|-------|-------------|--------|
| `cell_diameter` | 3.0 | um | Neutrophil/macrophage diameter | Convention |
| `neutrophil_spawn_delay_h` | 2 | hours | Hours post-wound | Kim et al. 2008 ([DOI](https://doi.org/10.1038/sj.jid.5701223)) |
| `neutrophil_spawn_waves` | 3 | count | Recruitment waves | Calibrated |
| `neutrophil_spawn_window_h` | 24 | hours | Time span for all waves | Calibrated |
| `neutrophil_lifespan_h` | 48 | hours | Hard ceiling | Wilgus et al. 2013 |
| `neutrophil_min_survival_h` | 12 | hours | Min before apoptosis possible | Calibrated |
| `neutrophil_apoptosis_rate` | 0.0029 | per step | Death probability (half-life ~24h) | Wilgus et al. 2013 |
| `macrophage_spawn_delay_h` | 12 | hours | Hours post-wound | Rodero & Khosrotehrani 2010 |
| `macrophage_spawn_threshold` | 0.002 | a.u. | Inflammation for recruitment | Calibrated |
| `macrophage_spawn_rate` | 0.8 | - | Recruitment probability scaling | Calibrated |
| `macrophage_spawn_taper` | 0.02 | per hour | Recruitment probability decay | Calibrated |
| `macrophage_m1_duration_h` | 48 | hours | M1 duration ceiling | Krzyszczyk et al. 2018 ([DOI](https://doi.org/10.3389/fphys.2018.00419)) |
| `m1_transition_threshold` | 0.003 | a.u. | Inflammation below this triggers M1 to M2 | Calibrated |
| `m1_transition_min_age_h` | 36 | hours | Min M1 time before transition | Krzyszczyk 2018 |
| `macrophage_lifespan_h` | 672 | hours | Hard ceiling (28 days) | Convention |
| `macrophage_min_survival_h` | 24 | hours | Min before apoptosis possible | Calibrated |
| `macrophage_apoptosis_rate` | 0.0008 | per step | Death probability (half-life ~87h) | Lucas et al. 2010 ([DOI](https://doi.org/10.4049/jimmunol.0903356)) |
| `cytokine_rate` | 0.004 | per step | Inflammation per cell (before taper) | Eming et al. 2007 |
| `resolution_rate` | 0.020 | per step | Inflammation consumed by M2 | Koh & DiPietro 2011 ([DOI](https://doi.org/10.1017/S1462399411001943)) |
| `migration_speed` | 1.5 | - | Tractor force magnitude | Lammermann et al. 2008 ([DOI](https://doi.org/10.1038/nature06887)) |
| `efferocytosis_enabled` | true | bool | Proximity-based M1 to M2 | Convention |
| `efferocytosis_radius` | 5.0 | um | Search radius | Ravichandran 2010 ([DOI](https://doi.org/10.1084/jem.20101157)) |
| `efferocytosis_age_fraction` | 0.8 | fraction | Neutrophil age for PS exposure | Savill et al. 1989 ([DOI](https://doi.org/10.1172/JCI113970)) |
| `chemotaxis_enabled` | true | bool | Follow inflammation gradient | Convention |
| `chemotaxis_speed_scale` | 1.0 | - | Gradient speed scaling | Calibrated |
| `neutrophil_cytokine_taper` | 0.004 | - | Exponential decay constant | Calibrated |
| `m1_cytokine_taper` | 0.003 | - | Exponential decay constant | Calibrated |
| `mech_immune_recruitment` | false | bool | Gradient-driven recruitment (test mode) | Convention |
| `mech_recruit_gradient_scale` | 0.12 | - | Gradient to probability scaling | Calibrated |
| `mech_recruit_saturation_k` | 0.005 | a.u. | Half-saturation for gradient | Calibrated |
| `mech_m1_m2_transition` | false | bool | Efferocytosis-count M1 to M2 (test mode) | Convention |
| `mech_efferocytosis_quota` | 1 | count | Engulfments to trigger M2 program | Ravichandran 2010 ([DOI](https://doi.org/10.1084/jem.20101157)) |
| `mech_m2_transition_rate` | 0.25 | per step | Transition probability once quota met | Calibrated |

## Coupling
### Reads
| Field | Source module | How used |
|-------|-------------|----------|
| Inflammation | inflammation | Macrophage recruitment threshold, M1 to M2 gating |
| Biofilm | biofilm | Blocks M1 to M2 when above threshold |

### Writes
| Field | Consumer modules | What is written |
|-------|-----------------|-----------------|
| Inflammation | inflammation, scar, biofilm | Cytokine production (age-tapered) |
| ImmunePressure | tissue (migration/proliferation gating) | Cytokine production (mirrors inflammation, excludes DAMPs) |
| TGF-beta | fibroblast | M2 macrophage production |
| VEGF | angiogenesis | M2 macrophage production |
| MMP | mmp | M1 macrophage MMP-9 production |

## Validation
| Dataset | Observable | Sources | Notes |
|---------|-----------|--------|-------|
| `immune_cell_kinetics` | Neutrophil and macrophage counts over time | Kim 2008, Wilgus 2013, Krzyszczyk 2018, Rodero 2010, Lucas 2010 | Peak day 2, decline by day 5 |
| `inflammation_timecourse` | Emergent cytokine curve shape | Eming 2007, Koh 2011 | Peak-and-decay from cell aging |

## Literature data

Reference curves for validation (full citations in [SOURCES.yaml](SOURCES.yaml)):

| Dataset | File | Normalization |
|---------|------|---------------|
| Immune cell kinetics | [immune_cell_kinetics.csv](data/immune_cell_kinetics.csv) | Peak = 1.0 per type |

<details>
<summary>Raw digitized data (5 papers)</summary>

| File | Source |
|------|--------|
| [kim2008_neutrophils.csv](data/raw/neutrophils/kim2008_neutrophils.csv) | Kim et al. 2008 |
| [wilgus2013_neutrophils.csv](data/raw/neutrophils/wilgus2013_neutrophils.csv) | Wilgus et al. 2013 |
| [krzyszczyk2018_macrophages.csv](data/raw/macrophages/krzyszczyk2018_macrophages.csv) | Krzyszczyk et al. 2018 |
| [rodero2010_macrophages.csv](data/raw/macrophages/rodero2010_macrophages.csv) | Rodero et al. 2010 |
| [lucas2010_macrophages.csv](data/raw/macrophages/lucas2010_macrophages.csv) | Lucas et al. 2010 |
</details>

## Metrics
| Column | Units | Description |
|--------|-------|-------------|
| `n_neutrophils` | count | Active neutrophil agents |
| `n_macrophages` | count | Active macrophage agents |

## Source files
| File | Purpose |
|------|---------|
| `immune_cell.h` | ImmuneCell agent with type (Neutrophil/Macrophage) and state (M1/M2) |
| `neutrophil_behavior.h` | Neutrophil lifecycle: spawn, cytokine production, apoptosis |
| `macrophage_behavior.h` | Macrophage lifecycle: M1/M2 polarization, efferocytosis, TGF-beta/VEGF |
| `immune_response.h` | Immune response orchestrator: spawn timing, wave scheduling |
| `immune_helpers.h` | Cytokine production and resolution helper functions |
| `config.toml` | Module configuration |
