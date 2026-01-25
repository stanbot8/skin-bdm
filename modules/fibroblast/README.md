> [Home](../../README.md) / [Modules](../README.md) / Fibroblast

# Fibroblast

TGF-beta driven fibroblast lifecycle with myofibroblast differentiation and collagen deposition.

## Biology

Dermal fibroblasts are the primary effectors of wound repair in the dermis. After wounding, resident fibroblasts near the wound margin become activated by TGF-beta signaling from M2 macrophages. Activated fibroblasts migrate toward the wound center and differentiate into myofibroblasts, which are contractile cells that deposit collagen (the structural basis of scar tissue) and produce additional TGF-beta, creating a positive feedback loop.

The loop sustains itself via myofibroblast TGF-beta production. It breaks when the wound closes, M2 macrophages die, and TGF-beta decays below the apoptosis threshold. Myofibroblasts undergo a stochastic apoptosis program (Desmouliere et al. 1995) after approximately 7 days, producing a gradual decline from peak density rather than a sharp cutoff.

Staggered recruitment spreads fibroblasts across 4 waves over 96 hours, producing a gradual rise in myofibroblast density rather than a synchronized spike.

## Model

**Fibroblast state machine:**
```
Quiescent  [TGF-beta > activation_threshold OR timeout]  Activated
Activated  [TGF-beta > myofibroblast_threshold AND delay]  Myofibroblast
Myofibroblast  [TGF-beta < apoptosis_threshold AND age > min_lifespan]  Removed
Myofibroblast  [state_age > apoptosis_onset, stochastic]  Removed
Any state  [age > lifespan]  Removed
```

Only myofibroblasts produce TGF-beta and deposit collagen. All non-quiescent states migrate via chemotaxis on TGF-beta gradient (or geometric fallback toward wound center).

TGF-beta production uses an exponential taper: `tgfb_rate * exp(-taper_rate * state_age)`, which gradually reduces the positive feedback as the wound matures.

**TGF-beta PDE:** diffusion 0.03, decay 0.008. Sources: M2 macrophages and myofibroblasts.

**Collagen PDE:** no diffusion (structural deposit), optional MMP decay. Deposited by myofibroblasts proportional to local TGF-beta concentration.

**Feedback loop:**
```
M2 Macrophages  [m2_tgfb_rate]  TGF-beta field
                                        |
                                   Fibroblast activation
                                        |
                                   Myofibroblast differentiation
                                     /        \
        [fibroblast_tgfb_rate]                  [collagen_deposition_rate * local_tgfb]
                  |                                      |
           TGF-beta field                        Collagen field
```

**Timeline** (default parameters; Desmouliere et al. 1995, Darby et al. 2014):
- Hour 1: first wave of resident dermal fibroblasts (quiescent)
- Day 1 to 4: remaining waves arrive (4 waves over 96h)
- Day 2 to 4: fibroblasts activate (TGF-beta threshold or auto-timeout at ~72h)
- Day 5: myofibroblast differentiation begins
- Day 5+: collagen deposition proportional to local TGF-beta
- Day 12+: stochastic apoptosis begins (after 7 days as myofibroblast)
- Day 28: hard lifespan limit

## Parameters

From modules/fibroblast/config.toml:

| Parameter | Default | Units | Description | Source |
|-----------|---------|-------|-------------|--------|
| `enabled` | true | bool | Master switch | Convention |
| `spawn_delay_h` | 1 | hours | Hours post-wound for first wave | Convention |
| `spawn_waves` | 4 | count | Recruitment pulses | Calibrated |
| `spawn_window_h` | 96 | hours | Time span for all waves (~4 days) | Calibrated |
| `diameter` | 4.0 | um | Cell diameter | Convention |
| `density_factor` | 1.0 | multiplier | Spawn count multiplier | Convention |
| `activation_delay_h` | 72 | hours | Auto-activation timeout | Calibrated |
| `activation_threshold` | 0.005 | a.u. | TGF-beta for quiescent to activated | Tomasek et al. 2002 ([DOI](https://doi.org/10.1038/nrm809)) |
| `myofibroblast_delay_h` | 120 | hours | Min hours in activated state (~5 days) | Van De Water et al. 2013 ([DOI](https://doi.org/10.1089/wound.2012.0393)) |
| `myofibroblast_threshold` | 0.015 | a.u. | TGF-beta for activated to myofibroblast | Tomasek et al. 2002 |
| `apoptosis_threshold` | 0.003 | a.u. | TGF-beta below this triggers removal | Hinz 2007 ([DOI](https://doi.org/10.1038/sj.jid.5700613)) |
| `apoptosis_onset_h` | 168 | hours | Hours as myofibroblast before stochastic death (~7d) | Desmouliere et al. 1995 |
| `apoptosis_rate` | 0.0008 | per step | Removal probability once eligible | Desmouliere et al. 1995 |
| `min_lifespan_h` | 168 | hours | Min before apoptosis eligible (7 days) | Desmouliere et al. 1995 |
| `lifespan_h` | 672 | hours | Hard max lifespan (28 days) | Desmouliere et al. 1995 |
| `migration_speed` | 1.0 | - | Tractor force magnitude | McDougall et al. 2006 ([DOI](https://doi.org/10.1098/rsta.2006.1773)) |
| `dermal_depth` | -1.0 | um | Z-coordinate for seeding (papillary dermis) | Convention |
| `dermal_margin` | 1.0 | um | Extra radius beyond wound for seeding ring | Calibrated |
| `tgfb_diffusion` | 0.03 | - | TGF-beta diffusion coefficient | Murphy et al. 2012 ([DOI](https://doi.org/10.1007/s11538-011-9712-y)) |
| `tgfb_decay` | 0.008 | per step | TGF-beta decay rate | Murphy et al. 2012 |
| `tgfb_rate` | 0.005 | per step | TGF-beta per myofibroblast | Calibrated |
| `tgfb_taper_rate` | 0.0012 | - | Exponential taper for TGF-beta production | Tomasek 2002 |
| `m2_tgfb_rate` | 0.001 | per step | TGF-beta per M2 macrophage | Koh & DiPietro 2011 ([DOI](https://doi.org/10.1017/S1462399411001943)) |
| `collagen_deposition_rate` | 0.0005 | per step | Collagen per myofibroblast | Murphy et al. 2012 |
| `collagen_decay` | 0.0 | per step | MMP remodeling (0 = permanent) | Convention |

## Coupling

### Reads
| Field | Source module | How used |
|-------|-------------|----------|
| TGF-beta | (self), immune (M2) | Activation and differentiation thresholds |
| Inflammation | inflammation | Indirectly, via wound state |

### Writes
| Field | Consumer modules | What is written |
|-------|-----------------|-----------------|
| TGF-beta | (self), fibroblast activation | Myofibroblast and M2 macrophage production |
| Collagen | scar, dermis, mmp | Structural deposit from myofibroblasts |

## Validation

| Dataset | Observable | Sources | Notes |
|---------|-----------|--------|-------|
| `myofibroblast_kinetics` | Myofibroblast count over time | Desmouliere 1995, Darby 2014, Tomasek 2002, Hinz 2007 | Peak day 10 to 14, decline by day 21 |
| `collagen_deposition` | Collagen accumulation curve | Gonzalez 2016, Zhou 2013, Mathew-Steiner 2021, Caetano 2016 | Monotonic rise, plateau by day 20 to 25 |
| `scar_formation_parameters` | TGF-beta/collagen PDE calibration | Murphy 2012, Dallon 2001, Koh 2011, Ogawa 2017 | Diffusion/decay derivation |

## Literature data

Reference curves for validation (full citations in [SOURCES.yaml](SOURCES.yaml)):

| Dataset | File | Normalization |
|---------|------|---------------|
| Myofibroblast kinetics | [myofibroblast_kinetics.csv](data/myofibroblast_kinetics.csv) | Peak = 1.0 |
| Collagen deposition | [collagen_deposition.csv](data/collagen_deposition.csv) | End = 1.0 |
| TGF-beta kinetics | [tgfb_kinetics.csv](data/tgfb_kinetics.csv) | Peak = 1.0 |
| Fibroblast kinetics | [fibroblast_kinetics.csv](data/fibroblast_kinetics.csv) | Peak = 1.0 |

<details>
<summary>Raw digitized data (8 papers)</summary>

| File | Source |
|------|--------|
| [desmouliere1995_myofibroblasts.csv](data/raw/myofibroblasts/desmouliere1995_myofibroblasts.csv) | Desmouliere et al. 1995 |
| [darby2014_myofibroblasts.csv](data/raw/myofibroblasts/darby2014_myofibroblasts.csv) | Darby et al. 2014 |
| [tomasek2002_myofibroblasts.csv](data/raw/myofibroblasts/tomasek2002_myofibroblasts.csv) | Tomasek et al. 2002 |
| [hinz2007_myofibroblasts.csv](data/raw/myofibroblasts/hinz2007_myofibroblasts.csv) | Hinz et al. 2007 |
| [gonzalez2016_collagen.csv](data/raw/collagen/gonzalez2016_collagen.csv) | Gonzalez et al. 2016 |
| [zhou2013_collagen.csv](data/raw/collagen/zhou2013_collagen.csv) | Zhou et al. 2013 |
| [mathew2021_collagen.csv](data/raw/collagen/mathew2021_collagen.csv) | Mathew et al. 2021 |
| [caetano2016_collagen.csv](data/raw/collagen/caetano2016_collagen.csv) | Caetano et al. 2016 |
</details>

## Metrics

| Column | Units | Description |
|--------|-------|-------------|
| `n_fibroblasts` | count | Active fibroblast agents (all states) |
| `n_myofibroblasts` | count | Myofibroblast-state fibroblasts |
| `mean_tgfb_wound` | a.u. | Mean TGF-beta in wound |
| `mean_collagen_wound` | a.u. | Mean collagen in wound |

## Source files

| File | Purpose |
|------|---------|
| `fibroblast.h` | Fibroblast agent with activation states |
| `fibroblast_behavior.h` | State machine, TGF-beta production, collagen deposition |
| `fibroblast_recruitment.h` | Staggered wave recruitment logic |
| `tgfbeta_pde.h` | TGF-beta diffusion field |
| `collagen_pde.h` | Collagen structural deposit field |
| `config.toml` | Module configuration |
