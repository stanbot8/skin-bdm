> [Home](../../README.md) / [Modules](../README.md) / Bioelectric

# Bioelectric

Transepithelial potential (TEP), endogenous wound electric fields, and galvanotaxis. Intact epithelium maintains a lateral voltage gradient at wound edges that guides cell migration toward the wound center.

## Biology

Intact skin epithelium maintains a transepithelial potential (TEP) of approximately 40 to 70 mV through Na+/K+ ATPase ion pumps in epidermal cells (Barker et al. 1982). When a wound disrupts the epithelial barrier, the TEP collapses locally, creating a lateral voltage gradient of 100 to 200 mV/mm at the wound edge. This endogenous electric field (EF) persists until re-epithelialization restores the barrier.

The lateral voltage gradient drives galvanotaxis, the directed migration of cells along the electric field. Epithelial cells, fibroblasts, and immune cells orient their migration cathodally (toward the wound center) through PI3K-gamma and PTEN signaling pathways (Zhao et al. 2006). Genetic disruption of this signaling delays wound closure, demonstrating that bioelectric cues are essential guidance signals that complement chemical gradients.

In diabetic wounds, impaired Na+/K+ ATPase activity and ion channel dysfunction reduce the endogenous wound current, contributing to delayed re-epithelialization (Reid and Zhao 2014).

## Model

**Voltage PDE:** A diffusing field representing the transepithelial potential. Fast diffusion (0.1) models rapid ionic current spread through extracellular fluid (McCaig et al. 2005). Decay (0.01) represents charge leakage through tissue resistance.

```
Source: intact epithelium   [voltage_epithelial_source * stratum_gate]
Sink:   tissue resistance   [voltage_decay]
Spread: ionic current       [voltage_diffusion = 0.1]
```

**Wound edge gradient:** Intact epithelium (stratum > 0) produces voltage; wound center (stratum = 0) does not. The diffusion solver creates a natural gradient that peaks at the wound margin and decays inward, matching measured wound EF profiles.

**Galvanotaxis:** Cell migration is biased along the voltage gradient with strength `galvanotaxis_strength` (default 0.4). This provides a directional cue that supplements chemotactic signals (TGF-beta, PDGF).

**Diabetic impairment:** In diabetic mode, the voltage source is multiplied by `diabetic_voltage_factor` (default 0.6), representing reduced Na+/K+ ATPase activity.

## Parameters

From modules/bioelectric/config.toml:

| Parameter | Default | Units | Description | Source |
|-----------|---------|-------|-------------|--------|
| `enabled` | true | bool | Master switch | Convention |
| `voltage_diffusion` | 0.1 | - | Fast ionic current spread through extracellular fluid | McCaig et al. 2005 ([DOI](https://doi.org/10.1152/physrev.00020.2004)) |
| `voltage_decay` | 0.01 | per step | Charge leakage through tissue resistance | Calibrated |
| `voltage_epithelial_source` | 0.05 | per step | Na+/K+ ATPase maintains TEP in intact epithelium | Barker et al. 1982 ([DOI](https://doi.org/10.1152/ajpregu.1982.242.3.R358)) |
| `galvanotaxis_strength` | 0.4 | - | Voltage gradient bias on cell migration | Zhao et al. 2006 ([DOI](https://doi.org/10.1038/nature04925)) |
| `diabetic_voltage_factor` | 0.6 | multiplier | Impaired endogenous wound currents in diabetes | Reid and Zhao 2014 ([DOI](https://doi.org/10.1089/wound.2013.0442)) |

## Coupling

### Reads
| Field | Source module | How used |
|-------|-------------|----------|
| Stratum | epithelial | Gates voltage source (intact epithelium = source, wound = no source) |

### Writes
| Field | Consumer modules | What is written |
|-------|-----------------|-----------------|
| Voltage | fibroblast, epithelial, immune | TEP field; gradient drives galvanotaxis |

## Source files

| File | Purpose |
|------|---------|
| `voltage_pde.h` | Voltage diffusion field (fast diffusion, tissue resistance decay) |
| `source_hook.h` | Epithelial voltage source gated by stratum, with diabetic impairment |
| `params.h` | Parameter struct (BioelectricParams) |
| `config.toml` | Module configuration |
