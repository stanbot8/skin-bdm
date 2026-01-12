> [Home](../../README.md) / [Modules](../README.md) / Inflammation

# Inflammation

Cytokine signaling field representing pro-inflammatory mediators (TNF-alpha, IL-1, IL-6) and their effect on keratinocyte behavior.

## Biology

When tissue is damaged, dying cells release damage-associated molecular patterns (DAMPs) and keratinocytes secrete constitutive IL-1 alpha from disrupted membranes. These early signals recruit neutrophils within hours, which in turn release TNF-alpha and IL-1 beta to amplify the inflammatory response. Monocytes arrive shortly after and differentiate into M1 macrophages that sustain high cytokine output while phagocytosing debris and pathogens. The resulting inflammatory milieu peaks around day 2 to 3 post-wounding.

Inflammation plays a dual role in wound healing. At moderate levels, cytokines recruit additional immune cells and prime fibroblasts for later matrix deposition. However, elevated inflammation actively suppresses keratinocyte proliferation and migration through Hill-function gating, meaning that re-epithelialization cannot proceed efficiently until inflammatory signals resolve. Resolution depends on macrophage polarization from M1 (pro-inflammatory) to M2 (anti-inflammatory), which typically occurs by day 3 to 5. In pathological conditions such as diabetic wounds, this transition stalls and chronic inflammation delays or prevents closure entirely.

The model introduces an ImmunePressure field that mirrors the Inflammation PDE in diffusion and decay but is written only by immune cell agents, excluding wound-derived DAMPs. Keratinocytes read ImmunePressure rather than raw Inflammation for their proliferation and migration gates. This separation prevents a biologically implausible feedback loop in which tissue damage signals suppress the very cells responsible for repairing the wound. Without this split, early DAMP release would immediately inhibit marginal keratinocytes and stall closure from the start.

An optional split pro/anti-inflammatory mode decomposes the single Inflammation field into separate ProInflammatory (TNF-alpha, IL-1) and AntiInflammatory (IL-10, TGF-beta) channels. The net inflammatory signal is computed as `max(0, pro - weight * anti)`, allowing explicit modeling of anti-inflammatory cytokine dynamics. This mode is used by the diabetic module, where impaired M2 polarization results in insufficient anti-inflammatory output and persistent net inflammation.

## Model

- Inflammation PDE with diffusion coefficient 0.30 and decay rate 0.010 per step
- Wound source term: DAMPs from open wound surface (`wound_source_rate * barrier_deficit`) with temporal taper via `wound_source_taper`
- ImmunePressure PDE: identical diffusion and decay, written only by immune agents (neutrophils and macrophages)
- Optional split: ProInflammatory and AntiInflammatory channels, with net inflammation calculated as `max(0, pro - weight * anti)`
- Hill-function gating of keratinocyte behavior: migration suppression with K = `migration_threshold`, proliferation suppression with K = `prolif_threshold`

## Parameters

From `modules/inflammation/config.toml`:

| Parameter | Default | Units | Description | Source |
|-----------|---------|-------|-------------|--------|
| `diffusion` | 0.30 | - | Lateral cytokine spread | Oyler-Yaniv et al. 2017 ([DOI](https://doi.org/10.1016/j.immuni.2017.03.011)) |
| `decay` | 0.010 | per step | Background clearance | Bemelmans et al. 1996 ([DOI](https://doi.org/10.1615/CritRevImmunol.v16.i1.10)) |
| `migration_threshold` | 0.1 | a.u. | Hill K for 50% migration suppression | Calibrated |
| `prolif_threshold` | 0.15 | a.u. | Hill K for 50% proliferation suppression | Calibrated |
| `wound_source_rate` | 0.008 | per step | DAMPs from open wound surface | Werner & Grose 2003 ([DOI](https://doi.org/10.1152/physrev.00031.2002)) |
| `wound_source_taper` | 0.001 | per step | DAMP clearance rate | Calibrated |
| `split_enabled` | false | bool | Use separate pro/anti fields | Mosser & Zhang 2008 ([DOI](https://doi.org/10.1111/j.1600-065X.2008.00706.x)) |
| `anti_diffusion` | 0.05 | - | Anti-inflammatory diffusion | Calibrated |
| `anti_decay` | 0.002 | per step | Anti-inflammatory decay | Huhn et al. 1996 ([DOI](https://doi.org/10.1182/blood.V87.2.699.bloodjournal872699)) |
| `anti_weight` | 1.0 | - | Weight for net = pro - weight*anti | Calibrated |

## Coupling

### Reads

| Field | Source module | How used |
|-------|-------------|----------|
| Stratum | tissue | Wound area detection for DAMP source term |

### Writes

| Field | Consumer modules | What is written |
|-------|-----------------|-----------------|
| Inflammation | immune (macrophage recruitment, M1 to M2 gating), scar, biofilm | Cumulative cytokine concentration |
| ImmunePressure | tissue (migration/proliferation gating) | Immune-only cytokine concentration (excludes DAMPs) |
| ProInflammatory | (when split enabled) | Pro-inflammatory cytokines |
| AntiInflammatory | (when split enabled) | Anti-inflammatory cytokines |

## Validation

| Dataset | Observable | Sources | Notes |
|---------|-----------|---------|-------|
| `inflammation_timecourse` | Peak and decay of cytokine levels | Eming 2007, Kondo 1996, Hubner 1996, Koh 2011, Canedo-Dorantes 2019 | Peak day 2 to 3, resolution by day 5 to 7 |
| `inflammation_parameters` | Diffusion/decay calibration | Oyler-Yaniv 2017, Bemelmans 1996, Mosser 2008 | PDE coefficient derivation |

## Literature data

Reference curves for validation (full citations in [SOURCES.yaml](SOURCES.yaml)):

| Dataset | File | Normalization |
|---------|------|---------------|
| Inflammation timecourse | [inflammation_timecourse.csv](data/inflammation_timecourse.csv) | Peak = 1.0 |

<details>
<summary>Raw digitized data (3 papers)</summary>

| File | Source |
|------|--------|
| [eming2007_inflammation.csv](data/raw/inflammation/eming2007_inflammation.csv) | Eming et al. 2007 |
| [kondo1996_inflammation.csv](data/raw/inflammation/kondo1996_inflammation.csv) | Kondo et al. 1996 |
| [hubner1996_inflammation.csv](data/raw/inflammation/hubner1996_inflammation.csv) | Hubner et al. 1996 |
</details>

## Metrics

| Column | Units | Description |
|--------|-------|-------------|
| `mean_infl_wound` | a.u. | Mean inflammation in wound cylinder |
| `mean_anti_infl_wound` | a.u. | Mean anti-inflammatory in wound (when split enabled) |

## Source files

| File | Purpose |
|------|---------|
| `inflammation_pde.h` | Inflammation PDE, ImmunePressure PDE, split pro/anti PDEs |
| `config.toml` | Module configuration |
