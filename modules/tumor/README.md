> [Home](../../README.md) / [Modules](../README.md) / Tumor

# Tumor

Basal cell carcinoma (BCC) neoplastic growth model with soft contact inhibition and UWYN continuum handoff.

## Biology

Basal cell carcinoma is the most common human malignancy, arising from the basal layer of the epidermis. BCC cells have an altered cell cycle (approximately 56 hours vs 17 hours for normal keratinocytes) and exhibit soft contact inhibition mediated by the YAP/TAZ mechanotransduction pathway. Unlike hard contact inhibition where proliferation stops above a neighbor threshold, soft CI produces a gradual quadratic suppression that allows continued growth at reduced rates in crowded environments.

BCC growth is slow compared to many cancers, with volume doubling times of approximately 148 days. The proliferation index (Ki-67 positive fraction) is typically around 27%, with a cell loss factor of approximately 0.112 per day from stochastic apoptosis.

## Model

**TumorCell agent:** inherits from Cell with BCC-specific properties:
- Cell cycle phases scaled by `cycle_factor` (3.3x normal)
- Soft contact inhibition: division probability = `max(0, 1 - (n_neighbors / max_neighbors)^ci_steepness)`
- Stochastic apoptosis at `apoptosis_rate` per step
- No differentiation (BCC cells remain basal-like)

**Seeding:** at `seed_time_h`, a cluster of `seed_count` cells is placed at (seed_x, seed_y, seed_z) in a sphere distribution.

**UWYN handoff:** when a tumor cell enters G0 and remains quiescent for `handoff_delay_h`, it converts to the Tumor field (binary marker, `stratum_value` = 10.0) and the agent is removed. This allows the growing tumor periphery to use agents while the stable interior reverts to continuum.

**Tumor PDE:** binary marker field with no diffusion and no decay. Records the spatial footprint of quiescent tumor cells after handoff.

**Scale-aware validation:** the simulation operates at a representative tissue scale (30x30 um domain). Growth metrics are compared against clinical data using surface-to-volume ratio corrections.

**Scenario composability:** the tumor module can run alone (`--preset=tumor`) or combined with wound healing (`--preset=tumor_wound`), where a tumor grows at t=0 and a wound is introduced at day 8. Both cellular populations coexist, sharing the same field infrastructure.

## Parameters

From modules/tumor/config.toml:

| Parameter | Default | Units | Description | Source |
|-----------|---------|-------|-------------|--------|
| `enabled` | false | bool | Master switch | Convention |
| `seed_time_h` | 0 | hours | Time to seed tumor | Convention |
| `seed_x` | 15.0 | um | Cluster center X | Convention |
| `seed_y` | 15.0 | um | Cluster center Y | Convention |
| `seed_z` | 2.0 | um | Cluster center Z (basal layer) | Convention |
| `seed_count` | 5 | count | Initial cluster size | Calibrated |
| `diameter` | 4.0 | um | Cell diameter | Convention |
| `cycle_factor` | 3.3 | multiplier | BCC cycle ~56h | Khoo et al. 2019 ([DOI](https://doi.org/10.2340/00015555-3325)) |
| `max_neighbors` | 12 | count | Contact inhibition threshold | Convention |
| `ci_steepness` | 2.0 | exponent | Soft CI (0=hard, 2=quadratic YAP/TAZ) | Calibrated |
| `growth_rate` | 5 | - | Volume growth rate | Convention |
| `max_cells` | 5000 | count | Carrying capacity (0 = unlimited) | Convention |
| `handoff_delay_h` | 20 | hours | G0 hours before continuum conversion | Calibrated |
| `stratum_value` | 10.0 | - | Stratum field value for tumor tissue | Convention |
| `apoptosis_rate` | 0.00047 | per step | Cell loss 0.112/day | Khoo et al. 2019 |

## Coupling

### Reads
| Field | Source module | How used |
|-------|-------------|----------|
| O2 | tissue | Proliferation gating (optional) |

### Writes
| Field | Consumer modules | What is written |
|-------|-----------------|-----------------|
| Tumor | (visualization) | Binary marker from UWYN handoff |
| Stratum | tissue | Tumor tissue value (10.0) at handoff |

## Validation

| Dataset | Observable | Sources | Notes |
|---------|-----------|--------|-------|
| `tumor_growth_rate` | Volume growth over time | Fijalkowska 2023, Sykes 2020, Kricker 2014 | Slow BCC growth |
| `tumor_doubling_time` | Volume doubling time ~148 days | Khoo 2019, Tejera 2023 | Clinical measurements |
| `tumor_proliferation_index` | Ki-67 ~27% | Toth 2012, Alferraly 2019, al-Sader 1996 | Immunohistochemistry |

## Literature data

Reference curves for validation (full citations in [SOURCES.yaml](SOURCES.yaml)):

| Dataset | File | Normalization |
|---------|------|---------------|
| Tumor doubling time | [tumor_doubling_time.csv](data/tumor_doubling_time.csv) | Absolute days |
| Tumor growth rate | [tumor_growth_rate.csv](data/tumor_growth_rate.csv) | Absolute mm |
| Tumor proliferation index | [tumor_proliferation_index.csv](data/tumor_proliferation_index.csv) | Absolute % |

<details>
<summary>Raw digitized data (8 papers)</summary>

| File | Source |
|------|--------|
| [fijalkowska2023_bcc_growth_meta.csv](data/raw/tumor_growth/fijalkowska2023_bcc_growth_meta.csv) | Fijalkowska et al. 2023 |
| [sykes2020_superficial_bcc.csv](data/raw/tumor_growth/sykes2020_superficial_bcc.csv) | Sykes et al. 2020 |
| [kricker2014_community_bcc.csv](data/raw/tumor_growth/kricker2014_community_bcc.csv) | Kricker et al. 2014 |
| [al_qahtani2020_bcc_doubling.csv](data/raw/tumor_doubling/al_qahtani2020_bcc_doubling.csv) | Al-Qahtani et al. 2020 |
| [tejera2023_skin_cancer_tdt.csv](data/raw/tumor_doubling/tejera2023_skin_cancer_tdt.csv) | Tejera et al. 2023 |
| [toth2012_bcc_ki67.csv](data/raw/tumor_proliferation/toth2012_bcc_ki67.csv) | Toth et al. 2012 |
| [alferraly2019_scc_ki67.csv](data/raw/tumor_proliferation/alferraly2019_scc_ki67.csv) | Alferraly et al. 2019 |
| [alsader1996_bcc_scc_indices.csv](data/raw/tumor_proliferation/alsader1996_bcc_scc_indices.csv) | Al-Sader et al. 1996 |
</details>

## Metrics

| Column | Units | Description |
|--------|-------|-------------|
| `n_tumor_cells` | count | Active tumor cell agents |
| `n_tumor_cycling` | count | Tumor agents not in G0 (Ki-67 proxy) |
| `tumor_field_cells` | count | Tumor field voxels > 0.5 (handoff footprint) |

## Source files

| File | Purpose |
|------|---------|
| `tumor_cell.h` | TumorCell agent with BCC-specific cell cycle and soft CI |
| `tumor_behavior.h` | Division, apoptosis, and UWYN handoff logic |
| `tumor_initiation.h` | Sphere seeding at configured time and location |
| `tumor_pde.h` | Tumor binary marker field |
| `config.toml` | Module configuration |
