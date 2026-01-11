> [Home](../../README.md) / [Modules](../README.md) / Wound

# Wound

Circular punch biopsy event that initiates the wound healing cascade.

## Biology

A punch biopsy creates a full-thickness circular wound that destroys the epithelium and damages the underlying dermis and vasculature. This technique is the standard clinical and experimental model for studying wound healing because it produces a reproducible, well-defined lesion with clean margins. Removing a cylinder of tissue eliminates the protective stratum corneum, disrupts calcium gradients that maintain barrier homeostasis, and severs dermal capillaries to create a hypoxic wound bed.

The wound event serves as the trigger for the entire healing cascade. Loss of the epithelial barrier exposes underlying tissue to the environment, releasing damage-associated molecular patterns (DAMPs) that recruit neutrophils and macrophages. Keratinocytes at the wound margin activate, lose their cell-cell adhesion, and begin migrating inward to re-cover the wound surface. These marginal cells also proliferate to supply additional daughter cells, which are biased toward the wound center. As the migrating sheet advances and cells stack vertically through crowding-driven extrusion, the wound gradually closes with a stratified epithelium that matures into a functional barrier.

## Model

- WoundEvent zeros Stratum, Calcium, O2 fields in the wound cylinder at the configured trigger step
- Spawns basal keratinocytes at the wound margin, spacing them by one cell diameter around the perimeter
- Inward division bias (`wound_inward_bias`) preferentially directs daughter cells toward the wound center, reproducing the centripetal migration sheet observed in vivo
- Crowding-driven vertical extrusion builds stratified layers as the wound bed fills
- WoundResolution provides a safety timeout that removes residual agents 200 steps before simulation end or at greater than 90% Stratum coverage, whichever comes first

## Parameters

From `modules/wound/config.toml`:

| Parameter | Default | Units | Description | Source |
|-----------|---------|-------|-------------|--------|
| `enabled` | true | bool | Master switch | Convention |
| `center_x` | 15.0 | um | X-center of punch biopsy | Convention |
| `center_y` | 15.0 | um | Y-center of punch biopsy | Convention |
| `radius` | 5.0 | um | Radius of circular wound | Convention |
| `trigger_h` | 0 | hours | Hours at which wound fires | Convention |
| `inward_bias` | 0.3 | 0 to 1 | Division bias toward wound center | Calibrated |
| `vascular_damage` | 0.5 | fraction | Legacy vascular damage (superseded by perfusion module) | Calibrated |

## Coupling

### Reads

| Field | Source module | How used |
|-------|-------------|----------|
| Stratum | tissue | Coverage check for wound closure percentage |

### Writes

| Field | Consumer modules | What is written |
|-------|-----------------|-----------------|
| Stratum | tissue, scar | Zeros in wound cylinder |
| Calcium | tissue | Zeros in wound cylinder |
| O2 | tissue | Zeros in wound cylinder |
| Vascular | perfusion | Damage proportional to `vascular_damage` |

## Validation

| Dataset | Observable | Sources | Notes |
|---------|-----------|---------|-------|
| `closure_kinetics_punch_biopsy` | Wound closure percentage over time | Eaglstein 1978, Cukjati 2000, Pastar 2014 | Sigmoid closure curve, approximately 22 days to 90% |

## Literature data

Reference curves for validation (full citations in [SOURCES.yaml](SOURCES.yaml)):

| Dataset | File | Normalization |
|---------|------|---------------|
| Closure kinetics (punch biopsy) | [closure_kinetics_punch_biopsy.csv](data/closure_kinetics_punch_biopsy.csv) | Absolute 0 to 100% |

<details>
<summary>Raw digitized data (3 papers)</summary>

| File | Source |
|------|--------|
| [eaglstein1978_closure.csv](data/raw/closure/eaglstein1978_closure.csv) | Eaglstein et al. 1978 |
| [cukjati2000_closure.csv](data/raw/closure/cukjati2000_closure.csv) | Cukjati et al. 2000 |
| [pastar2014_closure.csv](data/raw/closure/pastar2014_closure.csv) | Pastar et al. 2014 |
</details>

## Metrics

| Column | Units | Description |
|--------|-------|-------------|
| `wound_closure_pct` | % | Fraction of wound voxels with Stratum > 0.5 |

## Source files

| File | Purpose |
|------|---------|
| `wound_event.h` | Wound creation, field zeroing, agent spawning, resolution |
| `config.toml` | Module configuration |
