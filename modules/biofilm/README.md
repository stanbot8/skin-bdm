> [Home](../../README.md) / [Modules](../README.md) / Biofilm

# Biofilm

Bacterial biofilm colonization of wound beds with immune clearance and inflammatory feedback.

## Biology

Chronic wounds are frequently colonized by bacterial biofilms, with Staphylococcus aureus and Pseudomonas aeruginosa being the most common pathogens. Biofilm formation follows a delayed time course, with bacterial colonization typically beginning 48 hours after wounding when the exposed wound bed provides a suitable environment. Bacteria grow logistically within the wound matrix, and the biofilm produces pathogen-associated molecular patterns (PAMPs) that sustain pro-inflammatory signaling.

At high biofilm densities, the M1 to M2 macrophage transition is blocked regardless of inflammation level, trapping the immune response in a destructive pro-inflammatory state. In diabetic wounds, impaired neutrophil phagocytosis reduces biofilm clearance, creating a chronic non-healing feedback loop between persistent infection and unresolved inflammation.

## Model

**Biofilm PDE:** logistic growth with no diffusion (biofilm is sessile). Growth rate follows `growth_rate * biofilm * (1 - biofilm / carrying_capacity)`.

**Seeding:** at `seed_delay_h` post-wound (~48h), small inoculum (`seed_amount`) seeds wound-bed voxels where stratum < 1.0 (open wound surface).

**Immune clearance:** neutrophils clear at `neutrophil_clearance` per cell per step; macrophages clear at `macrophage_clearance` per cell per step. In diabetic mode, clearance rates are multiplied by `biofilm_clearance_factor` (0.5x).

**Inflammatory feedback:** biofilm produces PAMPs at `inflammation_rate * local_biofilm` per step, sustaining inflammation.

**M1 block:** when local biofilm exceeds `m1_block_threshold`, M1 to M2 transition is blocked, preventing immune resolution.

**Feedback loop:**
```
Wound --> Biofilm seeding --> Logistic growth --> PAMPs --> Sustained M1
                                                             |
                            Blocked M1 to M2 <-- biofilm > threshold
                                                             |
                       Neutrophils + Macrophages --> Clearance
                            (diabetic: impaired)
```

## Parameters

From modules/biofilm/config.toml:

| Parameter | Default | Units | Description | Source |
|-----------|---------|-------|-------------|--------|
| `enabled` | false | bool | Master switch | Convention |
| `growth_rate` | 0.005 | per step | Logistic growth (~2h doubling) | Gibson et al. 2018 ([DOI](https://doi.org/10.1098/rspb.2018.0789)) |
| `carrying_capacity` | 1.0 | normalized | Max density | Robson 1997 ([DOI](https://doi.org/10.1016/S0039-6109(05)70572-7)) |
| `seed_delay_h` | 48 | hours | Hours post-wound before inoculation | Davis et al. 2008 ([DOI](https://doi.org/10.1111/j.1524-475X.2007.00303.x)) |
| `seed_amount` | 0.01 | normalized | Initial inoculum density | Calibrated |
| `neutrophil_clearance` | 0.003 | per step | Clearance per neutrophil | Jesaitis et al. 2003 ([DOI](https://doi.org/10.4049/jimmunol.171.8.4329)) |
| `macrophage_clearance` | 0.001 | per step | Clearance per macrophage | Calibrated |
| `m1_block_threshold` | 0.3 | normalized | Biofilm level blocking M1 to M2 | Bjarnsholt et al. 2008 ([DOI](https://doi.org/10.1111/j.1524-475X.2007.00283.x)) |
| `inflammation_rate` | 0.001 | per step | PAMP-driven inflammation per biofilm | Calibrated |

## Coupling

### Reads
| Field | Source module | How used |
|-------|-------------|----------|
| Stratum | tissue | Open wound detection for seeding (stratum < 1.0) |

### Writes
| Field | Consumer modules | What is written |
|-------|-----------------|-----------------|
| Biofilm | immune (M1 block when above threshold) | Biofilm density per wound voxel |
| Inflammation | inflammation | PAMP-driven pro-inflammatory signal |

## Validation

| Dataset | Observable | Sources | Notes |
|---------|-----------|--------|-------|
| `biofilm_parameters` | Growth kinetics, immune clearance rates | James 2008, Bjarnsholt 2008, Robson 1997, Gibson 2018, Davis 2008, Jesaitis 2003 | Parameter derivation |

## Metrics

| Column | Units | Description |
|--------|-------|-------------|
| `mean_biofilm_wound` | a.u. | Mean biofilm in wound |

## Source files

| File | Purpose |
|------|---------|
| `biofilm_pde.h` | Biofilm PDE with logistic growth |
| `biofilm_op.h` | Seeding, immune clearance, inflammation feedback, M1 block |
| `config.toml` | Module configuration |
