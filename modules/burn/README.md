> [Home](../../README.md) / [Modules](../README.md) / Burn

# Burn

Jackson's three-zone burn wound model (1953) with zone-specific tissue fate, eschar formation, transepidermal water loss, and contracture risk. Disabled by default; activated for burn injury studies.

## Biology

Thermal injury produces three concentric zones of tissue damage (Jackson 1953):

1. **Zone of coagulation** (center): Irreversible protein denaturation and cell death. Tissue is non-viable and forms eschar (dry necrotic crust) that blocks inflammatory cell access.

2. **Zone of stasis** (intermediate): Tissue with impaired perfusion but initially viable cells. Without adequate resuscitation, progressive ischemia converts this zone to coagulation over 24 to 48 hours. This zone is the primary target of early burn treatment.

3. **Zone of hyperemia** (peripheral): Viable tissue with enhanced inflammatory response and increased blood flow. Generally recovers fully within 7 to 10 days.

Burn depth determines clinical severity. Superficial burns (epidermal only) heal rapidly from intact appendages. Partial-thickness burns damage the dermis but retain some regenerative structures. Full-thickness burns destroy all epidermal appendages and require surgical intervention (Herndon 2012). Barrier destruction causes massive transepidermal water loss (TEWL), and the compromised barrier greatly increases infection susceptibility (Singer and Clark 1999).

## Model

**Zone classification:** Each wound voxel is classified by normalized depth (z-coordinate relative to wound depth) and the `depth_fraction` parameter:

```
Zone 0 (Coagulation):  norm_z < depth_fraction * 0.5     Irreversible necrosis
Zone 1 (Stasis):       norm_z < depth_fraction            Deteriorating viability
Zone 2 (Hyperemia):    norm_z < depth_fraction + 0.2      Enhanced inflammation
Zone 3 (Unburned):     norm_z >= depth_fraction + 0.2     Normal tissue
```

**Zone of coagulation:** Perfusion is driven toward zero by progressively killing vascular supply (necrosis fraction 0.95). Eschar formation gradually suppresses local inflammation.

**Zone of stasis:** A runtime viability state (initialized at `stasis_initial_viability`) deteriorates at `stasis_deterioration_rate` per hour. Local perfusion scales with remaining viability. Without treatment, viability reaches zero over 24 to 48 hours, converting the zone to functional coagulation.

**Zone of hyperemia:** Inflammation is boosted at `hyperemia_inflammation_boost` rate, modeling the enhanced inflammatory response in viable peripheral tissue.

**Barrier effects:** Destroyed epithelial barrier causes fluid loss at `fluid_loss_rate * tewl_multiplier * (1 - barrier)`, and infection susceptibility is amplified by `infection_susceptibility` factor.

**Contracture:** In full-thickness burns (depth_fraction >= 0.7), scar accumulates in deep zones at `contracture_rate`, representing the excessive contraction tendency of deep burn scars.

## Parameters

From modules/burn/config.toml:

| Parameter | Default | Units | Description | Source |
|-----------|---------|-------|-------------|--------|
| `enabled` | false | bool | Master switch (opt-in for burn studies) | Convention |
| `coagulation_necrosis` | 0.95 | fraction | Irreversible cell death fraction in coagulation zone | Jackson 1953 ([DOI](https://doi.org/10.1002/bjs.18004016413)) |
| `stasis_initial_viability` | 0.5 | fraction | Zone of stasis initial viability (salvageable) | Jackson 1953 ([DOI](https://doi.org/10.1002/bjs.18004016413)) |
| `stasis_deterioration_rate` | 0.02 | per hour | Ischemic deterioration without treatment | Jackson 1953 ([DOI](https://doi.org/10.1002/bjs.18004016413)) |
| `hyperemia_inflammation_boost` | 0.3 | a.u. | Inflammatory enhancement in hyperemia zone | Jackson 1953 ([DOI](https://doi.org/10.1002/bjs.18004016413)) |
| `depth_fraction` | 0.5 | fraction | Burn depth (0 to 0.3 superficial, 0.3 to 0.7 partial, 0.7+ full) | Singer and Clark 1999 ([DOI](https://doi.org/10.1056/NEJM199909023411006)) |
| `eschar_rate` | 0.01 | per step | Eschar formation from necrotic tissue | Herndon 2012 (ISBN 978-1-4377-2786-9) |
| `tewl_multiplier` | 5.0 | multiplier | Transepidermal water loss through destroyed barrier | Herndon 2012 (ISBN 978-1-4377-2786-9) |
| `fluid_loss_rate` | 0.05 | per step | Fluid loss rate through burn surface | Herndon 2012 (ISBN 978-1-4377-2786-9) |
| `infection_susceptibility` | 3.0 | multiplier | Barrier breach infection risk | Guo and DiPietro 2010 ([DOI](https://doi.org/10.1177/0022034509359125)) |
| `contracture_rate` | 0.002 | per step | Excessive scar contraction tendency in deep burns | Herndon 2012 (ISBN 978-1-4377-2786-9) |

## Coupling

### Reads
| Field | Source module | How used |
|-------|-------------|----------|
| Vascular (perfusion) | angiogenesis | Zone-specific perfusion modification |
| Barrier | epithelial | Gates TEWL and fluid loss |

### Writes
| Field | Consumer modules | What is written |
|-------|-----------------|-----------------|
| Vascular (perfusion) | angiogenesis, oxygen | Reduced in coagulation/stasis zones |
| Inflammation | immune, senescence | Boosted in hyperemia, suppressed by eschar in coagulation |
| Water | hydration | Fluid loss through destroyed barrier (TEWL) |
| Scar | scar | Contracture risk accumulation in deep burns |

## Source files

| File | Purpose |
|------|---------|
| `source_hook.h` | Three-zone burn model with perfusion, inflammation, TEWL, and contracture logic |
| `params.h` | Parameter struct (BurnParams) |
| `config.toml` | Module configuration |
