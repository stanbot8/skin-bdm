> [Home](../../README.md) / [Modules](../README.md) / Neuropathy

# Neuropathy

Cutaneous nerve density, neuropeptide signaling (substance P and CGRP), nerve regeneration from wound margins, and diabetic denervation. Nerve density modulates neurogenic inflammation and vasodilation, linking innervation to healing capacity.

## Biology

Skin is densely innervated by sensory nerve fibers (C-fibers and A-delta fibers) that release neuropeptides upon activation. The two principal neuropeptides in wound healing are substance P (SP) and calcitonin gene-related peptide (CGRP).

Substance P promotes wound healing through multiple pathways: it triggers mast cell degranulation (neurogenic inflammation), boosts keratinocyte proliferation via NK1 receptor signaling, stimulates fibroblast migration, and promotes angiogenesis (Suvas 2017). CGRP is a potent vasodilator that increases local blood flow and perfusion, supporting oxygen and nutrient delivery to the wound bed.

Wounding causes Wallerian degeneration of severed nerve fibers within the wound area, creating a denervated zone. Nerve regeneration occurs slowly (weeks to months) from the wound margins, guided by Schwann cells and neurotrophic factors. VEGF promotes nerve sprouting through neurovascular coupling, while TGF-beta and scar tissue inhibit neurite extension.

In diabetic neuropathy, chronic hyperglycemia causes progressive loss of sensory fibers (approximately 70% reduction) through polyol pathway activation, AGE formation, PKC activation, and oxidative stress (Boulton et al. 2005; Callaghan et al. 2012). The resulting denervation reduces neuropeptide signaling, impairing neurogenic inflammation, vasodilation, and re-epithelialization, which contributes significantly to delayed diabetic wound healing.

## Model

**Nerve density field:** A structural field (slow diffusion 0.01, no decay) representing cutaneous sensory innervation density. Initialized at basal level in dermal tissue and zeroed in the wound area at wounding (Wallerian degeneration). Regenerates slowly from wound margins.

```
Init:   basal_nerve_density in dermis (reduced by diabetic_nerve_factor in diabetic mode)
Wound:  zeroed in wound area (Wallerian degeneration)
Regen:  regeneration_rate * VEGF_boost * TGF-beta_inhibition
        -> capped at basal_nerve_density (or diabetic target)
```

**Nerve regeneration (source hook):** Applied per dermal wound voxel. Regeneration rate is modulated by:
- **VEGF neurovascular coupling:** Local VEGF boosts nerve sprouting by `regeneration_vegf_boost`
- **TGF-beta scar barrier:** Local TGF-beta inhibits neurite extension by `regeneration_tgfb_inhibit`
- **Diabetic impairment:** Regeneration rate is reduced by `diabetic_regeneration_factor` (0.4x)

**Neuropeptide effects (post hook):** Applied per epidermal wound voxel, proportional to local nerve density:
- **Substance P neurogenic inflammation:** Deposits to Inflammation field at `substance_p_inflammation * nerve_density`
- **CGRP vasodilation:** Boosts local perfusion at `cgrp_vasodilation * nerve_density`, capped at basal perfusion

## Parameters

From modules/neuropathy/config.toml:

| Parameter | Default | Units | Description | Source |
|-----------|---------|-------|-------------|--------|
| `enabled` | true | bool | Master switch | Convention |
| `diffusion` | 0.01 | - | Slow neurite extension from wound margins | Calibrated |
| `decay` | 0.0 | per step | Structural (no spontaneous loss) | Convention |
| `basal_nerve_density` | 1.0 | normalized | Healthy dermal innervation density | Convention |
| `regeneration_rate` | 0.0002 | per step | Schwann cell-guided nerve regrowth (slow, weeks to months) | Calibrated |
| `regeneration_vegf_boost` | 0.5 | multiplier | VEGF promotes nerve sprouting (neurovascular coupling) | Calibrated |
| `regeneration_tgfb_inhibit` | 0.3 | multiplier | TGF-beta inhibits neurite extension (scar barrier) | Calibrated |
| `substance_p_inflammation` | 0.0001 | per step | SP promotes neurogenic inflammation (mast cell degranulation) | Suvas 2017 ([DOI](https://doi.org/10.4049/jimmunol.1601751)) |
| `substance_p_proliferation` | 0.15 | multiplier | SP boosts keratinocyte proliferation | Suvas 2017 ([DOI](https://doi.org/10.4049/jimmunol.1601751)) |
| `cgrp_vasodilation` | 0.01 | per step | CGRP-mediated vasodilation boosts local perfusion | Calibrated |
| `diabetic_nerve_factor` | 0.3 | fraction | Residual nerve density in diabetic tissue (70% loss) | Boulton et al. 2005 ([DOI](https://doi.org/10.2337/diacare.28.4.956)) |
| `diabetic_regeneration_factor` | 0.4 | fraction | Impaired Schwann cell function and neurotrophic support | Callaghan et al. 2012 ([DOI](https://doi.org/10.1016/s1474-4422(12)70065-0)) |

## Coupling

### Reads
| Field | Source module | How used |
|-------|-------------|----------|
| VEGF | angiogenesis | Neurovascular coupling boosts nerve sprouting |
| TGF-beta | fibroblast | Scar barrier inhibits neurite extension |
| Vascular (perfusion) | angiogenesis | Target cap for CGRP vasodilation |

### Writes
| Field | Consumer modules | What is written |
|-------|-----------------|-----------------|
| Nerve | (self) | Nerve density: initialized, zeroed at wound, regenerated |
| Inflammation | immune | Substance P neurogenic inflammation |
| Vascular (perfusion) | angiogenesis, oxygen | CGRP vasodilation boost |

## Source files

| File | Purpose |
|------|---------|
| `nerve_pde.h` | Nerve density structural field with initialization and Wallerian degeneration |
| `source_hook.h` | Nerve regeneration with VEGF boost and TGF-beta inhibition |
| `post_hook.h` | Neuropeptide effects: substance P inflammation, CGRP vasodilation |
| `params.h` | Parameter struct (NeuropathyParams) |
| `config.toml` | Module configuration |
