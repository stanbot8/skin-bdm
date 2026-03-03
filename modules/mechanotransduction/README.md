> [Home](../../README.md) / [Modules](../README.md) / Mechanotransduction

# Mechanotransduction

Tissue stiffness derived from ECM composition (collagen and elastin), YAP/TAZ mechanosensing for myofibroblast differentiation, and mechanical tension effects on scar formation.

## Biology

Tissue stiffness is a mechanical property determined primarily by extracellular matrix composition. Collagen is the dominant structural protein that increases tissue rigidity, while elastin provides compliance and reduces stiffness. During wound healing, progressive collagen deposition by myofibroblasts increases local tissue stiffness from the soft provisional matrix (~0.1 kPa) toward scar tissue values (~100 kPa).

Cells sense substrate stiffness through mechanotransduction pathways. On stiff substrates, YAP/TAZ transcriptional co-activators translocate to the nucleus where they drive expression of pro-fibrotic genes including alpha-SMA, promoting myofibroblast differentiation (Dupont et al. 2011). Myofibroblasts generate traction forces that contract the wound (Hinz 2007), and mechanical tension on healing wounds amplifies scar formation by reducing myofibroblast apoptosis (Aarabi et al. 2007). This creates a positive feedback loop: collagen deposition increases stiffness, which promotes further myofibroblast differentiation, which deposits more collagen.

## Model

**Stiffness field:** A derived structural property (non-diffusing, no decay) computed each step from local collagen and elastin concentrations:

```
stiffness = clamp(collagen * 0.8 - elastin * 0.2, 0, 1)
```

Collagen dominates stiffness (weight 0.8), while elastin provides compliance (weight -0.2). The field is recalculated in each step rather than evolved as a PDE.

**YAP/TAZ threshold:** When local stiffness exceeds `stiffness_yap_threshold` (default 0.3), nuclear YAP/TAZ signaling promotes myofibroblast differentiation and sustains existing myofibroblasts.

**Wound contraction:** Myofibroblast traction forces contract the wound at a rate governed by `stiffness_contraction_rate`.

**Scar amplification:** Mechanical tension amplifies scar formation by the `stiffness_scar_factor`, representing reduced apoptosis under load.

## Parameters

From modules/mechanotransduction/config.toml:

| Parameter | Default | Units | Description | Source |
|-----------|---------|-------|-------------|--------|
| `enabled` | true | bool | Master switch | Convention |
| `stiffness_yap_threshold` | 0.3 | a.u. | Stiffness above this promotes myofibroblast via YAP/TAZ | Dupont et al. 2011 ([DOI](https://doi.org/10.1038/nature10137)) |
| `stiffness_contraction_rate` | 0.0003 | per step | Wound contraction from myofibroblast traction force | Hinz 2007 ([DOI](https://doi.org/10.1038/sj.jid.5700613)) |
| `stiffness_scar_factor` | 0.5 | multiplier | Scar amplification from mechanical tension | Aarabi et al. 2007 ([DOI](https://doi.org/10.1096/fj.07-8218com)) |

## Coupling

### Reads
| Field | Source module | How used |
|-------|-------------|----------|
| Collagen | fibroblast | Primary stiffness component (weight 0.8) |
| Elastin | elastin | Compliance component (weight -0.2) |

### Writes
| Field | Consumer modules | What is written |
|-------|-----------------|-----------------|
| Stiffness | fibroblast, scar | Derived tissue stiffness from ECM composition |

## Source files

| File | Purpose |
|------|---------|
| `stiffness_pde.h` | Stiffness field definition (derived structural property, D=0, decay=0) |
| `source_hook.h` | Computes stiffness from collagen and elastin each step |
| `params.h` | Parameter struct (MechanotransductionParams) |
| `config.toml` | Module configuration |
