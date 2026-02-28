# Treatment Interventions for Diabetic Wounds

Computational treatment study modeling 8 therapeutic interventions for diabetic foot ulcers. Each treatment is a TOML overlay applied on top of the diabetic wound profile, modifying specific biological parameters based on published mechanisms of action.

## Available treatments

| Treatment | Targets | Key references |
|-----------|---------|----------------|
| `anti_inflammatory` | M1/M2 transition, TNF-alpha, baseline inflammation | Mirza & Koh 2011, Louiselle et al. 2021 |
| `hbo` | Tissue oxygenation, VEGF, angiogenesis, collagen | Catrina et al. 2004, Thangarajah et al. 2009 |
| `npwt` | Blood flow, granulation, moisture, wound contraction | Morykwas et al. 1997, Armstrong & Lavery 2005 |
| `doxycycline` | MMP inhibition, collagen/fibronectin preservation | Lobmann et al. 2002, Smith et al. 1999 |
| `growth_factor` | Fibroblast activation, PDGF-BB mitogenesis | Steed 2006, Smiell 1998 |
| `msc` | Multi-target paracrine (immune, vascular, matrix) | Cao et al. 2017, Li et al. 2024 |
| `moisture` | TEWL reduction, hydration, exudate management | Winter 1962, Junker et al. 2013 |
| `combination` | All four dysfunction axes simultaneously | Multi-target rational combination |

## Running treatments

```bash
# Single treatment
./run.sh --preset=diabetic_wound --skin=diabetic --treatment=hbo

# Full comparison study (all 8 treatments + baseline)
python3 scripts/study/treatment_study.py

# Specific treatments only
python3 scripts/study/treatment_study.py --treatments=hbo,msc,combination

# List available treatments
./run.sh --list-treatments
```

## How treatments work

Treatments are TOML overlays in `treatments/` applied after the diabetic profile and preset:

```
bdm.core.toml + modules/*/config.toml -> bdm.toml (merge)
  + profiles/diabetic.toml          (biology)
  + presets/diabetic_wound.toml     (scenario)
  + treatments/hbo.toml             (intervention)
```

Each treatment modifies the diabetic dysfunction parameters back toward healthy values, proportional to the treatment's clinical effect size from published literature.

## Treatment mechanisms

### Anti-inflammatory (anti-TNF-alpha)

The core diabetic wound pathology is failed M1-to-M2 macrophage transition. TNF-alpha neutralization directly addresses this by:
- Releasing the M1-to-M2 brake (m1_duration_factor: 3.0 -> 1.5)
- Restoring M2 anti-inflammatory resolution (resolution_factor: 0.3 -> 0.7)
- Reducing AGE/RAGE-driven baseline inflammation (0.001 -> 0.0002)
- Restoring keratinocyte proliferation and migration

### Hyperbaric oxygen (HBO)

HBO corrects the tissue hypoxia that drives much of diabetic wound dysfunction:
- Restores HIF-1alpha/VEGF cycling (vegf_factor: 0.4 -> 1.0)
- Improves microcirculation (perfusion basal: 0.7 -> 0.85)
- Accelerates angiogenesis (angio_rate: 0.002 -> 0.008)
- Partially rescues collagen synthesis via prolyl hydroxylase O2 supply

### Negative pressure wound therapy (NPWT)

Multi-modal mechanical intervention:
- 4x blood flow increase from sub-atmospheric pressure
- Microdeformation stimulates cell proliferation
- Macrodeformation assists wound edge advancement (inward_bias: 0.3 -> 0.5)
- Sealed dressing controls moisture and clears inflammatory exudate

### Doxycycline (sub-antimicrobial MMP inhibitor)

Targets the MMP/TIMP imbalance that destroys ECM in diabetic wounds:
- Halves excess MMP activity (mmp_factor: 3.0 -> 1.5)
- Restores TIMP balance (timp_factor: 0.5 -> 0.8)
- Preserves collagen and fibronectin scaffolds
- Mild anti-inflammatory effect via TACE inhibition

### Growth factor (PDGF-BB / becaplermin)

FDA-approved recombinant growth factor for DFU:
- Potent fibroblast mitogen (density_factor: 1.0 -> 1.8)
- Rapid fibroblast activation (activation_factor: 2.0 -> 1.0)
- Does not directly increase collagen synthesis (correct per literature)

### Mesenchymal stem cell (MSC) therapy

Broadest-acting intervention via paracrine secretome:
- M1-to-M2 macrophage reprogramming (IL-6, PGE2)
- VEGF production (dominant paracrine factor, up to 100x from spheroids)
- Fibroblast activation (PDGF, TGF-beta)
- Keratinocyte proliferation/migration (EGF, bFGF, HGF)

### Moisture dressings

Maintains optimal wound hydration:
- 83% TEWL reduction (surface_loss_rate: 0.03 -> 0.005)
- Hydrogel moisture donation (recovery_rate: 0.02 -> 0.06)
- Exudate absorption clears inflammatory mediators

### Combination therapy

Rational multi-target approach addressing all four diabetic dysfunction axes:
1. **Immune**: anti-TNF-alpha (M1/M2 transition)
2. **Vascular**: HBO (O2/VEGF/angiogenesis)
3. **ECM**: doxycycline (MMP/collagen preservation)
4. **Moisture**: advanced dressing (hydration)

## Creating custom treatments

Create a new TOML file in `treatments/`:

```toml
# Treatment: my_therapy
# Literature: Author et al. Year
# Mechanism: what it does biologically

[skin.diabetic]
m1_duration_factor = 2.0    # override diabetic dysfunction
prolif_factor = 0.7          # partially restore proliferation

[skin.perfusion]
angio_rate = 0.005           # override vascular params
```

The treatment file only needs to contain the parameters it modifies. All other parameters remain at their diabetic profile values.

## Output

The treatment study script produces:
- `output/treatment_study/metrics_<treatment>.csv` contains full metrics for each run
- `output/treatment_study/treatment_comparison.csv` contains the summary comparison table
- Console output with closure rates, inflammation peaks, healing times

## Biological validity

Treatment parameter mappings are derived from published quantitative data:
- Clinical effect sizes (fold-changes, percentage improvements)
- In vitro mechanistic studies (MMP inhibition %, fibroblast proliferation fold-change)
- Animal model data (wound closure rates, histological measures)

See individual treatment TOML files for DOI-linked references.
