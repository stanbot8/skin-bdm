# Treatment Interventions

Computational treatment studies modeling therapeutic interventions across six clinical domains. Each treatment is a TOML overlay modifying specific biological parameters based on published mechanisms of action.

30 treatments are available as TOML overlays in study-scoped `treatments/` directories (plus shared treatments in `studies/shared/treatments/`). For summary tables, see [studies/README.md#treatments](../studies/README.md#treatments).

## Running treatments

```bash
# Single treatment
./run.sh --study=diabetic-wound --skin=diabetic --treatment=hbo

# Full comparison study (all 9 treatments + baseline)
python3 scripts/study/treatment_study.py

# Specific treatments only
python3 scripts/study/treatment_study.py --treatments=hbo,msc,combination

# List available treatments
./run.sh --list-treatments
```

## How treatments work

Treatments are TOML overlays applied after the profile and study config. The search order is: study-scoped (`studies/{study}/treatments/`), then shared (`studies/shared/treatments/`), then all studies:

```
bdm.core.toml + modules/*/config.toml           (merge)
  + profiles/diabetic.toml                        (biology)
  + studies/diabetic-wound/preset.toml            (experiment)
  + studies/diabetic-wound/treatments/hbo.toml    (intervention)
  = bdm.toml                                      (runtime)
```

Each treatment modifies the diabetic dysfunction parameters back toward healthy values, proportional to the treatment's clinical effect size from published literature.

## Treatment mechanisms

### Anti-inflammatory (anti-TNF-alpha)

The core diabetic wound pathology is failed M1-to-M2 macrophage transition. TNF-alpha neutralization directly addresses this by:
- Releasing the M1-to-M2 brake (m1_duration_factor: 3.0 to 1.5)
- Restoring M2 anti-inflammatory resolution (resolution_factor: 0.3 to 0.7)
- Reducing AGE/RAGE-driven baseline inflammation (0.001 to 0.0002)
- Restoring keratinocyte proliferation and migration

### Hyperbaric oxygen (HBO)

HBO corrects the tissue hypoxia that drives much of diabetic wound dysfunction:
- Restores HIF-1alpha/VEGF cycling (vegf_factor: 0.4 to 1.0)
- Improves microcirculation (perfusion basal: 0.7 to 0.85)
- Accelerates angiogenesis (angio_rate: 0.002 to 0.008)
- Partially rescues collagen synthesis via prolyl hydroxylase O2 supply

### Negative pressure wound therapy (NPWT)

Multi-modal mechanical intervention:
- 4x blood flow increase from sub-atmospheric pressure
- Microdeformation stimulates cell proliferation
- Macrodeformation assists wound edge advancement (inward_bias: 0.3 to 0.5)
- Sealed dressing controls moisture and clears inflammatory exudate

### Doxycycline (sub-antimicrobial MMP inhibitor)

Targets the MMP/TIMP imbalance that destroys ECM in diabetic wounds:
- Halves excess MMP activity (mmp_factor: 3.0 to 1.5)
- Restores TIMP production (timp_production_factor: 0.4 to 0.76)
- Preserves collagen and fibronectin scaffolds
- Mild anti-inflammatory effect via TACE inhibition

### Growth factor (PDGF-BB / becaplermin)

FDA-approved recombinant growth factor for DFU:
- Potent fibroblast mitogen (density_factor: 1.0 to 1.8)
- Rapid fibroblast activation (activation_factor: 2.0 to 1.0)
- Does not directly increase collagen synthesis (correct per literature)

### Mesenchymal stem cell (MSC) therapy

Broadest-acting intervention via paracrine secretome:
- M1-to-M2 macrophage reprogramming (IL-6, PGE2)
- VEGF production (dominant paracrine factor, up to 100x from spheroids)
- Fibroblast activation (PDGF, TGF-beta)
- Keratinocyte proliferation/migration (EGF, bFGF, HGF)

### Moisture dressings

Maintains optimal wound hydration:
- 83% TEWL reduction (surface_loss_rate: 0.03 to 0.005)
- Hydrogel moisture donation (recovery_rate: 0.02 to 0.06)
- Exudate absorption clears inflammatory mediators

### Senolytic therapy (dasatinib + quercetin)

Targeted clearance of senescent cells that accumulate in chronic diabetic wounds:
- Accelerated senescent cell apoptosis (senolytic_clearance: 0.0 to 0.01)
- Reduced SASP burden (inflammation, MMP, TGF-beta output from senescent cells)
- Restored proliferative capacity in surrounding tissue
- Lower baseline inflammation from SASP-driven chronic signaling

### Combination therapy

Rational multi-target approach addressing all four diabetic dysfunction axes:
1. **Immune**: anti-TNF-alpha (M1/M2 transition)
2. **Vascular**: HBO (O2/VEGF/angiogenesis)
3. **ECM**: doxycycline (MMP/collagen preservation)
4. **Moisture**: advanced dressing (hydration)

## Treatment scheduling

The schedule generator creates experiment TOMLs that model delayed treatment start. A treatment started at day 0 has full effect; later starts have proportionally reduced effect based on the remaining healing window (42 days default).

```bash
# Single treatment schedules (NPWT/HBO/MSC at days 0, 7, 14, 21)
python3 scripts/study/gen_treatment_schedule.py --treatments npwt,hbo,msc --days 0,7,14,21

# Pairwise combinatorial scheduling (cross-product of start days)
python3 scripts/study/gen_treatment_schedule.py --treatments npwt,hbo,msc --days 0,7,14,21 --combo

# Run all generated schedules
python3 studies/run_experiments.py studies/diabetic-wound/experiments/schedule_*.toml
```

Output goes to `studies/diabetic-wound/experiments/schedule_*.toml`. Each experiment includes an untreated diabetic baseline and configs for each start day, with parameters interpolated between the diabetic profile and full treatment values.

## Creating custom treatments

Create a new TOML file in the study's `treatments/` directory (or `studies/shared/treatments/` for cross-study treatments):

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

## RA treatment mechanisms

Five biologic and DMARD treatments target the dual TNF-alpha/IL-6 inflammatory axis in rheumatoid arthritis.

### Anti-TNF (infliximab/adalimumab/etanercept)

Monoclonal antibody neutralization of soluble and membrane-bound TNF-alpha:
- Direct TNF clearance (anti_tnf_clearance: 0.0 to 0.05)
- Reduced NF-kB inflammatory amplification (tnf_inflammation_coupling: 0.01 to 0.004)
- Lower MMP transcription (tnf_mmp_boost: 0.008 to 0.003)
- Pannus regression via VEGF reduction (tnf_vegf_boost: 0.005 to 0.002)
- Partial M1 duration normalization (m1_prolongation: 2.0 to 1.4)

### Tocilizumab (anti-IL-6R)

Humanized monoclonal antibody blocking IL-6 receptor signaling:
- IL-6R blockade (anti_il6r_clearance: 0.0 to 0.04)
- Suppressed JAK/STAT3 amplification (il6_inflammation_coupling: 0.008 to 0.003)
- Reduced RANKL-mediated bone erosion (il6_cartilage_boost: 0.002 to 0.0008)

### Methotrexate (conventional DMARD)

Folate antagonist with broad immunosuppressive effects:
- Reduced autoimmune TNF/IL-6 source (autoimmune_source halved)
- Lower NF-kB and JAK/STAT3 amplification
- Moderate MMP reduction (tnf_mmp_boost: 0.008 to 0.005)
- Reduced FLS hyperplasia (pannus_fibroblast_boost: 2.0 to 1.5)

### JAK inhibitor (tofacitinib/baricitinib)

Small molecule inhibitor blocking JAK1/JAK3 intracellular signaling:
- Blocks JAK/STAT pathway downstream of multiple cytokine receptors
- Suppressed T cell proliferation and recruitment
- Reduced IL-6 amplification via STAT3 blockade
- Oral small molecule (faster onset than biologics)

### RA triple combination

Combined anti-TNF + tocilizumab + methotrexate:
- All three clearance and suppression mechanisms active simultaneously
- Near-normal M1 duration and FLS proliferation
- Strongest cartilage protection (all three erosion pathways suppressed)

Note: dual biologic therapy (anti-TNF + tocilizumab) is not standard clinical practice due to infection risk. This combination is modeled for mechanistic study of maximal cytokine blockade.

## Burn treatments

Five treatments in `studies/burn/treatments/` targeting thermal injury:

| Treatment | File | Mechanism |
|-----------|------|-----------|
| Cooling | `cooling.toml` | Immediate cold water first aid, preserves stasis zone |
| Debridement | `debridement.toml` | Surgical eschar removal, exposes viable wound bed |
| Silver sulfadiazine | `silver_sulfadiazine.toml` | Broad-spectrum topical antimicrobial (mildly cytotoxic to keratinocytes) |
| Skin substitute | `skin_substitute.toml` | Bioengineered collagen-GAG scaffold (Integra/Biobrane), reduces TEWL |
| Pressure garment | `pressure_garment.toml` | Compression therapy (15 to 25 mmHg), prevents hypertrophic scar |

## Pressure ulcer treatments

Five treatments in `studies/pressure-ulcer/treatments/` targeting ischemia-reperfusion injury:

| Treatment | File | Mechanism |
|-----------|------|-----------|
| Pressure redistribution | `pressure_redistribution.toml` | Low-air-loss mattress with 2h repositioning protocol |
| Wound VAC | `wound_vac.toml` | Negative pressure wound therapy for Stage III/IV |
| Nutrition | `nutrition.toml` | Protein, zinc, vitamin C supplementation for tissue repair |
| Silver dressing | `silver_dressing.toml` | Antimicrobial silver with moisture balance |
| Offloading | `offloading_protocol.toml` | Complete pressure elimination (gold standard prevention) |

## Surgical treatments

Five treatments in `studies/surgical/treatments/` targeting SSI prevention:

| Treatment | File | Mechanism |
|-----------|------|-----------|
| Prophylactic antibiotics | `prophylactic_antibiotics.toml` | Perioperative antimicrobial prophylaxis |
| Chlorhexidine | `chlorhexidine.toml` | Antiseptic skin preparation |
| Negative pressure | `negative_pressure.toml` | Incisional NPWT over closed surgical wounds |
| Enhanced recovery | `enhanced_recovery.toml` | ERAS protocol: nutrition, perfusion, reduced inflammation |
| Antimicrobial suture | `antimicrobial_suture.toml` | Triclosan-coated suture material |

## Shared treatments

Cross-study treatments in `studies/shared/treatments/`, searchable from any study:

| Treatment | File | Mechanism |
|-----------|------|-----------|
| Triple antibiotic | `triple_antibiotic.toml` | Bacitracin/neomycin/polymyxin B in petrolatum base; reduces infection, maintains moist wound environment |
