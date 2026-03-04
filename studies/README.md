# Studies

Computational studies and literature-validated simulation experiments.

## Running studies

Each study lives in its own subdirectory with a `preset.toml`, scripts, experiments, and results.

```bash
source ~/biodynamo/build/bin/thisbdm.sh

# Diabetic treatment comparison (baseline + 9 treatments)
python3 studies/diabetic-wound/treatment.py

# Adaptive combo search (surrogate-guided)
python3 studies/diabetic-wound/adaptive.py

# Skin type comparison (normal vs aged vs diabetic)
python3 studies/wound/skin_comparison.py

# Tumor growth kinetics
python3 studies/tumor/study.py

# RA treatment comparison (baseline + 5 treatments)
python3 studies/rheumatoid/study.py

# Run all experiments
python3 studies/run_experiments.py
```

## Study scripts

| Script | Backend | Description |
|--------|---------|-------------|
| [diabetic-wound/treatment.py](diabetic-wound/treatment.py) | [treatment_study.py](../scripts/study/treatment_study.py) | Baseline + 9 treatments, Excel workbook |
| [diabetic-wound/adaptive.py](diabetic-wound/adaptive.py) | [adaptive_study.py](../scripts/study/adaptive_study.py) | Surrogate-guided combo search with synergy detection |
| [wound/skin_comparison.py](wound/skin_comparison.py) | batch/batch.py | Normal / aged / diabetic skin profiles |
| [tumor/study.py](tumor/study.py) | batch/batch.py | BCC/SCC growth rate validation |
| [rheumatoid/study.py](rheumatoid/study.py) | batch/batch.py | RA treatment comparison (5 treatments) |
| [run_experiments.py](run_experiments.py) | [experiment_runner.py](../scripts/study/experiment_runner.py) | Run experiments across all studies |

## Experiments

### Diabetic wound experiments

TOML-defined experiments in `studies/diabetic-wound/experiments/`, each exploring a different axis of diabetic wound healing:

| Experiment | File | Configs | Runs |
|------------|------|---------|------|
| Severity spectrum | [severity_spectrum.toml](diabetic-wound/experiments/severity_spectrum.toml) | 5 levels (0.5x to 1.5x) | 5/config |
| Biofilm infection | [biofilm_infection.toml](diabetic-wound/experiments/biofilm_infection.toml) | 4 (control, early/late, +doxy) | 5/config |
| Treatment timing | [treatment_timing.toml](diabetic-wound/experiments/treatment_timing.toml) | 12 (3 treatments x 4 windows) | 5/config |
| Wound size | [wound_size.toml](diabetic-wound/experiments/wound_size.toml) | 4 (1.5mm to 8mm) | 5/config |
| Aged + diabetic | [aged_diabetic.toml](diabetic-wound/experiments/aged_diabetic.toml) | 4 phenotypes | 5/config |
| Chronic wound (90d) | [chronic_wound.toml](diabetic-wound/experiments/chronic_wound.toml) | 4 (severe + biofilm + rescue) | 3/config |
| Immune tuning | [immune_tuning.toml](diabetic-wound/experiments/immune_tuning.toml) | 4 (restore one axis) | 5/config |
| Cure protocol | [cure_protocol.toml](diabetic-wound/experiments/cure_protocol.toml) | 5 (untreated through 5-agent cure) | 5/config |

### Rheumatoid arthritis experiments

TOML-defined experiments in `studies/rheumatoid/experiments/`, exploring RA treatment strategies and disease progression:

| Experiment | File | Configs | Runs |
|------------|------|---------|------|
| Biologic comparison | [biologic_comparison.toml](rheumatoid/experiments/biologic_comparison.toml) | 5 (untreated + 4 treatments) | 5/config |
| Disease severity | [disease_severity.toml](rheumatoid/experiments/disease_severity.toml) | 5 (0.5x to 2.0x flare intensity) | 5/config |
| Treatment timing | [treatment_timing.toml](rheumatoid/experiments/treatment_timing.toml) | 8 (anti-TNF + tocilizumab x 4 windows) | 5/config |
| Cartilage protection | [cartilage_protection.toml](rheumatoid/experiments/cartilage_protection.toml) | 6 (untreated + 5 strategies) | 5/config |
| Cure protocol | [cure_protocol.toml](rheumatoid/experiments/cure_protocol.toml) | 5 (mono through triple+JAKi) | 5/config |

### Burn experiments

TOML-defined experiments in `studies/burn/experiments/`:

| Experiment | File | Configs | Runs |
|------------|------|---------|------|
| Depth spectrum | [depth_spectrum.toml](burn/experiments/depth_spectrum.toml) | 4 (superficial to full thickness) | 5/config |
| Stasis rescue | [stasis_rescue.toml](burn/experiments/stasis_rescue.toml) | 4 (untreated vs cooling at 0h, 6h, 24h) | 5/config |
| Treatment comparison | [treatment_comparison.toml](burn/experiments/treatment_comparison.toml) | 6 (untreated through full protocol) | 5/config |
| Scar prevention | [scar_prevention.toml](burn/experiments/scar_prevention.toml) | 3 (no prevention, pressure, pressure + silicone) | 5/config |
| Cure protocol | [cure_protocol.toml](burn/experiments/cure_protocol.toml) | 4 (untreated through full + antimicrobial) | 5/config |

### Pressure ulcer experiments

TOML-defined experiments in `studies/pressure-ulcer/experiments/`:

| Experiment | File | Configs | Runs |
|------------|------|---------|------|
| Stage progression | [stage_progression.toml](pressure-ulcer/experiments/stage_progression.toml) | 4 (Stage I through IV) | 5/config |
| Prevention vs treatment | [prevention_vs_treatment.toml](pressure-ulcer/experiments/prevention_vs_treatment.toml) | 3 (none, early, late) | 5/config |
| Treatment comparison | [treatment_comparison.toml](pressure-ulcer/experiments/treatment_comparison.toml) | 8 (untreated + 5 mono + 2 combos + all) | 5/config |
| Moisture and shear | [moisture_shear.toml](pressure-ulcer/experiments/moisture_shear.toml) | 4 (2x2 factorial) | 5/config |
| Cure protocol | [cure_protocol.toml](pressure-ulcer/experiments/cure_protocol.toml) | 4 (untreated through full + perfusion rescue) | 5/config |

### Surgical experiments

TOML-defined experiments in `studies/surgical/experiments/`:

| Experiment | File | Configs | Runs |
|------------|------|---------|------|
| Infection risk | [infection_risk.toml](surgical/experiments/infection_risk.toml) | 4 (CDC wound classes) | 5/config |
| Prevention bundle | [prevention_bundle.toml](surgical/experiments/prevention_bundle.toml) | 4 (additive SSI prevention) | 5/config |
| Comorbidity risk | [comorbidity_risk.toml](surgical/experiments/comorbidity_risk.toml) | 4 (obesity, diabetes, immunosuppression, age) | 5/config |
| Timing optimization | [timing_optimization.toml](surgical/experiments/timing_optimization.toml) | 4 (antibiotic prophylaxis windows) | 5/config |

```bash
python3 studies/run_experiments.py                                               # all experiments
python3 studies/run_experiments.py studies/diabetic-wound/experiments/wound_size.toml   # single experiment
python3 studies/run_experiments.py --quick                                       # all with 2 runs
```

## Literature-validated observables

Timeseries validated against published consensus curves (RMSE < 15%).

| Observable | Consensus CSV | Key sources | Condition |
|-----------|--------------|-------------|-----------|
| Wound closure | [closure_kinetics_punch_biopsy.csv](../modules/wound/data/closure_kinetics_punch_biopsy.csv) | Cukjati 2001, Gonzalez 2016 | Normal |
| Inflammation | [inflammation_timecourse.csv](../modules/inflammation/data/inflammation_timecourse.csv) | Eming 2007, Koh 2011 | Normal |
| Immune cells | [immune_cell_kinetics.csv](../modules/immune/data/immune_cell_kinetics.csv) | Kim 2008, Rodero 2010 | Normal |
| Myofibroblasts | [myofibroblast_kinetics.csv](../modules/fibroblast/data/myofibroblast_kinetics.csv) | Darby 2014, Desmouliere 1995 | Normal |
| Collagen | [collagen_deposition.csv](../modules/fibroblast/data/collagen_deposition.csv) | Zhou 2013, Murphy 2012 | Normal |
| Diabetic closure | [diabetic_closure_kinetics.csv](../modules/diabetic/data/diabetic_closure_kinetics.csv) | Mirza 2011, Louiselle 2021 | Diabetic |
| Diabetic inflammation | [diabetic_inflammation_timecourse.csv](../modules/diabetic/data/diabetic_inflammation_timecourse.csv) | Wetzler 2000, Mirza 2011 | Diabetic |
| Diabetic immune cells | [diabetic_immune_cell_kinetics.csv](../modules/diabetic/data/diabetic_immune_cell_kinetics.csv) | Khanna 2010, Wang 2020 | Diabetic |
| Tumor growth | [tumor_growth_rate.csv](../modules/tumor/data/tumor_growth_rate.csv) | Kricker 2014, Fijalkowska 2023, Sykes 2020 | Tumor |

## Parameter-validated modules

Modules with parameter values sourced from literature (not timeseries-validated).

| Module | Status | Parameters | Key sources |
|--------|--------|-----------|-------------|
| Angiogenesis | Enabled | VEGF diffusion, sprout rate, capillary density | Schugart 2008 ([10.1016/j.jtbi.2008.06.042](https://doi.org/10.1016/j.jtbi.2008.06.042)), Flegg 2012 |
| MMP | Enabled | Decay, production rates, TIMP interaction | Nagase 1999, Ladwig 2002 ([10.1067/mjd.2002.124601](https://doi.org/10.1067/mjd.2002.124601)) |
| Fibronectin | Enabled | Deposition rate, degradation | Clark 1990, Grinnell 1984 |
| Scar | Enabled | Collagen-based scar scoring, remodeling | Ogawa 2017 ([10.3390/ijms18030606](https://doi.org/10.3390/ijms18030606)), Gauglitz 2011 |
| Hemostasis | Disabled | Platelet plug, fibrin mesh, clotting cascade | Brass 2010, Reininger 2006 |
| Perfusion | Enabled | O2 transport, vascular delivery | Johnson 1971, Stucker 2002 |
| Dermis | Enabled | Dermal ECM, thickness, hydration | Braverman 2000, Singer 1999 |

## Treatments

### Diabetic wound treatments

| Treatment | Mechanism | Key sources | Parameters modified |
|-----------|----------|-------------|-------------------|
| Anti-inflammatory | Anti-TNF-alpha, accelerates M1-to-M2 | Goren 2007 | M1 decay, M2 transition |
| Combination | Multi-modal (anti-infl + HBO + doxy + moisture) | Composite | Multiple |
| Doxycycline | Sub-antimicrobial MMP inhibitor | Siqueira 2010, Smith 1999 | MMP production, collagen decay |
| Growth factor | Becaplermin (PDGF-BB) | Steed 2006, Smiell 1998 | KGF rate, chemotaxis |
| HBO | Hyperbaric oxygen therapy | Londahl 2010, Fedorko 2016 | O2 delivery, VEGF |
| Moisture | Advanced dressings (hydrogel/foam) | Junker 2013, Kannon 1995 | Water recovery, evaporation |
| MSC | Mesenchymal stem cell therapy | Cao 2017 | Immune modulation, growth factors |
| NPWT | Negative pressure wound therapy | Morykwas 1997, Armstrong 2005 | Perfusion, granulation |
| Senolytic | Dasatinib + quercetin senescent cell clearance | Zhu 2015, Hickson 2019 | Senolytic clearance, SASP reduction |

### Burn treatments

| Treatment | Mechanism | Key sources | Parameters modified |
|-----------|----------|-------------|-------------------|
| Cooling | Immediate cold water first aid | Jackson 1953 | Stasis preservation, inflammation |
| Debridement | Surgical eschar removal | Atiyeh 2005 | Eschar clearance, wound bed exposure |
| Silver sulfadiazine | Topical antimicrobial cream | Fox 1968, Atiyeh 2007 | Infection, mild keratinocyte toxicity |
| Skin substitute | Collagen-GAG scaffold (Integra) | Burke 1981, Heimbach 2003 | TEWL barrier, vascular ingrowth, migration |
| Pressure garment | Compression (15 to 25 mmHg) | Anzarut 2009, Engrav 2010 | Contracture, fibroblast apoptosis, collagen alignment |

### Pressure ulcer treatments

| Treatment | Mechanism | Key sources | Parameters modified |
|-----------|----------|-------------|-------------------|
| Pressure redistribution | Low-air-loss mattress + 2h repositioning | NPUAP 2019 | Ischemia, tissue damage, shear |
| Wound VAC | Negative pressure wound therapy | Morykwas 1997 | Tissue damage, moisture, granulation |
| Nutrition | Protein, zinc, vitamin C | Stechmiller 2010 | Proliferation, collagen, fibroblast density |
| Silver dressing | Antimicrobial + moisture balance | Lansdown 2002 | Infection, inflammation |
| Offloading | Complete pressure elimination | NPUAP 2019 | Compression, ischemia, migration |

### Surgical treatments

| Treatment | Mechanism | Key sources | Parameters modified |
|-----------|----------|-------------|-------------------|
| Prophylactic antibiotics | Perioperative antimicrobial prophylaxis | Bratzler 2013 | Infection, inflammation |
| Chlorhexidine | Antiseptic skin preparation | Darouiche 2010 | Biofilm, infection risk |
| Negative pressure | Incisional NPWT | Webster 2019 | Tissue damage, moisture, perfusion |
| Enhanced recovery | ERAS protocol (nutrition + perfusion) | Ljungqvist 2017 | Proliferation, collagen, perfusion |
| Antimicrobial suture | Triclosan-coated suture | Leaper 2017 | Infection, local inflammation |

### Shared treatments

Cross-study treatments in `studies/shared/treatments/`, usable with any study via `--treatment`:

| Treatment | Mechanism | Parameters modified |
|-----------|----------|-------------------|
| Triple antibiotic | Bacitracin/neomycin/polymyxin B ointment | Infection, moisture, proliferation, migration |

### RA treatments

| Treatment | Mechanism | Key sources | Parameters modified |
|-----------|----------|-------------|-------------------|
| Anti-TNF | Infliximab/adalimumab TNF-alpha neutralization | Feldmann & Maini 2003 | TNF clearance, NF-kB, MMP, M1 duration |
| Tocilizumab | Anti-IL-6R JAK/STAT3 blockade | Kishimoto 2005 | IL-6R clearance, RANKL bone erosion |
| Methotrexate | Folate antagonist DMARD | Smolen 2018 | Autoimmune source, FLS hyperplasia |
| JAK inhibitor | Tofacitinib/baricitinib JAK1/3 blockade | O'Shea 2013, Fleischmann 2017 | JAK/STAT signaling, T cell proliferation |
| RA combination | Anti-TNF + tocilizumab + methotrexate | Smolen 2018, McInnes 2011 | All three axes |

## Coded but disabled

Modules with implementation complete but disabled by default (awaiting calibration or validation data).

| Module | Config key | Sources | Notes |
|--------|-----------|---------|-------|
| Biofilm | `skin.biofilm.enabled` | James 2008, Bjarnsholt 2008, Davis 2008 | Bacterial colonization, immune evasion |
| pH | `skin.ph.enabled` | Schneider 2007, Gethin 2007 | Wound bed pH gradient, enzyme activity |
| Hyaluronan | `skin.hyaluronan.enabled` | Toole 2004, Stern 2006 | HAS2-driven HA synthesis, hydration |
| Elastin | `skin.elastin.enabled` | Kielty 2002, Almine 2012 | Tropoelastin deposition, cross-linking |
| Hemostasis | `skin.hemostasis.enabled` | Brass 2010, Reininger 2006 | Platelet activation, fibrin scaffold |

## Potential expansions

Areas with published mechanistic data that could extend the model.

| Feature | Biological basis | Candidate sources | Complexity |
|---------|-----------------|------------------|------------|
| Temperature therapy | Accelerated enzymatic rates, vasodilation | Ikeda 2005, Kloth 2002 | Low |
| Oxygen therapy (topical) | Direct O2 application vs HBO | Gordillo 2007 | Low |
| pH treatment | Acidic dressings shift enzyme optimum | Schneider 2007 | Medium |
| HA treatment | Exogenous hyaluronan scaffolds | Tolg 2014, Voigt 2012 | Medium |
| Basement membrane | Laminin/collagen IV reassembly | Rousselle 2019 | High |

## Mechanistic test study

The `full-model` study enables mechanistic toggles for validation testing. These replace simplified parametric models with biophysically grounded alternatives (Michaelis-Menten kinetics, HIF-1alpha stabilization, efferocytosis counting).

```bash
# 10-run consensus with mechanistic toggles
python3 batch/batch.py -n 10 --skin normal --study full-model --validate
```

Enabled toggles: M1 to M2 efferocytosis count, constitutive + TGF-beta collagen, HIF-1alpha VEGF. See [docs/architecture.md](../docs/architecture.md#mechanistic-toggles) for details.

## Validation pipeline

```bash
# Run literature validation on any metrics CSV
python3 literature/validate_all.py output/skibidy/metrics.csv

# Check source integrity (DOIs referenced in configs vs SOURCES.yaml)
python3 literature/check_sources.py

# 10-run batch consensus with validation
python3 batch/batch.py -n 10 --skin normal --study wound --validate
python3 batch/batch.py -n 10 --skin diabetic --study diabetic-wound --validate
```

## Example outputs

Pre-generated results are included so you can inspect output without running simulations.

| Study | Example file | Description |
|-------|-------------|-------------|
| Diabetic treatments | [diabetic-wound/example/](diabetic-wound/example/) | 8-treatment comparison Excel workbook |
| Adaptive combos | [diabetic-wound/adaptive-example/](diabetic-wound/adaptive-example/) | Surrogate predictions + synergy analysis |
