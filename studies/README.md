# Studies

Computational studies and literature-validated simulation experiments.

## Running studies

Each `.sh` script in this directory is self-contained: it builds (if needed),
runs simulations, collects metrics, and generates output.

```bash
source ~/biodynamo/build/bin/thisbdm.sh

# Diabetic treatment comparison (baseline + 8 treatments, ~3-4h)
./studies/diabetic-study.sh

# Adaptive combo search (surrogate-guided, ~30-45min)
./studies/adaptive-study.sh

# Skin type comparison (normal vs aged vs diabetic)
./studies/skin-comparison-study.sh

# Tumor growth kinetics
./studies/tumor-study.sh

# Novel experiment scenarios (7 experiments)
./studies/run-scenarios.sh
```

## Study scripts

| Script | Python backend | Description |
|--------|---------------|-------------|
| [diabetic-study.sh](diabetic-study.sh) | [treatment_study.py](../scripts/study/treatment_study.py) | Baseline + 8 treatments, Excel workbook |
| [adaptive-study.sh](adaptive-study.sh) | [adaptive_study.py](../scripts/study/adaptive_study.py) | Surrogate-guided combo search with synergy detection |
| [skin-comparison-study.sh](skin-comparison-study.sh) | [treatment_study.py](../scripts/study/treatment_study.py) | Normal / aged / diabetic skin profiles |
| [tumor-study.sh](tumor-study.sh) | [treatment_study.py](../scripts/study/treatment_study.py) | BCC/SCC growth rate validation |
| [run-scenarios.sh](run-scenarios.sh) | [scenario_runner.py](../scripts/study/scenario_runner.py) | Novel experiment scenarios (7 TOMLs) |

## Experiment scenarios

Seven TOML-defined experiments in [`scenarios/`](../scenarios/), each exploring a different axis of diabetic wound healing:

| Scenario | File | Configs | Runs |
|----------|------|---------|------|
| Severity spectrum | [severity_spectrum.toml](../scenarios/severity_spectrum.toml) | 5 levels (0.5x to 1.5x) | 5/config |
| Biofilm infection | [biofilm_infection.toml](../scenarios/biofilm_infection.toml) | 4 (control, early/late, +doxy) | 5/config |
| Treatment timing | [treatment_timing.toml](../scenarios/treatment_timing.toml) | 12 (3 treatments x 4 windows) | 5/config |
| Wound size | [wound_size.toml](../scenarios/wound_size.toml) | 4 (1.5mm to 8mm) | 5/config |
| Aged + diabetic | [aged_diabetic.toml](../scenarios/aged_diabetic.toml) | 4 phenotypes | 5/config |
| Chronic wound (90d) | [chronic_wound.toml](../scenarios/chronic_wound.toml) | 4 (severe + biofilm + rescue) | 3/config |
| Immune tuning | [immune_tuning.toml](../scenarios/immune_tuning.toml) | 4 (restore one axis) | 5/config |

```bash
bash studies/run-scenarios.sh                                # all scenarios
bash studies/run-scenarios.sh scenarios/wound_size.toml      # single scenario
bash studies/run-scenarios.sh --quick                        # all with 2 runs
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
| Cellular senescence | SASP-driven chronic inflammation | Demaria 2014, Wilkinson 2019 | High |
| Oxidative stress | ROS-mediated tissue damage | Schafer 2008, Sen 2009 | High |
| Basement membrane | Laminin/collagen IV reassembly | Rousselle 2019 | High |

## Validation pipeline

```bash
# Run literature validation on any metrics CSV
python3 literature/validate_all.py output/skibidy/metrics.csv

# Check source integrity (DOIs referenced in configs vs SOURCES.yaml)
python3 literature/check_sources.py

# 10-run batch consensus with validation
python3 batch/batch.py -n 10 --skin normal --preset wound --validate
```

## Example outputs

Pre-generated results are included so you can inspect output without running simulations.

| Study | Example file | Description |
|-------|-------------|-------------|
| Diabetic treatments | [diabetic-example/](diabetic-example/) | 8-treatment comparison Excel workbook |
| Adaptive combos | [adaptive-example/](adaptive-example/) | Surrogate predictions + synergy analysis |
