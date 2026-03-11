> [Home](../../README.md) / [Modules](../README.md) / Blood

# Blood

Hemorrhage, coagulation cascade (intrinsic and extrinsic pathways), blood volume tracking, and systemic effects on perfusion and oxygen delivery. Disabled by default; activated for hemorrhagic or surgical wound studies.

## Biology

Wound healing begins with hemostasis, the arrest of bleeding. When a vessel is injured, two coagulation pathways activate in parallel: the intrinsic pathway (contact activation of factor XII by exposed subendothelial collagen) and the extrinsic pathway (tissue factor exposure from damaged cells). Both converge on thrombin generation, which converts fibrinogen to fibrin and amplifies platelet aggregation to form a stable clot (Versteeg et al. 2013).

Platelets aggregating at the wound site release growth factors (PDGF, TGF-beta) from their alpha-granules, initiating the inflammatory cascade that recruits neutrophils and macrophages (Singer and Clark 1999). The provisional fibrin clot matures over approximately 4 hours and begins fibrinolytic remodeling after about 48 hours.

Significant blood loss reduces circulating volume and hemoglobin, impairing tissue perfusion and oxygen delivery to the wound bed. Below the hemorrhagic shock threshold (approximately 40% volume loss), perfusion drops sharply and wound healing is severely compromised (Guo and DiPietro 2010). Tissue oxygen tension is a direct predictor of wound infection risk (Hopf et al. 1997).

Clinical modifiers include anticoagulant therapy (which slows clotting but increases bleeding) and thrombocytopenia (low platelet count impairing clot formation).

## Model

**Three coupled subsystems:**

```
1. Hemorrhage
   bleed_rate * (1 - barrier) * vascularity * anticoagulant_boost
   -> reduces blood_volume over time
   -> volume_recovery partially compensates

2. Coagulation Cascade
   (intrinsic_rate + extrinsic_rate) * anticoag_penalty * platelet_factor
   -> coag_state (0 to 1, activation level)
   -> thrombin generation from converged pathways
   -> thrombin enhances fibrin deposition (hemostasis module coupling)
   -> platelet cytokine release (PDGF, TGF-beta) to inflammation

3. Blood Volume and Systemic Effects
   blood_volume tracks cumulative hemorrhage vs recovery
   Below shock_threshold: perfusion * shock_perfusion_penalty
   Below anemia_threshold: O2 delivery * anemia_o2_penalty
```

**Hemorrhagic shock:** When blood volume drops below `shock_threshold` (default 0.6), perfusion is progressively penalized down to `shock_perfusion_penalty` (0.3x), modeling circulatory failure.

**Anemia:** When hemoglobin drops below `anemia_threshold` (default 0.7), oxygen delivery is reduced by `anemia_o2_penalty`, impairing aerobic wound healing.

**Anticoagulant therapy:** Therapeutic anticoagulation (warfarin, heparin) reduces coagulation rates by `anticoag_coag_penalty` while increasing bleeding by `anticoag_bleed_boost`.

**Thrombocytopenia:** Platelet counts below `thrombocytopenia_threshold` proportionally reduce coagulation cascade activation.

## Parameters

From modules/blood/config.toml:

| Parameter | Default | Units | Description | Source |
|-----------|---------|-------|-------------|--------|
| `enabled` | false | bool | Master switch (opt-in) | Convention |
| `bleed_rate` | 0.02 | normalized | Baseline blood loss rate from wound | Calibrated |
| `depth_bleed_factor` | 1.5 | multiplier | Deeper wounds bleed more | Calibrated |
| `vascularity_coupling` | true | bool | Bleed rate scales with local perfusion | Convention |
| `intrinsic_rate` | 0.01 | per step | Contact activation (factor XII) | Versteeg et al. 2013 ([DOI](https://doi.org/10.1152/physrev.00016.2011)) |
| `extrinsic_rate` | 0.05 | per step | Tissue factor exposure from injury | Versteeg et al. 2013 ([DOI](https://doi.org/10.1152/physrev.00016.2011)) |
| `thrombin_rate` | 0.03 | per step | Thrombin generation from converged pathways | Versteeg et al. 2013 ([DOI](https://doi.org/10.1152/physrev.00016.2011)) |
| `thrombin_fibrin_coupling` | 0.02 | per step | Thrombin to fibrin conversion | Versteeg et al. 2013 ([DOI](https://doi.org/10.1152/physrev.00016.2011)) |
| `platelet_aggregation_rate` | 0.04 | per step | Platelet aggregation at wound site | Singer and Clark 1999 ([DOI](https://doi.org/10.1056/NEJM199909023411006)) |
| `platelet_cytokine_coupling` | 0.001 | per step | PDGF/TGF-beta release to inflammation | Singer and Clark 1999 ([DOI](https://doi.org/10.1056/NEJM199909023411006)) |
| `clot_maturation_h` | 4.0 | hours | Clot stabilization time | Calibrated |
| `fibrinolysis_delay_h` | 48.0 | hours | Fibrinolysis onset delay after clot formation | Calibrated |
| `initial_volume` | 1.0 | normalized | Normal blood volume | Convention |
| `volume_loss_rate` | 0.01 | per hour | Volume loss fraction at bleed_rate=1.0 | Calibrated |
| `volume_recovery_rate` | 0.005 | per hour | Physiological compensation and resuscitation | Guo and DiPietro 2010 ([DOI](https://doi.org/10.1177/0022034509359125)) |
| `shock_threshold` | 0.6 | fraction | Below this: hemorrhagic shock | Guo and DiPietro 2010 ([DOI](https://doi.org/10.1177/0022034509359125)) |
| `shock_perfusion_penalty` | 0.3 | multiplier | Perfusion reduction in shock | Guo and DiPietro 2010 ([DOI](https://doi.org/10.1177/0022034509359125)) |
| `hemoglobin` | 1.0 | normalized | Hemoglobin level (1.0 = 14 g/dL) | Convention |
| `anemia_threshold` | 0.7 | fraction | Below this: impaired O2 delivery | Hopf et al. 1997 ([DOI](https://doi.org/10.1001/archsurg.1997.01430330063010)) |
| `anemia_o2_penalty` | 0.5 | multiplier | O2 delivery reduction factor | Hopf et al. 1997 ([DOI](https://doi.org/10.1001/archsurg.1997.01430330063010)) |
| `anticoagulant_level` | 0.0 | 0 to 1 | Anticoagulant therapy level (0 = none, 1 = therapeutic) | Convention |
| `anticoag_coag_penalty` | 0.4 | multiplier | Coagulation rate reduction from therapy | Calibrated |
| `anticoag_bleed_boost` | 1.8 | multiplier | Bleed rate increase from therapy | Calibrated |
| `platelet_count` | 1.0 | normalized | Platelet count (1.0 = 250k/uL) | Convention |
| `thrombocytopenia_threshold` | 0.4 | fraction | Below this: impaired clotting | Calibrated |

## Coupling

### Reads
| Field | Source module | How used |
|-------|-------------|----------|
| Vascular (perfusion) | angiogenesis | Vascularity coupling for hemorrhage rate; perfusion modification target |
| Oxygen | oxygen | Anemia penalty applied to O2 delivery |
| Barrier | epithelial | Gates hemorrhage (open wound = bleeding) |

### Writes
| Field | Consumer modules | What is written |
|-------|-----------------|-----------------|
| Vascular (perfusion) | angiogenesis, oxygen | Reduced by hemorrhagic shock |
| Oxygen | metabolism | Reduced by anemia |
| Fibrin | hemostasis | Enhanced by thrombin generation |
| Inflammation | immune | Platelet cytokine release (PDGF, TGF-beta) |

## Source files

| File | Purpose |
|------|---------|
| `source_hook.h` | Hemorrhage, coagulation cascade, blood volume tracking, systemic perfusion and O2 effects |
| `params.h` | Parameter struct (BloodParams) |
| `config.toml` | Module configuration |
