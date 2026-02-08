# Glucose Module

Interstitial glucose field representing metabolic substrate availability.
In healthy tissue, glucose is supplied by vascular perfusion at a normalized
level of 1.0. In diabetic tissue, hyperglycemia elevates the basal level
to ~1.8, driving AGE (advanced glycation end-product) formation that
activates RAGE receptors and sustains NF-kB mediated inflammation.

## Fields

| Field | ID | D | Decay | Unit |
|-------|-----|-----|-------|------|
| Glucose | 24 | 0.005 | 0.001 | normalized (fasting = 1.0) |

## Dynamics

- **Perfusion supply**: dermal voxels receive glucose proportional to local
  vascular perfusion (Casciari et al. 1992)
- **Cellular uptake**: background decay represents metabolic consumption
- **AGE formation**: in diabetic mode, elevated glucose drives AGE formation
  which triggers RAGE mediated inflammation (Bierhaus et al. 2005)
- **Bacterial consumption**: biofilm consumes glucose (Bowler et al. 2001)

## Key parameters

| Parameter | Default | Source |
|-----------|---------|--------|
| basal_conc | 1.0 | normalized |
| age_inflammation | 0.002 | Bierhaus et al. 2005 |
| prolif_threshold | 0.2 | estimated |
