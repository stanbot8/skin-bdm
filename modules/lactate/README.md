# Lactate Module

Tissue lactate produced by anaerobic glycolysis in hypoxic wound regions.
Wound lactate levels reach 4-12 mM, 3-5x higher than arterial blood.
Lactate acts as an angiogenesis signal by stabilizing HIF-1a (boosting
VEGF production) and directly stimulates collagen synthesis by fibroblasts.

## Fields

| Field | ID | D | Decay | Unit |
|-------|-----|-----|-------|------|
| Lactate | 25 | 0.01 | 0.002 | normalized |

## Dynamics

- **Hypoxia production**: produced where O2 < threshold (anaerobic glycolysis)
- **VEGF boost**: stabilizes HIF-1a, boosting VEGF production rate
  (Constant et al. 2000)
- **Collagen boost**: enhances fibroblast collagen synthesis ~2x
  (Hunt et al. 2007)
- **Perfusion clearance**: washed out by restored vascular perfusion

## Key parameters

| Parameter | Default | Source |
|-----------|---------|--------|
| o2_threshold | 0.3 | Gordillo & Sen 2003 |
| vegf_boost | 0.3 | Constant et al. 2000 |
| collagen_boost | 0.2 | Hunt et al. 2007 |
