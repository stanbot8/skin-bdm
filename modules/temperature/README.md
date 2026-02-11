# Temperature Module

Wound bed thermal regulation. Healthy skin is maintained at core body
temperature (~37C) by vascular perfusion. Wounding exposes tissue to the
environment, cooling the wound surface to ~33C. Temperature affects cell
migration speed, enzyme kinetics (MMP activity), and bacterial growth via
Q10 coefficients.

## Fields

| Field | ID | D | Decay | Unit |
|-------|-----|-----|-------|------|
| Temperature | 23 | 0.001 | 0 | normalized (37C = 1.0) |

## Dynamics

- **Perfusion warming**: dermal voxels warm toward 1.0 proportional to local
  vascular perfusion (McGuiness et al. 2004)
- **Surface cooling**: exposed wound surface loses heat via evaporation and
  radiation (Fierheller & Sibbald 2010)
- **Q10 coupling**: cell migration, MMP activity, and biofilm growth rates
  scale as Q10^((T-37)/10)

## Key parameters

| Parameter | Default | Source |
|-----------|---------|--------|
| wound_surface | 0.892 | Fierheller & Sibbald 2010 |
| q10_migration | 2.0 | Kanokwan & Bhattacharya 2012 |
| q10_mmp | 1.5 | Fields 2001 |
| q10_biofilm | 2.5 | Ratkowsky et al. 1982 |
