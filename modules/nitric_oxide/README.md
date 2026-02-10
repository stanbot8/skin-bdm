# Nitric Oxide Module

Nitric oxide produced by iNOS in M1 macrophages and neutrophils during
the inflammatory phase. Very short half-life (seconds to minutes),
modeled as high decay rate. Three functional roles: vasodilation
(increases local perfusion), antimicrobial (kills bacteria via reactive
nitrogen species), and anti-fibrotic (suppresses excessive collagen).

## Fields

| Field | ID | D | Decay | Unit |
|-------|-----|-----|-------|------|
| NitricOxide | 26 | 0.02 | 0.05 | normalized |

## Dynamics

- **iNOS production**: M1 macrophages and neutrophils produce NO via iNOS
  (Witte & Barbul 2002). M2 macrophages do not produce NO.
- **Vasodilation**: NO boosts effective perfusion target by
  (1 + factor * NO) (Luo & Chen 2005)
- **Antimicrobial**: NO suppresses biofilm growth via reactive nitrogen
  species (Fang 1997)
- **Anti-fibrotic**: NO reduces collagen deposition rate
  (Schaffer et al. 1996)
- **Diabetic impairment**: iNOS expression reduced ~40% in diabetic tissue
  (Witte et al. 2002)

## Key parameters

| Parameter | Default | Source |
|-----------|---------|--------|
| m1_production | 0.01 | Witte & Barbul 2002 |
| vasodilation_factor | 0.3 | Luo & Chen 2005 |
| antimicrobial_factor | 0.5 | Fang 1997 |
| collagen_suppression | 0.15 | Schaffer et al. 1996 |
