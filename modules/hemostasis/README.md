> [Home](../../README.md) / [Modules](../README.md) / Hemostasis

# Hemostasis

Fibrin scaffold from platelet aggregation at the wound surface.

## Biology

Blood coagulation forms a fibrin clot within minutes of wounding. The clot
serves as a provisional matrix scaffold for keratinocyte and fibroblast
migration, and releases platelet-derived growth factors (PDGF, TGF-beta).
During remodeling, MMP-mediated fibrinolysis gradually replaces fibrin with
collagen-based tissue.

## Model

Fibrin PDE with D=0 (structural, immobile) and manual decay. Wound event
seeds fibrin at `wound_seed`. Coupling: fibrin presence produces TGF-beta
and fibronectin at configurable rates. MMP degrades fibrin (fibrinolysis).
Fibrin scaffold boosts keratinocyte migration speed.

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `enabled` | false | Master switch |
| `diffusion` | 0.0 | Immobile structural ECM |
| `decay` | 0.0003 | Slow turnover |
| `wound_seed` | 1.0 | Initial clot density |
| `tgfb_coupling` | 0.0005 | Fibrin to TGF-beta rate |
| `fibronectin_coupling` | 0.001 | Fibrin to fibronectin rate |
| `mmp_degradation` | 0.002 | MMP fibrinolysis rate |
| `migration_boost` | 0.3 | Keratinocyte migration boost |

## Source files

| File | Purpose |
|------|---------|
| `hemostasis_pde.h` | Fibrin field definition |
| `config.toml` | Module configuration |
