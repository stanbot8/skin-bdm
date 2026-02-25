# Validation

Compares simulation output against digitized literature data and checks source integrity. Reference data lives in each module under `modules/<module>/data/`, with citations in per-module `SOURCES.yaml` files.

## Running

```bash
# Full pipeline (source check + sim validation + plots)
python3 literature/validate_all.py output/skibidy/metrics.csv

# Source integrity only (no sim data needed)
python3 literature/check_sources.py

# 10-run batch consensus with validation
python3 batch/batch.py -n 10 --preset wound --validate
```

`validate_all.py` runs the source check first, then loads simulation metrics and validates whichever modules are present (wound, fibroblast, tumor).

## Scripts

| Script | Purpose |
|--------|---------|
| `validate_all.py` | Single entry point: source check, compute, print, plot |
| `check_sources.py` | SOURCES.yaml structural integrity + DOI cross-check |
| `lib.py` | Shared utilities: CSV loading, interpolation, RMSE, plotting |

## Source integrity checks

`check_sources.py` parses per-module `SOURCES.yaml` files and validates:

| Check | Level | What |
|-------|-------|------|
| YAML parse | ERROR | File loads cleanly |
| `description` field | ERROR | Every dataset has a description |
| `sources` list | ERROR | Every dataset has at least one citation |
| Citation fields | ERROR | Each source has `id`, `authors`, `year`, `title` |
| File references | ERROR | `consensus` and `raw_files` paths resolve to existing files |
| DOI/PMC/URL | WARN | At least one locator per source |
| Config DOI cross-check | WARN | Inline `doi:10.xxx` comments in configs match SOURCES.yaml |

## Output

Plots are saved to `output/plots/`:

| File | Contents |
|------|----------|
| `validation_dashboard.png` | Combined multi-panel dashboard |
| `wound_validation.png` | Closure, inflammation, immune cells, stratification |
| `fibroblast_validation.png` | Myofibroblast count + collagen accumulation |
| `tumor_validation.png` | Growth rate + doubling time |

For detailed parameter sources and validation datasets, see each module's README under `modules/*/README.md`.
