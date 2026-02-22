> [Home](../README.md) / Batch

# Batch System

Multi-run consensus and parameter sweeps for Skibidy. Self-contained: all scripts, configs, and results live here.

## Multi-run consensus

Run N simulations with the same config but different random seeds, compute mean/std, and validate against literature:

```bash
python3 batch/batch.py                         # 20 runs, default config
python3 batch/batch.py -n 5                    # 5 runs (faster)
python3 batch/batch.py -n 10 --preset wound    # wound preset, 10 runs
python3 batch/batch.py --skin aged --preset wound --validate
```

Output: `batch/results/consensus_<skin>_<preset>_<timestamp>/`
- `raw/run_NNN.csv`: individual run metrics
- `metrics.csv`: mean values (validation-compatible)
- `consensus.csv`: mean + std for each column

## Parameter sweeps

Sweep one or two parameters across a range of values with replicates at each point:

```bash
python3 batch/sweep.py batch/configs/cytokine_rate.toml       # single-param sweep
python3 batch/sweep.py batch/configs/collagen_vs_scar.toml     # two-param heatmap
python3 batch/sweep.py batch/configs/cytokine_rate.toml -n 3   # override replicates
python3 batch/sweep.py batch/configs/diabetic_factors.toml --analyze
```

Output: `batch/results/<name>_<timestamp>/`
- `config.toml`: copy of sweep config
- `raw/ptNNN_val*_runNNN.csv`: individual run metrics
- `summary.csv`: param values + outcome mean/std at each point
- `sweep.png`: auto-generated plot (single-param with error bars)
- `heatmap.png`: auto-generated plot (two-param colored grid)

## Post-hoc analysis

Analyze any results directory after the fact:

```bash
python3 batch/analyze.py batch/results/cytokine_rate_20260225_143000/
python3 batch/analyze.py batch/results/consensus_default_wound_20260225/ --validate
```

## Sweep config format

Configs are TOML files in `batch/configs/`. Single parameter:

```toml
[sweep]
name = "cytokine_rate"
param = "skin.immune.cytokine_rate"
values = [0.0005, 0.001, 0.0015, 0.002, 0.003]
runs_per_value = 5
preset = "wound"          # optional

[outcomes]
primary = "wound_closure_pct"
secondary = ["mean_infl_wound", "n_macrophages"]
measure = "final"         # final, peak, auc, time_to_90
```

Two parameters (produces NxM grid):

```toml
[sweep]
name = "collagen_vs_scar"
runs_per_value = 3

[[sweep.params]]
param = "skin.fibroblast.collagen_deposition_rate"
values = [0.004, 0.008, 0.012, 0.016]

[[sweep.params]]
param = "skin.scar.collagen_threshold"
values = [0.05, 0.1, 0.2, 0.5]

[outcomes]
primary = "scar_magnitude"
measure = "final"
```

## Included configs

| Config | Type | What it tests |
|--------|------|---------------|
| `cytokine_rate.toml` | 1D | Inflammation production rate vs healing outcome |
| `resolution_rate.toml` | 1D | M2 resolution speed vs inflammation duration |
| `collagen_vs_scar.toml` | 2D | Collagen deposition vs scar threshold (heatmap) |
| `diabetic_factors.toml` | 1D | M1 duration factor vs diabetic closure delay |

## Outcome measures

| Measure | Description |
|---------|-------------|
| `final` | Last value in the timeseries |
| `peak` | Maximum value |
| `auc` | Area under curve (trapezoidal, using time_h) |
| `time_to_N` | First time value exceeds N (e.g. `time_to_90` for 90% closure) |

## Directory structure

```
batch/
  batch.py           Multi-run consensus runner
  sweep.py           Parameter sweep runner
  analyze.py         Post-hoc analysis and plotting
  lib.py             Shared helpers (config, simulation, CSV, outcomes)
  configs/           Sweep config files
  results/           Output (gitignored)
  README.md
```
