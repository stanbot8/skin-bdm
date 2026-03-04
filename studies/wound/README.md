# Wound Healing Study

Punch biopsy wound healing with full repair cascade: hemostasis, inflammation, proliferation, and remodeling. Primary validation target against published wound closure kinetics.

## Configuration

- Duration: 30 days
- Wound: enabled (t=0)
- pH, fibroblasts, hemostasis: enabled

## Usage

```bash
# Single skin phenotype (10-run consensus with literature validation)
python3 batch/batch.py -n 10 --study wound --skin normal --validate

# Cross-phenotype comparison (normal, aged, diabetic)
python3 studies/wound/skin_comparison.py
python3 studies/wound/skin_comparison.py --quick  # 5 runs per phenotype
```
