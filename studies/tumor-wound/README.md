# Tumor-Wound Study

Coupled tumor growth and wound healing with sequential activation. Tumor seeds at t=0, wound triggers at 200h (~8 days). Tests interaction between neoplastic progression and repair processes.

## Configuration

- Duration: 75 days
- Wound: enabled (trigger at 200h)
- Tumor: enabled (200 seed cells at t=0)
- Fibroblasts: enabled

## Usage

```bash
python3 batch/batch.py -n 10 --study tumor-wound --skin normal
```
