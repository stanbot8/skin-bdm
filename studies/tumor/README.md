# Tumor Study

Neoplastic growth dynamics without wound interference. Seeds 200 tumor cells at t=0 and tracks proliferation, hypoxia suppression, and contact inhibition override.

## Configuration

- Duration: 30 days
- Wound: disabled
- Tumor: enabled (200 seed cells, apoptosis_rate=0.00047, steepness=2.0)
- Fibroblasts: disabled

## Usage

```bash
# Tumor only
python3 batch/batch.py -n 10 --study tumor --skin normal

# Cross-experiment comparison (tumor only, tumor+wound, aged+tumor)
python3 studies/tumor/study.py
```
