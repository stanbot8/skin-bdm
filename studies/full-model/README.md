# Full Model

Wound healing with all mechanistic toggles enabled. Validates that coupling feedback loops (collagen deposition, M1/M2 transition, VEGF production) improve biological accuracy over simplified models.

## Configuration

- Duration: 30 days
- Wound: enabled (t=0)
- pH: enabled
- Fibroblasts: enabled (mechanistic collagen deposition)
- Hemostasis: enabled
- Immune: mechanistic M1/M2 transition
- Angiogenesis: mechanistic VEGF production

## Usage

```bash
python3 batch/batch.py -n 10 --study full-model --skin normal --validate
```
