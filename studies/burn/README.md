# Burn Wound Study

Thermal injury simulation using Jackson's burn wound model (1953). Burns are the fourth most common type of trauma worldwide, killing approximately 180,000 people per year and causing significant morbidity from scarring and contracture.

Jackson's model describes three concentric zones of injury:

- **Zone of coagulation** (center): irreversible necrosis from direct thermal damage, protein denaturation
- **Zone of stasis** (middle): ischemic and potentially salvageable tissue; converts to necrosis without treatment over 24 to 48 hours
- **Zone of hyperemia** (outer): inflammatory response, typically recovers spontaneously

Burn depth classification:

| Classification | depth_fraction | Layers affected |
|---|---|---|
| Superficial | < 0.3 | Epidermis only |
| Superficial partial | 0.3 to 0.5 | Upper dermis |
| Deep partial | 0.5 to 0.7 | Deep dermis |
| Full thickness | >= 0.7 | Through dermis |

## Running

```bash
# Baseline burn (untreated partial thickness, 10-run consensus)
python3 batch/batch.py -n 10 --skin burn --study burn

# With cooling treatment
python3 batch/batch.py -n 10 --skin burn --study burn --treatment cooling

# With silver sulfadiazine
python3 batch/batch.py -n 10 --skin burn --study burn --treatment silver_sulfadiazine

# Run experiments
python3 studies/run_experiments.py studies/burn/experiments/depth_spectrum.toml
python3 studies/run_experiments.py studies/burn/experiments/stasis_rescue.toml
python3 studies/run_experiments.py studies/burn/experiments/treatment_comparison.toml
python3 studies/run_experiments.py studies/burn/experiments/scar_prevention.toml
```

## Configuration

- Duration: 60 days
- Wound: enabled (t=0)
- Burn module: enabled
- Scar module: enabled
- Profile: `profiles/burn.toml` (hypertrophic scar risk, massive immune response, barrier loss)

## Treatments

| Treatment | File | Mechanism |
|---|---|---|
| Cooling | [cooling.toml](treatments/cooling.toml) | Immediate cold water first aid, preserves stasis zone |
| Debridement | [debridement.toml](treatments/debridement.toml) | Surgical eschar removal, exposes viable wound bed |
| Silver sulfadiazine | [silver_sulfadiazine.toml](treatments/silver_sulfadiazine.toml) | Topical antimicrobial, standard burn care |
| Skin substitute | [skin_substitute.toml](treatments/skin_substitute.toml) | Bioengineered scaffold (Integra/Biobrane), reduces water loss |
| Pressure garment | [pressure_garment.toml](treatments/pressure_garment.toml) | Compression therapy, prevents hypertrophic scar |

## Experiments

| Experiment | File | Configs | Runs |
|---|---|---|---|
| Depth spectrum | [depth_spectrum.toml](experiments/depth_spectrum.toml) | 4 (superficial to full thickness) | 5/config |
| Stasis rescue | [stasis_rescue.toml](experiments/stasis_rescue.toml) | 4 (untreated vs cooling at 0h, 6h, 24h) | 5/config |
| Treatment comparison | [treatment_comparison.toml](experiments/treatment_comparison.toml) | 6 (untreated through full protocol) | 5/config |
| Scar prevention | [scar_prevention.toml](experiments/scar_prevention.toml) | 3 (no prevention, pressure, pressure + silicone) | 5/config |
