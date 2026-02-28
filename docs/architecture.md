> [Home](../README.md) / [Docs](README.md) / Architecture

# Architecture

Technical architecture of Skibidy, a hybrid agent-continuum skin tissue simulation. See also: [Guide](guide.md) | [Parameters](parameters.md) | [Treatments](treatments.md)

## UWYN (Use What You Need)

Skibidy uses a hybrid agent-continuum architecture:

- **Continuum at rest**: healthy skin is represented by diffusion fields (Calcium, KGF, O2, Water, Vascular, Inflammation, ImmunePressure, Stratum, Scar, Dermis, Elastin, Hyaluronan). No agents exist during homeostasis.
- **Event-driven agents**: when something happens (wound, infection, tumor), agents spawn to model the cellular response. Keratinocytes bootstrap from local field state at spawn time; immune cells arrive at configured post-wound delays and drive the inflammation field.
- **Per-cell handoff**: stable cornified cells dissolve back into the continuum field individually.
- **LOD toggles** per volume via `agents_enabled`: event modules escalate resolution at points of interest.

```
Corneum    [continuum]  barrier, desquamation
Granulosum [continuum]  keratohyalin, tight junctions
Spinosum   [continuum]  desmosomes; agents on event
Basale     [continuum]  agents on event (stem/TA cells, wound repair)
---------- basement membrane ----------
Dermis     [continuum]  vasculature, KGF + O2 source; agents on event
```

## CompositeField: composable PDE channels

Each diffusion field is a `PDE` subclass that owns its initialization, per-step source terms, and wound response. A `CompositeField` registers channels conditionally based on enabled modules and runs them through one scheduled operation:

| Channel | Init profile | Source term | Wound response |
|---------|-------------|-------------|----------------|
| Vascular | 1.0 in dermis, 0 in epidermis | Angiogenesis recovery (Stratum-gated) | Zero in wound cylinder (dermis) |
| Calcium | Sigmoid (basal-to-surface) | Stratum-gated recovery | Zero in cylinder |
| KGF | Exponential decay from dermis | None (static) | None (dermal signal) |
| O2 | Exponential decay from dermis | Dermal pinning (reads Vascular) | Zero in cylinder |
| Stratum | VolumeManager profile | None (agents write directly) | Zero above basement membrane |
| Water | Exponential decay from dermis | Dermal pinning (reads Vascular) + TEWL evaporation | Zero in cylinder |
| Inflammation | Zero (healthy) | Immune agents write directly (age-tapered) | None (starts at zero) |
| Scar | Empty | Collagen-driven (emergent) or inflammation integral | None (starts at zero) |
| TGF-beta | Zero | M2 macrophages + myofibroblasts write | None (starts at zero) |
| Collagen | Zero | Myofibroblasts deposit (proportional to local TGF-beta) | None (starts at zero) |
| MMP | Zero | M1 macrophages + fibroblasts produce | None (starts at zero) |
| Fibronectin | Zero | Fibroblasts deposit + serum leakage | None (starts at zero) |
| VEGF | Zero | Hypoxic tissue + M2 macrophages | None (starts at zero) |
| Biofilm | Zero | Logistic growth (when seeded) | None (starts at zero) |
| Dermis | Sub-layer profile (pap=1.0, ret=0.8, hyp=0.5) | Collagen-driven recovery, MMP degradation | Zero in wound dermis |
| ImmunePressure | Zero (healthy) | Immune cell cytokines (excludes wound DAMPs) | None (starts at zero) |
| Elastin | Dermal layer profile | Fibroblast production, MMP degradation | Wound disruption |
| Hyaluronan | Dermal layer profile | Fibroblast production | Wound disruption |
| Tumor | Zero | Binary marker from UWYN handoff (quiescent agent to field) | None (starts at zero) |

All PDE source terms share a `GridContext` helper that precomputes grid metadata (resolution, domain bounds, voxel coordinates) and wound geometry, eliminating repeated boilerplate across channels.

Cross-channel coupling happens inside `ApplySource()`: each PDE reads other channels via the `CompositeField` reference. For example, `OxygenPDE::ApplySource` reads the Vascular grid to determine local perfusion. New modules add a PDE subclass and register it; existing channels are untouched.

## Dynamic field coupling

Cross-channel source terms run each step via `CompositeField`:

- **Vascular perfusion** (`VascularPDE`): represents vessel density/integrity in the dermis. Healthy dermis = 1.0; wound destroys vessels to 0. After `perfusion_angio_delay` (~48h), angiogenesis gradually recovers perfusion, gated by tissue presence above (Stratum field as proxy for growth factor signals). Lateral capillary sprouting is handled by BioDynaMo's built-in diffusion. VascularPDE runs first in the CompositeField so O2 and Water read updated perfusion each step.
- **O2 source** (`OxygenPDE::ApplySource`): dermal vasculature maintains O2 below the basement membrane. Dermal voxels are pinned to `oxygen_basal_conc * local_perfusion`, so O2 supply tracks vessel recovery automatically.
- **Water source** (`WaterPDE::ApplySource`): dermal voxels are pinned to `water_basal_conc * local_perfusion`. Epidermal wound voxels lose water to TEWL evaporation (reduced by barrier integrity from Stratum) and gain water from dermal serum recovery.
- **Calcium recovery** (`CalciumPDE::ApplySource`): restores the calcium sigmoid in healed wound voxels, gated on local Stratum field value.
- **O2-modulated proliferation**: G1->S transition in `BasalDivision` scales with local O2, so margin cells (near intact vasculature) divide faster than wound-center cells.
- **Water-gated migration and proliferation**: tissue moisture gates both crawling speed (in `Migration`) and G1->S probability (in `BasalDivision`). After wounding, the wound bed must re-hydrate from dermal serum before effective re-epithelialization can begin.
- **Inflammation-gated migration and proliferation**: immune cell agents write pro-inflammatory cytokines into the Inflammation field with an age-dependent taper (neutrophils: `exp(-0.005 * age)`, M1 macrophages: `exp(-0.003 * state_age)`). High inflammation suppresses keratinocyte crawling and G1->S transition. Macrophages transition M1 to M2, actively consuming cytokines to resolve inflammation after ~3-5 days.

## Emergent behaviors

Several system-level behaviors emerge from local cell rules:

- **Inflammation curve**: the peak-and-decay shape emerges from neutrophil wave timing + individual cell age-dependent cytokine taper + M2 anti-inflammatory resolution, rather than being tuned via a single scalar rate parameter.
- **Scar formation**: `WriteScarValue` checks local collagen concentration when a keratinocyte dissolves. High collagen (above `scar_collagen_threshold`) produces scar tissue (stratum+5); low collagen produces normal stratum. This makes scar an emergent outcome of the fibroblast/collagen cascade.
- **Wound closure rate**: emerges from the interplay of O2 availability, water hydration, inflammation suppression, KGF growth factor signaling, and mechanical crowding.

## How to add a module

Follow these steps to add a new module. The fibroblast module is a good reference implementation.

### 1. Define field names

Add string constants to `src/core/field_names.h`:

```cpp
constexpr const char* kMyField = "MyField";
```

### 2. Create a PDE subclass

For simple agent-written fields (diffusion + decay, no profile), use `SimplePDE`:

```cpp
struct MyFieldPDE : public SimplePDE {
  explicit MyFieldPDE(const SimParam* sp)
      : SimplePDE(fields::kMyField, fields::kMyFieldId,
                  sp->my_diffusion, sp->my_decay) {}
};
```

For fields with custom initialization or source terms, inherit from `PDE` directly:

```cpp
struct MyFieldPDE : public PDE {
  const char* GetName() const override { return fields::kMyField; }
  int GetId() const override { return fields::kMyFieldId; }

  void Init(Simulation* sim) override { /* define + initialize grid */ }
  void ApplySource(Simulation* sim, const CompositeField& fields) override { /* per-step */ }
  void ApplyWound(Simulation* sim, real_t cx, real_t cy, real_t r) override { /* wound */ }
};
```

Override only the methods you need. `Init` runs once, `ApplySource` runs every step, `ApplyWound` runs once on wound event. Fields with `diffusion=0` and `decay=0` are automatically marked as prescribed (solver skipped).

### 3. Add parameters

Add parameters to `src/infra/sim_param.h`:

```cpp
real_t my_diffusion = 0.05;
real_t my_decay = 0.01;
```

### 4. Bind parameters to TOML

Add `BDM_ASSIGN_CONFIG_VALUE` entries in `src/skibidy.cc` with dotted section paths:

```cpp
BDM_ASSIGN_CONFIG_VALUE(my_diffusion, "skin.mymodule.diffusion");
BDM_ASSIGN_CONFIG_VALUE(my_decay, "skin.mymodule.decay");
```

### 5. Create module config file

Create `modules/mymodule/config.toml` with a `[skin.mymodule]` section:

```toml
[skin.mymodule]
enabled = false
diffusion = 0.05
decay = 0.01
```

The merge script (`scripts/config/merge_config.py`) discovers `modules/*/config.toml` automatically.

### 6. Register the PDE channel

In `src/core/registration.h`, register your PDE conditionally in `RegisterFields()`:

```cpp
if (sp->my_module_enabled) {
  cf.Add(new MyFieldPDE(sp));
}
```

### 7. Add agents and behaviors (if needed)

If the module needs agents, add them to the module directory:
- Create agent class in `modules/mymodule/my_agent.h` (inherit from `Cell`, use `BDM_AGENT_HEADER`)
- Create behavior in `modules/mymodule/my_behavior.h` (inherit from `Behavior`, use `BDM_BEHAVIOR_HEADER`)
- Create event in `modules/mymodule/my_event.h` (spawn agents at configured time/location)
- Register the event in `skibidy.h`

### 8. Add metrics columns

In `src/core/metrics.h`:
- Add counter variables in the agent loop
- Add column names to the CSV header
- Write values in the CSV row

### 9. Add tests

In `tests/test-suite.cc`:
- Unit test for the PDE (init profile, source terms, wound response)
- Unit test for agent/behavior logic
- Integration test (run N steps, verify expected outcome)
