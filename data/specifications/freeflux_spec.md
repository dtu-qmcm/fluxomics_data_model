# FreeFlux Format Specification

## Overview

FreeFlux uses tabular files (TSV/CSV/XLSX) for 13C metabolic flux analysis models.

| Property | Value |
|----------|-------|
| **File Formats** | TSV, CSV, XLSX |
| **Python Package** | `freeflux` |
| **Documentation** | https://freeflux.readthedocs.io |
| **GitHub** | https://github.com/Chaowu88/freeflux |

## File Types

| File | Required | Purpose |
|------|----------|---------|
| `reactions` | Yes | Metabolic network with atom mappings |
| `fluxes` | No | Simulated/reference flux values |
| `concentrations` | No | Metabolite pool sizes |
| `measured_MDVs` | No | Steady-state labeling measurements |
| `measured_fluxes` | No | Flux measurements with uncertainty |
| `measured_inst_MDVs` | No | Time-course labeling (non-stationary) |

Files can have `.tsv`, `.csv`, or `.xlsx` extensions.

## File Schemas

### reactions.tsv (Required)

Defines the metabolic reaction network with atom transition mappings.

| Column | Description |
|--------|-------------|
| `#reaction_ID` | Unique reaction identifier |
| `reactant_IDs(atom)` or `substrate_IDs(atom)` | Reactants with atom mapping |
| `product_IDs(atom)` | Products with atom mapping |
| `reversibility` | `0` = irreversible, `1` = reversible |

**Example:**
```
#reaction_ID	reactant_IDs(atom)	product_IDs(atom)	reversibility
v1	OAA(abcd)+AcCoA(ef)	Cit(dcbfea)	0
v2	Cit(abcdef)	AKG(abcde)+CO2(f)	0
v3	AKG(abcde)	SucCoA(bcde)+CO2(a)	0
v4	SucCoA(abcd)	Suc(abcd,dcba)	0
v5	Suc(abcd,dcba)	Fum(abcd,dcba)	0
v6	Fum(abcd,dcba)	OAA(abcd)	1
v7	AKG(abcde)+CO2(f)	Glu(abcdef)	0
```

#### Atom Mapping Conventions

- **Lowercase letters**: `a, b, c, d, e, f, ...`
- **Format**: `Metabolite(atom_mapping)`
- **Multiple reactants/products**: Connected with `+`
- **Symmetric molecules**: Comma-separated variants

**Atom Mapping Examples:**

| Pattern | Meaning |
|---------|---------|
| `G6P(abcdef)` | 6-carbon metabolite, atoms a-f |
| `Pyr(abc)` | 3-carbon metabolite |
| `OAA(abcd)+AcCoA(ef)` | Two reactants combined |
| `Cit(dcbfea)` | Product with rearranged atoms |
| `Suc(abcd,dcba)` | Symmetric 4-carbon (two equivalent forms) |
| `CO2(a)` | Single carbon released |

**Symmetric Metabolites:**

Molecules like succinate and fumarate have symmetric carbon skeletons:
```
Suc(abcd,dcba)    # Succinate: C1-C2-C3-C4 equivalent to C4-C3-C2-C1
Fum(abcd,dcba)    # Fumarate: same symmetry
```

Each variant has equal probability (0.5 for 2 variants).

### fluxes.tsv

Simulated or reference metabolic flux values.

| Column | Description |
|--------|-------------|
| `#flux_ID` | Flux identifier (with optional `_f`/`_b` suffix) |
| `value` | Numeric flux value |

**Example:**
```
#flux_ID	value
v1	10
v2	10
v3	7
v4	7
v5	7
v6_f	12.5
v6_b	7.5
v7	3
```

**Flux Naming Conventions:**

| Pattern | Meaning |
|---------|---------|
| `v1` | Net flux (irreversible reaction) |
| `v6_f` | Forward flux (reversible reaction) |
| `v6_b` | Backward flux (reversible reaction) |

**Flux Relationships:**
- Net flux = forward - backward
- For `v6`: net = 12.5 - 7.5 = 5.0

### concentrations.tsv

Metabolite pool sizes for kinetic simulations.

| Column | Description |
|--------|-------------|
| `#metab_ID` | Metabolite identifier |
| `value` | Concentration value |

**Example:**
```
#metab_ID	value
OAA	0.1
Cit	5
AKG	0.3
SucCoA	0.05
Suc	1
Fum	0.2
Glu	0.5
```

**Units:** Typically mM or mol/gCDW (cell dry weight).

### measured_MDVs.tsv

Steady-state mass distribution vector (MDV) measurements from MS.

| Column | Description |
|--------|-------------|
| `#fragment_ID` | Format: `Metabolite_Positions` |
| `mean` | Comma-separated MDV values (M+0, M+1, ...) |
| `sd` | Comma-separated standard deviations |

**Example:**
```
#fragment_ID	mean	sd
Glu_12345	0.328,0.276,0.274,0.088,0.03,0.004	0.01,0.01,0.01,0.01,0.01,0.01
Ala_23	0.42,0.35,0.23	0.01,0.01,0.01
```

**Fragment ID Format:**
- `Glu_12345` - Glutamate, measuring carbons 1,2,3,4,5
- `Ala_23` - Alanine, measuring carbons 2,3 only

**MDV Interpretation:**
- M+0 (0.328): 32.8% unlabeled
- M+1 (0.276): 27.6% with 1 labeled carbon
- M+2 (0.274): 27.4% with 2 labeled carbons
- Sum should approximately equal 1.0

### measured_fluxes.tsv

Experimentally measured flux values with uncertainty.

| Column | Description |
|--------|-------------|
| `#reaction_ID` | Reaction identifier |
| `mean` | Measured flux value |
| `sd` | Standard deviation |

**Example:**
```
#reaction_ID	mean	sd
v1	10	1
v7	3	0.5
```

### measured_inst_MDVs.tsv (Non-stationary)

Time-course MDV measurements for isotopically non-stationary experiments.

| Column | Description |
|--------|-------------|
| `#fragment_ID` | Fragment identifier |
| `time` | Time point (in experiment time units) |
| `mean` | Comma-separated MDV values |
| `sd` | Comma-separated standard deviations |

**Example:**
```
#fragment_ID	time	mean	sd
Glu_12345	0	0.948,0.051,0.001,0.0,0.0,0.0	0.01,0.01,0.01,0.01,0.01,0.01
Glu_12345	0.1	0.928,0.06,0.012,0.0,0.0,0.0	0.01,0.01,0.01,0.01,0.01,0.01
Glu_12345	0.2	0.873,0.085,0.041,0.001,0.0,0.0	0.01,0.01,0.01,0.01,0.01,0.01
Glu_12345	0.5	0.682,0.158,0.128,0.024,0.007,0.001	0.01,0.01,0.01,0.01,0.01,0.01
Glu_12345	1	0.487,0.233,0.212,0.048,0.017,0.003	0.01,0.01,0.01,0.01,0.01,0.01
Glu_12345	2	0.381,0.274,0.262,0.063,0.018,0.001	0.01,0.01,0.01,0.01,0.01,0.01
```

**Time Points:**
- t=0: Nearly unlabeled (start of labeling)
- t>0: Increasing labeling as 13C propagates through metabolism
- Steady-state: Labeling stabilizes

## File Format Conventions

### TSV Structure
- Tab-separated values
- First row contains column headers (may start with `#`)
- Comments may start with `#`
- Empty cells and `NaN` values are handled gracefully

### Column Name Variations

The parser should accept these variations:

| Canonical | Also Accepted |
|-----------|---------------|
| `reaction_ID` | `#reaction_ID`, `Reaction_ID`, `rxn_id` |
| `reactant_IDs(atom)` | `substrate_IDs(atom)`, `substrates` |
| `product_IDs(atom)` | `products` |
| `reversibility` | `reversible`, `rev`, `bidirectional` |
| `mean` | `value`, `avg`, `average` |
| `sd` | `std`, `stdev`, `std_dev`, `error` |

### Edge Cases and Error Handling

| Situation | Default Handling | Notes |
|-----------|------------------|-------|
| Missing `reversibility` column | Default to `1` (reversible) | Emit warning |
| Empty cell in `reversibility` | Default to `1` | Emit warning |
| MDV sum < 0.5 | Error | Data likely corrupted |
| MDV sum ≠ 1.0 (but > 0.5) | Normalize and warn | Divide by sum |
| `NaN` in measurement value | Skip row | Emit warning |
| `inf` in any numeric field | Error | Invalid data |
| Negative MDV value | Error | Physically impossible |
| Duplicate reaction ID | Error | IDs must be unique |
| Duplicate fragment ID | Merge with warning | Common in replicate data |
| Missing required column | Error | List missing columns |
| Extra columns | Ignore | Emit debug message |

## Data Relationships

```
Model Structure:
├── reactions.tsv          → Defines network topology
│   ├── Metabolites        (extracted from reactions)
│   └── Atom mappings      (carbon tracking)
├── fluxes.tsv             → Reference/simulated values
├── concentrations.tsv     → Pool sizes (non-stationary)
└── Measurements
    ├── measured_MDVs.tsv      → Stationary experiment
    ├── measured_inst_MDVs.tsv → Non-stationary experiment
    └── measured_fluxes.tsv    → Flux constraints
```

## Data Files in Repository

| Path | Description |
|------|-------------|
| `data/convert/FreeFlux/` | Conversion test files |
| `data/tests/freeflux_models/toy/` | Simple TCA model (TSV) |
| `data/tests/freeflux_models/synechocystis/` | Cyanobacterium (XLSX) |
| `data/tests/freeflux_models/ecoli/` | E. coli model (XLSX) |

## Example: Toy Model

A minimal TCA cycle model demonstrating all file types:

**reactions.tsv** (7 reactions):
```
v1: OAA + AcCoA → Cit (citrate synthase)
v2: Cit → AKG + CO2 (isocitrate dehydrogenase)
v3: AKG → SucCoA + CO2 (alpha-ketoglutarate dehydrogenase)
v4: SucCoA → Suc (succinyl-CoA synthetase)
v5: Suc → Fum (succinate dehydrogenase)
v6: Fum ↔ OAA (fumarase + malate dehydrogenase)
v7: AKG + CO2 → Glu (glutamate dehydrogenase)
```

**Experimental setup:**
- Feed 13C-labeled substrate
- Measure Glu MDV at steady-state or over time
- Estimate flux distribution

## Python API Example

```python
from freeflux import Model

# Load model from files
model = Model('toy')
model.read_from_file(
    reactions_file='reactions.tsv',
    measured_MDVs_file='measured_MDVs.tsv',
    measured_fluxes_file='measured_fluxes.tsv'
)

# Set labeling strategy
model.set_labeling_strategy(
    substrates=['AcCoA', 'OAA'],
    labeling_patterns=['11', '0000'],
    percentages=[1.0, 1.0]
)

# Simulate or estimate fluxes
simulator = model.get_simulator('stationary')
results = simulator.simulate(fluxes)
```

## Code-Based Constraints (Not File-Based)

FreeFlux does NOT support file-based constraints. All bounds and constraints must be set via Python API:

### Flux Bounds

```python
# Set bounds for all fluxes
fit.set_flux_bounds('all', bounds=[-100, 100])

# Set bounds for specific flux
fit.set_flux_bounds('v1', bounds=[0, 50])
```

### Measured Fluxes (Soft Equality Constraint)

```python
# Measured flux acts as fitting target with uncertainty
fit.set_measured_flux('v1', mean=10, sd=1)
```

### Concentration Bounds (INST Only)

```python
# NOTE: Only defines initial guess sampling range, NOT optimization constraint!
fit.set_concentration_bounds('all', bounds=[0.001, 10])
```

**IMPORTANT**: From FreeFlux docs: "Bounds specified by `set_concentration_bounds` are not used as constraints in the optimization, it just defines the sampling range of the initial guess of concentrations."

### Limitations

| Feature | Supported |
|---------|-----------|
| Flux bounds | ✓ (code only) |
| Measured flux constraint | ✓ (code only) |
| Linear equality (`v1 - v2 = 0`) | ✗ |
| Linear inequality | ✗ |
| File-based constraints | ✗ |
| Pool size optimization bounds | ✗ (only initial guess sampling) |
