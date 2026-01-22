# INCA Format Specification

## Overview

INCA (Isotopomer Network Compartmental Analysis) uses MATLAB script files for 13C-MFA models.

| Property | Value |
|----------|-------|
| **File Format** | MATLAB script (`.m`) |
| **Software** | INCA (commercial) |
| **Language** | MATLAB |

## Document Structure

INCA models are defined using MATLAB functions and object assignments:

```matlab
% ========================================
% Reactions
% ========================================
r = reaction({
  'Substrate (atoms) <-> Product (atoms)';  % v1
  ...
});

% ========================================
% Tracers (Isotope Labeling)
% ========================================
t = tracer({
  'tracer_name: metabolite @ labeled_positions';
  ...
});
t.frac = [fraction1, fraction2, ...];

% ========================================
% Flux Measurements
% ========================================
f = data('reaction_IDs');
f.val = [values];
f.std = [standard_deviations];

% ========================================
% MS Measurements
% ========================================
d = msdata({
  'measurement_name: metabolite @ atom_positions';
  ...
});
d.idvs = idv;  % Isotopomer Distribution Vectors

% ========================================
% Experiment Setup
% ========================================
x = experiment(t);
x.data_flx = f;
x.data_ms = d;

% ========================================
% Model Assembly
% ========================================
m = model(r);
m.expts = x;

% ========================================
% Symmetric Metabolites
% ========================================
m.mets{'Suc'}.sym = list('rotate180', atommap('1:4 2:3 3:2 4:1'));

% ========================================
% Initial Values and Bounds
% ========================================
m.rates.flx.val = [flux_values];
m.rates.flx.lb = lower_bound;
m.rates.flx.ub = upper_bound;

m.states.val = [concentration_values];
m.states.lb = lower_bound;
m.states.ub = upper_bound;

% ========================================
% Simulation Options
% ========================================
m.options.sim_ss = false;
m.options.int_tspan = [time_points];
```

## Elements

### reaction()

Defines the metabolic reaction network.

**Syntax:**
```matlab
r = reaction({
  'Reactant (atoms) <-> Product (atoms)';
  'Reactant (atoms) -> Product (atoms)';
  ...
});
```

**Conventions:**

| Symbol | Meaning |
|--------|---------|
| `<->` | Reversible reaction |
| `->` | Irreversible reaction |
| `(abc)` | Atom mapping (lowercase letters) |
| `+` | Multiple reactants/products |
| `0.5 Met` | Stoichiometric coefficient |

**Examples:**
```matlab
r = reaction({
  'G6P (abcdef) <-> F6P (abcdef)';               % v1 - isomerase
  'FBP (abcdef) <-> DHAP (cba) + GAP (def)';     % v3 - aldolase
  'Pyr (abc) -> AcCoA (bc) + CO2 (a)';           % v19 - PDH
  'OAA (abcd) + AcCoA (ef) -> Cit (dcbfea)';     % v20 - citrate synthase
  '0.488 Ala (abc) + 3.0 ATP () -> biomass';     % Biomass with stoich
});
```

**Atom Mapping:**
- Uses lowercase letters: `a, b, c, d, e, f, ...`
- Tracks carbon positions through reactions
- Order in product shows where each reactant atom ends up

### tracer()

Specifies isotope labeling strategy.

**Syntax:**
```matlab
t = tracer({
  'tracer_name: metabolite @ labeled_positions';
  ...
});
t.frac = [fractions];
```

**Examples:**
```matlab
% Single tracer
t = tracer({
  'U-13C_Gluc: Gluc.ext @ 1 2 3 4 5 6';  % Uniformly labeled
});
t.frac = [1.0];

% Multiple tracers (mixture)
t = tracer({
  '1-13C_Gluc: Gluc.ext @ 1';            % Position 1 labeled
  'U-13C_Gluc: Gluc.ext @ 1 2 3 4 5 6';  % All positions labeled
});
t.frac = [0.75, 0.25];  % 75% C1-labeled, 25% U-labeled

% CO2 labeling (photosynthesis)
t = tracer({
  '13C_CO2: CO2.ex @ 1';
  'CO2_unlabeled: CO2.ex @ ';  % Empty = unlabeled
});
t.frac = [0.99, 0.01];
```

**Position Notation:**
- `@ 1 2 3` - Positions 1, 2, 3 are labeled
- `@ 1` - Only position 1 labeled
- `@ ` (empty) - Unlabeled

### data() (Flux Measurements)

Measured flux values for fitting.

**Syntax:**
```matlab
f = data('reaction_ID1 reaction_ID2 ...');
f.val = [measured_values];
f.std = [standard_deviations];
```

**Example:**
```matlab
f = data('R61 R68 R64 R66 R67');
f.val = [92.5864, 1.983, 129.17, 177.54, 0.904];
f.std = f.val / 20;  % 5% relative error

% Add measurement noise (for simulation)
f.val = normrnd(f.val, f.std);
```

### msdata() (MS Measurements)

Mass spectrometry isotopologue measurements.

**Syntax:**
```matlab
d = msdata({
  'measurement_name: metabolite @ atom_positions';
  ...
});
d.idvs = idv;  % IDV data structure
```

**Examples:**
```matlab
d = msdata({
  'Ala2: Ala @ 2 3';           % 2-carbon fragment
  'Ala3: Ala @ 1 2 3';         % 3-carbon fragment
  'AKG5: AKG @ 1 2 3 4 5';     % Full molecule
  'Glu5: Glu @ 1 2 3 4 5';
  'Asp4: Asp @ 1 2 3 4';
  'Phe9: Phe @ 1 2 3 4 5 6 7 8 9';  % 9-carbon amino acid
});

% Technical replicates
d = msdata({
  'Asp2a: Asp @ 2 3';  % Replicate a
  'Asp2b: Asp @ 2 3';  % Replicate b
  'Asp2c: Asp @ 2 3';  % Replicate c
});
```

**Measurement Naming Convention:**
- `Ala2` - Alanine, 2-carbon fragment
- `Asp2a` - Aspartate 2-carbon, replicate 'a'

### experiment()

Combines tracers and measurements into an experiment.

```matlab
x = experiment(t);
x.data_flx = f;   % Flux measurements
x.data_ms = d;    % MS measurements
```

### model()

Assembles the complete model.

```matlab
m = model(r);
m.expts = x;      % Add experiment(s)
```

### Symmetric Metabolites

Define symmetry for molecules like succinate and fumarate.

**Syntax:**
```matlab
m.mets{'metabolite_name'}.sym = list('symmetry_name', atommap('mapping'));
```

**Examples:**
```matlab
% Succinate: 4-carbon symmetric (C1-C2-C3-C4 = C4-C3-C2-C1)
m.mets{'Suc'}.sym = list('rotate180', atommap('1:4 2:3 3:2 4:1'));

% Fumarate: same symmetry
m.mets{'Fum'}.sym = list('rotate180', atommap('1:4 2:3 3:2 4:1'));
```

**Mapping Format:** `original_pos:new_pos` pairs

### Flux Values and Bounds

Initial values and bounds for flux estimation.

```matlab
% Flux bounds (all reactions)
m.rates.flx.lb = 1e-7;    % Lower bound
m.rates.flx.ub = 1e7;     % Upper bound

% Flux initial values (ordered by reaction)
m.rates.flx.val = [
  10.5, ...   % v1
  8.2, ...    % v2
  ...
];
```

### State Variables (Pool Sizes)

Metabolite concentrations for non-stationary analysis.

```matlab
% Pool size bounds
m.states.lb = 1e-6;
m.states.ub = 1;

% Initial concentrations (ordered by metabolite)
m.states.val = [
  0.5, ...    % G6P
  0.3, ...    % F6P
  1.6e-3, ... % PEP
  ...
];
```

### Simulation Options

```matlab
m.options.int_tspan = [0 0.1 0.2 0.5 1 2];  % Time points (hours)
m.options.sim_tunit = 'h';                   % Time unit
m.options.fit_reinit = true;                 % Reinitialize
m.options.sim_ss = false;                    % Non-stationary
m.options.sim_sens = false;                  % Sensitivity analysis
m.options.int_reltol = 1e-10;               % Integration tolerance
```

## Naming Conventions

| Pattern | Meaning |
|---------|---------|
| `v1`, `v2` | Numbered reactions |
| `R1`, `R39` | Alternative reaction IDs |
| `.ext`, `.ex` | External/boundary metabolites |
| `Gluc.ext` | External glucose |
| `CO2.ext` | External CO2 |
| `Ala2`, `Glu5` | MS fragment (metabolite + # carbons) |
| `Ala2a` | Technical replicate 'a' |
| `biomass` | Biomass pseudo-metabolite |
| `Dummy` | Balancing species |

## Metabolic Pathways

Typical organization in INCA files:

```matlab
% Glycolysis (v1-v10)
'G6P (abcdef) <-> F6P (abcdef)';           % v1

% Pentose Phosphate Pathway (v11-v20)
'G6P (abcdef) -> Gnt6P (abcdef)';          % v11

% TCA Cycle (v21-v30)
'Pyr (abc) -> AcCoA (bc) + CO2 (a)';       % v21

% Amino Acid Biosynthesis (v31-v50)
'Pyr (abc) -> Ala (abc)';                  % v31

% Transport (v51-v60)
'Gluc.ext (abcdef) -> G6P (abcdef)';       % v51

% Biomass
'0.488 Ala + 0.32 Glu + ... -> biomass';   % v61
```

## Data Files in Repository

| Path | Description |
|------|-------------|
| `data/convert/INCA/EC.m` | E. coli model |
| `data/convert/INCA/EC_poolsizes.m` | E. coli with pool sizes |
| `data/convert/INCA/Syn.m` | Cyanobacterium model |

## Example: E. coli Model Summary

```matlab
% EC.m structure:
% - 68 reactions (glycolysis, PPP, TCA, amino acids, transport)
% - 2 tracers (1-13C glucose + U-13C glucose mix)
% - 7 measured fluxes
% - 30+ MS measurements
% - Symmetric metabolites: Suc, Fum
% - 66 metabolite pools
```

## Conversion Considerations

When converting INCA to other formats:

1. **Atom Mapping**: INCA uses lowercase letters (same as FreeFlux)
2. **Flux Representation**: Extract from `m.rates.flx.val` array
3. **Pool Sizes**: From `m.states.val` array
4. **Symmetry**: Must be translated to format-specific notation
5. **Reaction IDs**: Generate from comments or use `v1`, `v2`, etc.
6. **External Metabolites**: `.ext` suffix indicates boundary species

## Array Index Conventions

INCA uses MATLAB's 1-indexed arrays. The mapping works as follows:

| Array | Index Meaning |
|-------|---------------|
| `m.rates.flx.val(i)` | Flux value for i-th reaction (1-indexed) |
| `m.states.val(i)` | Concentration for i-th metabolite (1-indexed) |
| `t.frac(i)` | Fraction for i-th tracer isotopomer (1-indexed) |

### Reaction Order

Reactions are numbered in the order they appear in the `reaction({...})` definition:

```matlab
r = reaction({
  'A -> B';    % i=1
  'B -> C';    % i=2
  'C -> D';    % i=3
});

m.rates.flx.val(1)  % Flux for 'A -> B'
m.rates.flx.val(2)  % Flux for 'B -> C'
```

### Metabolite Order

Metabolites are ordered alphabetically by default, but may be reordered based on internal INCA processing. To get reliable ordering:

```matlab
% Get metabolite names in order
met_names = m.mets.id;  % Cell array of names

% Map index to name
met_name = m.mets{i}.id;
```

### Parsing Guidelines

When parsing INCA files:

1. Extract reaction order from `reaction({...})` block
2. Build index-to-ID mapping using reaction order
3. Use mapping when reading flux arrays
4. Handle comments to extract meaningful reaction IDs (e.g., `% v1 - isomerase`)
