# influx_si/MTF Format Specification

## Overview

influx_si uses Multiple TSV Files (MTF) format for 13C metabolic flux analysis.

| Property | Value |
|----------|-------|
| **File Format** | Tab-separated values (TSV) |
| **Documentation** | https://influx-si.readthedocs.io |
| **Legacy Format** | FTBL (still supported for computation) |
| **Introduced** | v6.0 |

## File Types

| Suffix | Required | Purpose |
|--------|----------|---------|
| `.netw` | Yes | Reaction network with atom transitions |
| `.linp` | Yes | Isotopic labeling inputs |
| `.miso` | Yes | MS/NMR isotopomer measurements |
| `.tvar` | No | Flux variable definitions |
| `.cnstr` | No | Linear constraints |
| `.mflux` | No | Measured fluxes |
| `.mmet` | No | Metabolite concentrations |
| `.opt` | No | Solver options |
| `.vmtf` | No | Variable MTF for batch experiments |

## File Schemas

### .netw (Network) - Required

Defines the reaction network with atom transition mappings.

**Format:**
```
reaction_name: Substrate (atoms) [arrow] Product (atoms)
```

**Arrow Types:**

| Arrow | Meaning |
|-------|---------|
| `->` | Irreversible (forward only) |
| `->>` | Non-negative net flux |
| `<->` | Reversible (bidirectional) |
| `<->>` | Reversible, non-negative net flux |

**Example:**
```
# Glycolysis
pgi: Glc6P (ABCDEF) <-> Fru6P (ABCDEF)
pfk: Fru6P (ABCDEF) -> FruBP (ABCDEF)
ald: FruBP (ABCDEF) <-> GA3P (CBA) + GA3P (DEF)

# Pentose Phosphate Pathway
zwf: Glc6P (ABCDEF) -> Gnt6P (ABCDEF)
gnd: Gnt6P (ABCDEF) -> Ru5P (BCDEF) + CO2 (A)

# TCA Cycle
pdh: Pyr (ABC) -> AcCoA (BC) + CO2 (A)
citsynth: AcCoA (AB) + OAA (abcd) -> ICit (dcbaBA)

# Scrambling reaction (symmetric metabolite)
fum_a: Mal (ABCD) <-> Fum (ABCD)
fum_b: Mal (ABCD) <-> Fum (DCBA)
```

**Conventions:**

| Convention | Description |
|------------|-------------|
| `#` | Comment line |
| `###` | Pathway header |
| UPPERCASE | Atom labels (A, B, C, ...) |
| `_ext` | External metabolite suffix |
| Stoichiometry | Precedes compound: `2.0*ATP ()` |

**Symmetric Metabolites:**
Instead of comma notation, influx_si uses separate reactions:
```
fum_a: Mal (ABCD) <-> Fum (ABCD)
fum_b: Mal (ABCD) <-> Fum (DCBA)
```
With constraint: `fum_a - fum_b == 0`

### .linp (Labeling Input) - Required

Specifies isotopic labeling of input substrates.

| Column | Description |
|--------|-------------|
| `Id` | Optional identifier |
| `Comment` | Optional notes |
| `Specie` | Substrate name (must match .netw) |
| `Isotopomer` | Binary pattern (`0`=unlabeled, `1`=labeled) |
| `Value` | Fraction (0-1) |

**Example:**
```
Id	Comment	Specie	Isotopomer	Value
		Gluc_ext	111111	0.25
		Gluc_ext	100000	0.69592
		Gluc_ext	000000	0.05408
		CO2_ext	0	1.0
```

**Rules:**
- Fractions for each species must sum to 1.0
- Missing isotopomers default to unlabeled
- Can specify natural abundance correction

### .miso (Isotopomer Measurements) - Required

MS and NMR isotopomer measurement data.

| Column | Description |
|--------|-------------|
| `Id` | Optional identifier |
| `Comment` | Optional notes |
| `Specie` | Metabolite name |
| `Fragment` | Atom positions (e.g., `1-3`, `1,2,3`, empty for whole) |
| `Dataset` | Measurement identifier (e.g., `MS-1`, `HSQC`) |
| `Isospecies` | Mass isotopomer or NMR pattern |
| `Value` | Measured value (empty/NA allowed) |
| `SD` | Standard deviation (required) |
| `Time` | Time point (empty for stationary) |

**Example (MS data):**
```
Id	Comment	Specie	Fragment	Dataset	Isospecies	Value	SD	Time
		Ala	2,3	MS-1	M0	0.35	0.01
		Ala	2,3	MS-1	M1	0.42	0.01
		Ala	2,3	MS-1	M2	0.23	0.01
		Glu	1-5	MS-2	M0	0.15	0.02
		Glu	1-5	MS-2	M1	0.28	0.02
		Glu	1-5	MS-2	M2	0.32	0.02
		Glu	1-5	MS-2	M3	0.18	0.02
		Glu	1-5	MS-2	M4	0.05	0.02
		Glu	1-5	MS-2	M5	0.02	0.02
```

**Example (Non-stationary):**
```
Id	Comment	Specie	Fragment	Dataset	Isospecies	Value	SD	Time
		AKG	1-5	inst	M0	0.95	0.01	0
		AKG	1-5	inst	M0	0.72	0.01	0.5
		AKG	1-5	inst	M0	0.45	0.01	1.0
		AKG	1-5	inst	M0	0.32	0.01	2.0
```

**NMR Notation:**
- `2->` : Singlet at position 2
- `2->1` : Doublet (coupled to position 1)
- `2->1,3` : Doublet of doublets
- Cumomer pattern: `0`, `1`, `x`

### .tvar (Variable Types)

Defines flux parameter classes and starting values.

| Column | Description |
|--------|-------------|
| `Name` | Variable identifier (reaction name) |
| `Kind` | `NET`, `XCH`, or `METAB` |
| `Type` | `F`=free, `D`=dependent, `C`=constrained |
| `Value` | Starting value (required for F and C) |

**Example:**
```
Id	Comment	Name	Kind	Type	Value
		v_upt	NET	F	1.0
		pgi	NET	D
		pgi	XCH	F	0.75
		pfk	NET	C	0.8
		Pyr	METAB	F	0.5
```

**Type Definitions:**

| Type | Description |
|------|-------------|
| `F` | Free parameter (to be estimated) |
| `D` | Dependent (calculated from constraints) |
| `C` | Constrained (fixed value) |

### .cnstr (Constraints)

Linear equality and inequality constraints.

| Column | Description |
|--------|-------------|
| `Kind` | `NET`, `XCH`, or `MET` |
| `Formula` | Linear expression |
| `Operator` | `==`, `>=`, `<=` |
| `Value` | Right-hand side (float or expression) |

**Example:**
```
Id	Comment	Kind	Formula	Operator	Value
		NET	v_upt	==	1
		NET	pgi+pfk	==	1
		NET	fum_a-fum_b	==	0
		NET	edd	>=	0.0001
		XCH	pgi	<=	0.99
```

**Formula Syntax:**
- Variables: reaction names
- Operators: `+`, `-`, `*`
- Coefficients precede variables: `0.5*v1`

### .mflux (Measured Fluxes)

Experimentally measured flux values.

| Column | Description |
|--------|-------------|
| `Flux` | Flux identifier |
| `Value` | Measured value |
| `SD` | Standard deviation |

**Example:**
```
Id	Comment	Flux	Value	SD
		v_upt	1.02	0.05
		out_Ac	0.213	0.02
```

### .mmet (Metabolite Concentrations)

Metabolite pool size measurements.

| Column | Description |
|--------|-------------|
| `Specie` | Metabolite identifier |
| `Value` | Concentration |
| `SD` | Standard deviation |

**Example:**
```
Id	Comment	Specie	Value	SD
		Pyr	0.426	0.05
		PEP	0.182	0.02
		AKG	0.41e-5	0.01e-5
```

### .opt (Options)

Solver configuration and execution parameters.

| Column | Description |
|--------|-------------|
| `Name` | Option name |
| `Value` | Option value |

**Example:**
```
Id	Comment	Name	Value
		commandArgs	--TIMEIT --noscale
		posttreat_R	plot_results.R
		optctrl.nlsic.maxiter	50
```

**Common Options:**

| Option | Description |
|--------|-------------|
| `--TIMEIT` | Enable timing |
| `--noscale` | Disable scaling |
| `--ln` | Log normalization |
| `posttreat_R` | Post-processing R script |

## File Naming Convention

All files share a base name:
```
mymodel.netw
mymodel.linp
mymodel.miso
mymodel.tvar
mymodel.cnstr
mymodel.mflux
mymodel.mmet
mymodel.opt
```

Non-stationary experiments may have separate measurement files:
```
mymodel_inst.miso   # Instantaneous measurements
mymodel_stat.miso   # Stationary measurements
```

## Output Files

influx_si generates output in MTF format:

| File | Description |
|------|-------------|
| `*.miso.sim` | Simulated isotope measurements |
| `*.mflux.sim` | Simulated fluxes |
| `*.mmet.sim` | Simulated concentrations |
| `*.tvar.sim` | Estimated parameters with SD |
| `*.stat` | Chi-square test results |
| `*.log` | Runtime log |

## Data Files in Repository

| Path | Description |
|------|-------------|
| `data/convert/influx_si/` | EC model conversion files |
| `data/tests/mtf_models/e_coli/` | E. coli model |
| `data/tests/mtf_models/Nouveau dossier/` | Acetate model |

## Example: E. coli Model

```
e_coli/
├── e_coli.netw      # 88 reactions
├── e_coli.linp      # Glucose labeling (U-13C + 1-13C mix)
├── e_coli.miso      # 9 MS datasets
├── e_coli.tvar      # 99 flux variables
├── e_coli.cnstr     # 14 constraints
├── e_coli.mflux     # 1 measured flux
└── e_coli.opt       # Solver options
```

## Relationship to FTBL

MTF replaced FTBL as the primary format in influx_si v6.0:

| Aspect | FTBL | MTF |
|--------|------|-----|
| Format | Single file, sections | Multiple TSV files |
| Editing | Complex | Easy (spreadsheet) |
| Version control | Difficult | Easy (diff-friendly) |
| Programmatic | Custom parser | Standard TSV |

FTBL is still used internally for computation.
