# Fluxomics Data Model

[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![Tests](https://github.com/dtu-qmcm/fluxomics_data_model/actions/workflows/run_tests.yml/badge.svg)](https://github.com/dtu-qmcm/fluxomics_data_model/actions/workflows/run_tests.yml)
[![Supported Python versions: 3.12 and newer](https://img.shields.io/badge/python->=3.12-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A universal Python library for converting and representing fluxomics data across multiple formats. This tool enables seamless interoperability between different 13C metabolic flux analysis (13C-MFA) software tools.

## Overview

Fluxomics Data Model provides a unified data structure for representing metabolic network models, atom mappings, experimental configurations, and measurement data. It serves as a bridge between different fluxomics file formats, allowing researchers to:

- **Parse** models from FluxML (13CFlux), SBML, CSV/TSV sheets, and MTF (influx_si) formats
- **Convert** between formats without data loss
- **Export** to formats compatible with major 13C-MFA tools:
  - [INCA](https://mfa.vueinnovations.com/) (via [incawrapper](https://github.com/biosustain/incawrapper))
  - [13CFlux/13CFlux2](https://www.13cflux.net/)
  - [influx_si](https://metasys.insa-toulouse.fr/software/influx/)

## Supported Formats

| Format | Software | Read | Write | Status |
|--------|----------|------|-------|--------|
| FluxML | 13CFlux2, 13CFlux3 | ✅ | ✅ | Fully implemented |
| MTF | influx_si | ✅ | ✅ | Implemented (cumomer, `.mmet`, METAB vars supported) |
| CSV/TSV/XLSX | Freeflux | ✅ | 🔲 | Read implemented (no tracers/constraints from files) |
| MATLAB | INCA | 🔲 | 🔲 | Planned |
| XML | SBML | 🔲 | 🔲 | Planned (dependency ready) |

## Installation

### Using pip

```bash
pip install fluxomics-data-model
```

### Using uv (recommended for development)

```bash
# Clone the repository
git clone https://github.com/dtu-qmcm/fluxomics_data_model.git
cd fluxomics_data_model

# Install with uv
uv sync
```

### Requirements

- Python >= 3.12
- Key dependencies: pydantic, jax, sympy, sbmlmath

## Quick Start

### Parsing files

The library provides parsers for multiple formats. All parsers return a `FluxomicsDataModel` object:

```python
from fluxomics_data_model.io import parse_fluxml_file, parse_mtf, parse_freeflux

# Parse a FluxML model (13CFlux2/13CFlux3)
model = parse_fluxml_file("path/to/model.fml")

# Parse an MTF model (influx_si) - provide base path without extension
model = parse_mtf("path/to/model")  # reads .netw, .linp, .miso, etc.

# Parse a Freeflux model - provide base path (directory or prefix)
model = parse_freeflux("path/to/freeflux_dir")  # reads reactions.tsv, fluxes.tsv, etc.

# Access model components
print(f"Model: {model.info.name}")
print(f"Metabolites: {len(model.model.metabolites)}")
print(f"Reactions: {len(model.model.reactions)}")
print(f"Experiments: {model.experiments_names}")
```

### Accessing model data

```python
# Get reaction by ID
reaction = model.model.reactions.get_by_id("PGI")
print(f"Reaction: {reaction.id}")
print(f"Reversible: {reaction.reversible}")

# Get atom mapping for a reaction
atom_mapping = model.model.atom_mappings["PGI"]
print(atom_mapping.to_letter_notation())

# Access experimental data
for experiment in model.experiments:
    print(f"Experiment: {experiment.name}")
    print(f"Stationary: {experiment.stationary}")
```

### Working with constraints

```python
from fluxomics_data_model.model import ConstraintEvaluator

# Create evaluator with reaction IDs
reaction_ids = model.model.reactions.ids
evaluator = ConstraintEvaluator(reaction_ids, parameters={"mu": 0.03})

# Parse and evaluate constraints
constraint_fn, operator, rhs = evaluator.parse_formula("uptGLC = 1.0")
```

## Project Structure

```
fluxomics_data_model/
├── src/fluxomics_data_model/
│   ├── core/                    # Core data structures
│   │   ├── core.py              # FluxomicsDataModel, Metadata, Model, Experiments
│   │   └── common.py            # Shared utilities (DictList, Annotation, etc.)
│   ├── model/                   # Metabolic network components
│   │   ├── metabolite.py        # Metabolite definitions with atomic composition
│   │   ├── reaction.py          # Reaction definitions with variant support
│   │   ├── atom_mapping.py      # Atom mapping classes (AtomMapping, AtomMap)
│   │   ├── constraint.py        # Constraint definitions
│   │   └── constraint_eval.py   # JAX-compatible constraint evaluation
│   ├── experiment/              # Experimental data structures
│   │   ├── tracer.py            # Tracer and label composition
│   │   └── measurement.py       # Measurement data types
│   ├── output/                  # Output and simulation structures
│   │   └── simulation.py        # Simulation settings, variables, bounds
│   └── io/                      # Input/Output handlers
│       ├── fluxml_parser.py     # FluxML format parser
│       ├── fluxml_writer.py     # FluxML format writer
│       ├── mtf_parser.py        # MTF (influx_si) format parser
│       ├── mtf_writer.py        # MTF (influx_si) format writer
│       └── freeflux_parser.py   # Freeflux tabular format parser
├── tests/                       # Test suite
│   └── test_unit/               # Unit tests
├── data/                        # Example data files
│   ├── benchmark/               # Cross-format benchmark models
│   │   ├── FluxML/              # FluxML benchmark models
│   │   ├── FreeFlux/            # FreeFlux benchmark models
│   │   ├── INCA/                # INCA MATLAB benchmark models
│   │   └── influx_si/           # influx_si MTF benchmark models
│   ├── MTF_tests/               # MTF format test cases
│   └── specifications/          # Format specification documents
└── docs/                        # Documentation
```

## Data Model Architecture

The library uses a hierarchical data model built on [Pydantic](https://docs.pydantic.dev/) for validation and serialization:

```
FluxomicsDataModel
├── info: Metadata (name, version, modeler, strain, etc.)
├── model: Model
│   ├── metabolites: DictList[Metabolite]
│   ├── reactions: DictList[Reaction]
│   ├── atom_mappings: Dict[str, AtomMapping]
│   ├── compartments: List[str]
│   └── parameters: Dict[str, float]     # (planned) Growth rate, biomass coefficients
├── constraints: Constraints (net, exchange, metabolite size)
└── experiments: List[Experiments]
    ├── tracers: List[Tracers]
    ├── measurement: Measurement (labeling, flux, pool sizes)
    └── simulation: Simulation (variables, bounds)
```

## Features

### Variant Atom Mapping Support

Handles symmetric metabolites with multiple possible atom mappings:

```python
# Reactions with symmetric compounds generate variants
atom_mapping = model.model.atom_mappings["SCS"]  # Succinate symmetry
print(atom_mapping.maps.keys())  # ['SCS___1', 'SCS___2']
print(atom_mapping.weights)      # {'SCS___1': 0.5, 'SCS___2': 0.5}
```

### JAX Integration

Built for high-performance numerical computing with JAX:

```python
import jax.numpy as jnp

# Constraint evaluation is JAX-differentiable
from jax import grad

def loss_fn(flux_vector):
    residuals = evaluator.evaluate_constraints(flux_vector)
    return jnp.sum(residuals**2)

gradient = grad(loss_fn)(flux_vector)
```

### Isotopomer Transformations

Transform isotopomer distributions through reactions:

```python
# Define reactant isotopomer distribution
reactant_dist = {"A": jnp.array([0.99, 0.01, 0.0, 0.0])}  # Natural abundance

# Transform through reaction
product_dist = atom_mapping.transform_isotopomers(
    reactant_dist,
    atom_map_id="RXN___1"
)
```

## Development

### Running tests

```bash
uv run pytest tests/
```

### Running tests with coverage

```bash
uv run pytest tests/ --cov=fluxomics_data_model --cov-report=html
```

### Linting

```bash
uv run ruff check src/
uv run ruff format src/
```

### Type checking

```bash
uv run mypy src/
```

## Roadmap: Missing Features & Implementation Strategy

### 1. Freeflux Tabular Parser

**Status**: ✅ Read implemented (with caveats)

The Freeflux parser supports reading tabular models (TSV/CSV/XLSX):
- `reactions.{tsv,csv,xlsx}` - Network definition with atom mappings (required)
- `fluxes.{tsv,csv,xlsx}` - Flux values (optional)
- `concentrations.{tsv,csv,xlsx}` - Metabolite pool sizes (optional)
- `measured_MDVs.{tsv,csv,xlsx}` - Steady-state mass distribution vectors (optional)
- `measured_fluxes.{tsv,csv,xlsx}` - Measured flux values with uncertainties (optional)
- `measured_inst_MDVs.{tsv,csv,xlsx}` - Time-course MDV measurements (optional)

**Known limitations**:
- No tracers are parsed (FreeFlux specifies tracers via Python API, not files)
- No constraints are parsed (FreeFlux defines constraints in code)
- Metabolite atom counts are not available from file data (must be inferred from atom mappings)
- Export to FreeFlux format not yet implemented

---

### 2. SBML Import/Export

**Status**: Dependency installed (`sbmlmath>=0.4.0`), not integrated

**Required functionality**:
- Parse SBML Level 2/3 models
- Extract reaction stoichiometry and metabolite definitions
- Handle compartments and species references
- Note: SBML does not include atom mappings (requires separate annotation)

**Implementation strategy**:
```python
# Proposed module: src/fluxomics_data_model/io/sbml_parser.py
from sbmlmath import SBMLMathModel

class SBMLParser:
    def parse(self, filepath: Path) -> FluxomicsDataModel
    def _extract_metabolites(self, sbml_model) -> DictList[Metabolite]
    def _extract_reactions(self, sbml_model) -> DictList[Reaction]
```

**Priority**: Medium - useful for importing existing genome-scale models

---

### 3. MTF (influx_si) Import/Export

**Status**: ✅ Implemented (read and write)

The MTF parser supports reading and writing multi-file models for influx_si:
- `.netw` - Network definition with reactions and atom mappings
- `.linp` - Label input (tracer specifications, binary isotopomer patterns)
- `.miso` - MS isotopomer measurements (MS mass isotopomers and NMR/cumomer patterns)
- `.mflux` - Flux measurements
- `.mmet` - Metabolite concentration measurements
- `.cnstr` - Constraints (NET and XCH)
- `.tvar` - Variable types and starting values (NET, XCH, METAB with F/D/C classification)
- `.opt` - Options/command arguments (parsed but not stored in model)

**Usage**:
```python
from fluxomics_data_model.io import parse_mtf, write_mtf

# Read
model = parse_mtf("path/to/model")

# Write
write_mtf(model, "path/to/output")
```

**Remaining work**:
- Handle FTBL format (older influx_si format)
- Better handling of biomass parameters (see section 8)
- `.mmet` writer support

---

### 4. INCA Import/Export

**Status**: Not implemented

**Required functionality**:
- Parse INCA MATLAB `.m` files (reaction network, tracers, MS data, pool sizes, symmetric metabolites)
- Generate INCA-compatible MATLAB structures
- Export measurement data for INCA
- Integration with [incawrapper](https://github.com/biosustain/incawrapper)

**Implementation strategy**:
```python
# Proposed modules:
# src/fluxomics_data_model/io/inca_parser.py
class INCAParser:
    def parse(self, filepath: Path) -> FluxomicsDataModel
    def _parse_reactions(self, m_content: str) -> DictList[Reaction]
    def _parse_tracers(self, m_content: str) -> List[Tracers]
    def _parse_ms_data(self, m_content: str) -> Measurement

# src/fluxomics_data_model/io/inca_writer.py
class INCAWriter:
    def to_inca_model(self, model: FluxomicsDataModel) -> dict
    def write_mat(self, model: FluxomicsDataModel, filepath: Path) -> None
    def to_incawrapper(self, model: FluxomicsDataModel) -> "INCAModel"
```

**Priority**: High - INCA is widely used in the field

---

### 5. FluxML Writer

**Status**: ✅ Implemented

The FluxML writer exports models to FluxML XML format:

**Usage**:
```python
from fluxomics_data_model.io import write_fluxml

# Write model to FluxML format
write_fluxml(model, "path/to/output.fml")
```

**Supported elements**:
- Metadata (info element)
- Reaction network with atom mappings
- Constraints (NET, XCH, metabolite size)
- Experimental configurations

---

### 6. Command-Line Interface (CLI)

**Status**: Not implemented

**Required functionality**:
- Convert between formats from command line
- Validate model files
- Display model summaries
- Batch conversion

**Implementation strategy**:
```python
# Proposed module: src/fluxomics_data_model/cli.py
import click

@click.group()
def cli():
    """Fluxomics Data Model - Universal format converter"""
    pass

@cli.command()
@click.argument("input_file")
@click.argument("output_file")
@click.option("--input-format", "-i", type=click.Choice(["fluxml", "sbml", "ftbl", "csv"]))
@click.option("--output-format", "-o", type=click.Choice(["fluxml", "ftbl", "inca", "csv"]))
def convert(input_file, output_file, input_format, output_format):
    """Convert between fluxomics file formats"""
    pass

@cli.command()
@click.argument("input_file")
def validate(input_file):
    """Validate a fluxomics model file"""
    pass

@cli.command()
@click.argument("input_file")
def info(input_file):
    """Display model summary information"""
    pass
```

**Priority**: Medium - improves usability for non-programmers

---

### 7. Model Validation & Consistency Checks

**Status**: Basic validation via Pydantic, no semantic validation

**Required functionality**:
- Validate atom balance in reactions
- Check stoichiometric matrix rank
- Validate measurement references against model
- Check constraint feasibility

**Implementation strategy**:
```python
# Proposed module: src/fluxomics_data_model/validation/
class ModelValidator:
    def validate_atom_balance(self, model: Model) -> List[ValidationError]
    def validate_stoichiometry(self, model: Model) -> List[ValidationError]
    def validate_measurements(self, model: FluxomicsDataModel) -> List[ValidationError]
    def validate_all(self, model: FluxomicsDataModel) -> ValidationReport
```

**Priority**: High - prevents silent errors in analysis

---

### 8. Biomass & Growth Rate Parameters

**Status**: Not implemented

**Background**:
Different formats handle biomass and growth rate differently:
- **FluxML**: Biomass reactions (`G6P_bm___1`, `bmALA`) are proper reactions with annotations (`pathway="Biomass"`). Growth rate (`mu`) is a flux used in constraints.
- **MTF**: `BM` is a parameter used in constraints like `bs_glc6P-0.2050*BM==0`. Individual `bs_*` reactions drain metabolites to biomass.

**Required functionality**:
- Add `parameters` field to `Model` for storing growth-related parameters (`mu`, `BM`)
- Parse biomass parameters from MTF `.tvar` files
- Link biomass reactions to growth rate in constraints
- Identify biomass reactions by pathway annotation in FluxML

**Implementation strategy**:
```python
# In Model class
parameters: Dict[str, float] = Field(
    default_factory=dict,
    description="Model parameters (growth rate, biomass coefficients)"
)

# Special parameter names
# - "mu" or "BM": Growth rate
# - "bm_*": Biomass composition coefficients
```

**Priority**: Medium - important for growth-coupled flux analysis

---

### Implementation Priority Summary

| Feature | Priority | Status | Dependencies |
|---------|----------|--------|--------------|
| FluxML Parser | High | ✅ Done | - |
| FluxML Writer | High | ✅ Done | - |
| MTF Parser | High | ✅ Done | - |
| MTF Writer | High | ✅ Done | - |
| Freeflux Parser | High | ✅ Done (caveats) | pandas (optional) |
| INCA Parser | High | 🔲 Pending | - |
| INCA Writer | High | 🔲 Pending | incawrapper (optional) |
| Model Validation | High | 🔲 Pending | None |
| Stoichiometric coefficients | High | 🔲 Pending | None |
| Exchange flux parsing | Medium | 🔲 Pending | None |
| ErrorModel population | Medium | 🔲 Pending | None |
| Biomass Parameters | Medium | 🔲 Pending | None |
| SBML Parser | Medium | 🔲 Pending | sbmlmath |
| CLI | Medium | 🔲 Pending | click |
| Freeflux Writer | Low | 🔲 Pending | pandas |

### Known Issues

- **36 failing unit tests** in `tests/test_unit/` — primarily in `test_freeflux_parser.py` (20 failures), `test_roundtrip.py` (6 failures), `test_atom_mapping.py` (6 failures)
- **Pydantic V1 deprecation** — all models use `class Config: frozen = True` instead of `model_config = ConfigDict(frozen=True)`
- **Variant reactions with reactant+product symmetry** — fixed for FluxML parser, but FreeFlux comma-separated symmetric notation and MTF symmetric products still only generate single-mapping AtomMappings
- **Exchange flux measurements** — `ExchangeFlux` class exists but no parser populates `xch_fluxes`
- **ErrorModel** — FluxML parser reads error models but does not store them
- **FreeFlux models lack tracers** — no `labeling_strategy` file to read from; users must add tracers programmatically
- **Stoichiometric coefficients not parsed** — biomass reactions like `0.488*Ala → biomass` have incorrect stoichiometry
- **No metabolite-level SymmetryDefinition** — symmetry only represented via reaction variants (FluxML-style), not as a property of the metabolite itself (INCA-style)

See [`.agent/analysis.md`](.agent/analysis.md) for the full capability analysis and [`.agent/change_log.md`](.agent/change_log.md) for detailed change history.

## Contributing

We welcome contributions! Please see [CONTRIBUTING.md](docs/contributing.md) for guidelines.

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/new-parser`)
3. Make your changes with tests
4. Run the test suite (`uv run pytest`)
5. Submit a pull request

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Credits

### Authors

- [Quantitative Modelling of Cell Metabolism (QMCM)](https://github.com/dtu-qmcm) - DTU Bioengineering

### Acknowledgments

- [13CFlux](https://www.13cflux.net/) team for the FluxML format specification
- [influx_si](https://metasys.insa-toulouse.fr/software/influx/) developers
- [incawrapper](https://github.com/biosustain/incawrapper) team

### Citing

If you use this software in your research, please cite:

```bibtex
@software{fluxomics_data_model,
  author = {QMCM Team},
  title = {Fluxomics Data Model: Universal format converter for 13C-MFA},
  year = {2025-2026},
  url = {https://github.com/dtu-qmcm/fluxomics_data_model}
}
```

## Related Projects

- [incawrapper](https://github.com/biosustain/incawrapper) - Python wrapper for INCA
- [influx_si](https://metasys.insa-toulouse.fr/software/influx/) - 13C flux analysis software
- [13CFlux2](https://www.13cflux.net/) - Comprehensive 13C-MFA platform
- [COBRApy](https://github.com/opencobra/cobrapy) - Constraint-based metabolic modeling
