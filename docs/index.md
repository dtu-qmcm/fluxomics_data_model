# Fluxomics Data Model

A universal Python library for converting and representing fluxomics data
across multiple formats. It enables seamless interoperability between
different 13C metabolic flux analysis (13C-MFA) software tools.

## Features

- **Parse** models from FluxML, FreeFlux, and MTF (influx_si) formats
- **Convert** between formats
- **Export** to formats compatible with major 13C-MFA tools
- **Validate** model data via Pydantic schemas
- **JAX-compatible** constraint evaluation for numerical computing

## Quick example

```python
from fluxomics_data_model.io import parse_fluxml_file, parse_mtf

# Parse a FluxML model
model = parse_fluxml_file("path/to/model.fml")

# Parse an MTF model (influx_si)
model = parse_mtf("path/to/model")  # reads .netw, .linp, .miso, etc.

# Access model components
print(f"Metabolites: {len(model.model.metabolites)}")
print(f"Reactions: {len(model.model.reactions)}")
print(f"Experiments: {model.experiments_names}")
```

## Installation

```bash
pip install fluxomics-data-model
```

Or for development:

```bash
git clone https://github.com/dtu-qmcm/fluxomics_data_model.git
cd fluxomics_data_model
uv sync
```
