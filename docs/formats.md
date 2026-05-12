# Supported Formats

| Format | Software | Read | Write | Notes |
|--------|----------|------|-------|-------|
| FluxML | 13CFlux2, 13CFlux3 | ✅ | ✅ | Fully implemented |
| MTF | influx_si | ✅ | ✅ | Cumomer, `.mmet`, METAB vars supported |
| CSV/TSV/XLSX | FreeFlux | ✅ | 🔲 | No tracers/constraints from files |
| MATLAB | INCA | 🔲 | 🔲 | Planned |
| XML | SBML | 🔲 | 🔲 | Planned (dependency ready) |

## FluxML

XML-based format used by 13CFlux2/13CFlux3. Supports the richest
feature set including multiple experiments, NMR data, error models,
and variant atom transitions.

```python
from fluxomics_data_converter.io import parse_fluxml_file, write_fluxml

model = parse_fluxml_file("model.fml")
write_fluxml(model, "output.fml")
```

## MTF (influx_si)

Multi-file text format. Each model consists of files sharing a
basename:

| File | Description |
|------|-------------|
| `.netw` | Network with reactions and atom transitions |
| `.linp` | Tracer specifications |
| `.miso` | MS/NMR isotopomer measurements |
| `.mflux` | Flux measurements |
| `.mmet` | Metabolite concentration measurements |
| `.cnstr` | Constraints (NET and XCH) |
| `.tvar` | Variable types and starting values |
| `.opt` | Options |

```python
from fluxomics_data_converter.io import parse_mtf, write_mtf

model = parse_mtf("path/to/model")  # reads all files
write_mtf(model, "path/to/output", experiment_name="exp1")
```

## FreeFlux

Tabular format (TSV/CSV/XLSX). Parsed files:

- `reactions.{tsv,csv,xlsx}` — network definition (required)
- `fluxes.{tsv,csv,xlsx}` — flux values
- `concentrations.{tsv,csv,xlsx}` — metabolite pool sizes
- `measured_MDVs.{tsv,csv,xlsx}` — steady-state MDVs
- `measured_fluxes.{tsv,csv,xlsx}` — measured fluxes
- `measured_inst_MDVs.{tsv,csv,xlsx}` — time-course MDVs

!!! note
    FreeFlux specifies tracers and constraints via Python API, not
    files. Parsed models will have empty tracers and no constraints.

```python
from fluxomics_data_converter.io import parse_freeflux

model = parse_freeflux("path/to/freeflux_dir")
```
