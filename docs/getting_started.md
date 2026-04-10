# Getting Started

## Parsing files

All parsers return a `FluxomicsDataModel` object:

```python
from fluxomics_data_model.io import (
    parse_fluxml_file,
    parse_mtf,
    parse_freeflux,
)

# FluxML (13CFlux2/13CFlux3)
model = parse_fluxml_file("path/to/model.fml")

# MTF (influx_si) - provide base path without extension
model = parse_mtf("path/to/model")

# FreeFlux - provide directory containing reaction files
model = parse_freeflux("path/to/freeflux_dir")
```

## Accessing model data

```python
# Model info
print(f"Model: {model.info.name}")

# Reactions
reaction = model.model.reactions.get_by_id("PGI")
print(f"Reaction: {reaction.id}")
print(f"Reversible: {reaction.reversibility}")

# Atom mapping
atom_mapping = model.model.atom_mappings["PGI"]
print(atom_mapping.to_letter_notation())

# Experiments
for experiment in model.experiments:
    print(f"Experiment: {experiment.name}")
    print(f"Stationary: {experiment.stationary}")
```

## Writing files

```python
from fluxomics_data_model.io import write_fluxml, write_mtf

# Write to FluxML
write_fluxml(model, "path/to/output.fml")

# Write to MTF (specify experiment for multi-experiment models)
write_mtf(model, "path/to/output", experiment_name="exp1")
```

## Working with constraints

```python
from fluxomics_data_model.model import ConstraintEvaluator

reaction_ids = model.model.reactions.ids
evaluator = ConstraintEvaluator(reaction_ids, parameters={"mu": 0.03})

constraint_fn, operator, rhs = evaluator.parse_formula("uptGLC = 1.0")
```
