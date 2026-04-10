# Contributing

We welcome contributions! Please see
[CONTRIBUTING.md](https://github.com/dtu-qmcm/fluxomics_data_model/blob/main/CONTRIBUTING.md)
for guidelines.

## Development setup

```bash
git clone https://github.com/dtu-qmcm/fluxomics_data_model.git
cd fluxomics_data_model
uv sync
```

## Running tests

```bash
uv run pytest tests/
```

## Linting

```bash
uv run ruff check src/
uv run ruff format src/
```
