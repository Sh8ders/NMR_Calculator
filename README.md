# NMR Calculator

NMR Calculator is a Python project for future 1H NMR, proton NMR, and hydrogen NMR spectrum prediction from molecular structures.

## Phase 10 Status

Phase 10 refactors prediction behind a clean predictor interface. The project now has a rule-based predictor implementation plus placeholder database and hybrid predictor classes for future nmrshiftdb2 work.

Prediction methods:

- `rules`: current supported rule-based predictor
- `database`: placeholder for future nmrshiftdb2 database prediction
- `hybrid`: placeholder for future database-first prediction with rule fallback

Database and hybrid modes are not implemented yet. nmrshiftdb2 support will be added in a later phase. Current predictions remain approximate rule-based estimates.

## Installation

Use Python 3.11 or newer.

Create and activate a virtual environment:

```bash
python -m venv .venv
source .venv/bin/activate
```

Install the project:

```bash
python -m pip install -e .
```

For development and tests:

```bash
python -m pip install -e ".[dev]"
```

RDKit is included as the `rdkit` pip dependency in `pyproject.toml`. If RDKit installation fails on a specific platform, install it in a fresh virtual environment first:

```bash
python -m pip install rdkit
```

## Run the CLI

Use the current rule-based predictor explicitly:

```bash
python -m nmr_calculator.cli predict "CCO" --method rules
python -m nmr_calculator.cli peaks "CCO" --method rules
python -m nmr_calculator.cli plot "CCO" --output ethanol_1h_nmr.png --method rules
python -m nmr_calculator.cli report "CCO" --output-dir reports/ethanol --method rules
```

The default method is `rules`, so existing commands still work without `--method`.

Other useful commands:

```bash
python -m nmr_calculator.cli molecule-image "CCO" --output ethanol_structure.png
python -m nmr_calculator.cli inspect "CCO"
python -m nmr_calculator.cli groups "CCO"
```

Generated reports and PNGs should usually not be committed.

## Run Tests

```bash
python -m pytest
```

## Roadmap

- Phase 2: RDKit molecule parsing
- Phase 3: equivalent proton grouping
- Phase 4: basic rule-based 1H shift prediction
- Phase 5: peak list generation
- Phase 6: spectrum plotting
- Phase 7: molecule image generation
- Phase 8: report generation
- Phase 9: improved rules
- Phase 10: predictor interface
- Phase 11: nmrshiftdb2 database support
- Phase 12: hybrid database plus rule predictor
