# NMR Calculator

NMR Calculator is a Python project for future 1H NMR, proton NMR, and hydrogen NMR spectrum prediction from molecular structures.

## Phase 11 Status

Phase 11 adds local 1H NMR shift database loading, validation, normalization, and CSV caching. This prepares the project for future nmrshiftdb2-style database prediction, but it does not perform database matching or replace the rule-based predictor yet.

Supported local CSV schema:

```text
smiles,molecule_name,atom_index,shift_ppm,assignment_label,source
```

The loader validates required columns, atom index types, numeric shifts, and parseable SMILES. nmrshiftdb2 SDF parsing and true database matching will come later.

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

Validate a local shift database CSV:

```bash
python -m nmr_calculator.cli database-check tests/fixtures/test_1h_shift_database.csv
```

Create a normalized database cache:

```bash
python -m nmr_calculator.cli database-cache tests/fixtures/test_1h_shift_database.csv --output data/processed/test_cache.csv
```

Use the current rule-based predictor:

```bash
python -m nmr_calculator.cli predict "CCO" --method rules
python -m nmr_calculator.cli peaks "CCO" --method rules
python -m nmr_calculator.cli plot "CCO" --output ethanol_1h_nmr.png --method rules
python -m nmr_calculator.cli report "CCO" --output-dir reports/ethanol --method rules
```

Generated reports, processed data caches, and PNGs should usually not be committed.

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
