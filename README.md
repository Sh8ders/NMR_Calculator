# NMR Calculator

NMR Calculator is a Python project for future 1H NMR, proton NMR, and hydrogen NMR spectrum prediction from molecular structures.

## Phase 12 Status

Phase 12 adds exact-match local CSV database prediction. The `database` prediction method loads a validated local 1H NMR shift database, canonicalizes query and database SMILES with RDKit, finds exact molecule matches, and maps database atom-index assignments onto predicted proton environment groups.

This is still not full nmrshiftdb2 SDF matching. Similarity matching, SDF parsing, ML prediction, and hybrid fallback are not implemented yet.

Supported local CSV schema:

```text
smiles,molecule_name,atom_index,shift_ppm,assignment_label,source
```

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

Use exact local database prediction:

```bash
python -m nmr_calculator.cli predict "CCO" --method database --database tests/fixtures/test_1h_shift_database.csv
python -m nmr_calculator.cli peaks "CCO" --method database --database tests/fixtures/test_1h_shift_database.csv
python -m nmr_calculator.cli plot "CCO" --output ethanol_database.png --method database --database tests/fixtures/test_1h_shift_database.csv
python -m nmr_calculator.cli report "CCO" --output-dir reports/ethanol_database --method database --database tests/fixtures/test_1h_shift_database.csv
```

Validate or cache a local database:

```bash
python -m nmr_calculator.cli database-check tests/fixtures/test_1h_shift_database.csv
python -m nmr_calculator.cli database-cache tests/fixtures/test_1h_shift_database.csv --output data/processed/test_cache.csv
```

Use the rule-based predictor:

```bash
python -m nmr_calculator.cli predict "CCO" --method rules
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
- Phase 12: exact local database predictor
