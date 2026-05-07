# NMR Calculator

NMR Calculator is a Python project for future 1H NMR, proton NMR, and hydrogen NMR spectrum prediction from molecular structures.

## Phase 7 Status

Phase 7 adds RDKit molecule structure image generation. The program can now parse SMILES, inspect and group proton environments, estimate rule-based 1H NMR shifts, generate peak lists, save simple spectrum plots, and save molecule structure PNGs.

Molecule images include RDKit atom index labels by default. These labels are useful because the `predict` and `groups` commands list atom indices for proton environments.

The NMR predictions and plots remain approximate and rule-based. Multiplicity and splitting are still placeholder-level, and database, ML, nmrshiftdb2, and true spin-spin splitting support are not implemented yet.

Generated PNG files should usually not be committed to the repository unless they are intentionally added as documentation examples.

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

Predict approximate 1H NMR shifts:

```bash
python -m nmr_calculator.cli predict "CCO"
```

Generate predicted peak lists and spectrum PNGs:

```bash
python -m nmr_calculator.cli peaks "CCO"
python -m nmr_calculator.cli plot "CCO" --output ethanol_1h_nmr.png
```

Generate molecule structure PNGs:

```bash
python -m nmr_calculator.cli molecule-image "CCO" --output ethanol_structure.png
python -m nmr_calculator.cli molecule-image "c1ccccc1" --output benzene_structure.png
python -m nmr_calculator.cli molecule-image "CC(=O)C" --output acetone_structure.png
```

To hide atom indices:

```bash
python -m nmr_calculator.cli molecule-image "CCO" --output ethanol_structure.png --no-atom-indices
```

Inspect proton-bearing atoms and groups:

```bash
python -m nmr_calculator.cli inspect "CCO"
python -m nmr_calculator.cli groups "CCO"
```

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
