# NMR Calculator

NMR Calculator is a Python project for future 1H NMR, proton NMR, and hydrogen NMR spectrum prediction from molecular structures.

## Phase 6 Status

Phase 6 generates a simple simulated 1H NMR spectrum plot from the predicted peak list. The program now parses SMILES with RDKit, groups equivalent proton environments, estimates ppm values, generates sorted peak lists, and saves basic PNG spectrum images.

The plot is still rule-based and approximate. Multiplicity and splitting remain placeholder-level; each predicted proton group is plotted as one vertical signal with height based on integration. Database, ML, nmrshiftdb2, and true spin-spin splitting support are not implemented yet.

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

Generate predicted peak lists:

```bash
python -m nmr_calculator.cli peaks "CCO"
```

Generate simulated spectrum PNGs:

```bash
python -m nmr_calculator.cli plot "CCO" --output ethanol_1h_nmr.png
python -m nmr_calculator.cli plot "c1ccccc1" --output benzene_1h_nmr.png
python -m nmr_calculator.cli plot "CC(=O)C" --output acetone_1h_nmr.png
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
