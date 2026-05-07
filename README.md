# NMR Calculator

NMR Calculator is a Python project for future 1H NMR, proton NMR, and hydrogen NMR spectrum prediction from molecular structures.

## Phase 4 Status

Phase 4 adds basic rule-based 1H NMR chemical shift estimates for the equivalent proton environment groups created in Phase 3. The program now parses SMILES with RDKit, adds explicit hydrogens, groups proton environments, and assigns approximate ppm values and ppm ranges.

These predictions are simple organic chemistry rule-based estimates. They are not database-backed, ML-based, or validated against nmrshiftdb2 yet, and plotted spectra are not generated in this phase.

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
python -m nmr_calculator.cli predict "c1ccccc1"
python -m nmr_calculator.cli predict "CC(=O)C"
```

Example ethanol prediction output:

```text
Input SMILES: CCO
Predicted 1H NMR chemical shifts:
  Group 1: atoms [0], 3H, alkyl CH3, ~1.2 ppm (0.8-1.8 ppm)
  Group 2: atoms [1], 2H, heteroatom-adjacent alkyl proton, ~3.5 ppm (3.0-4.5 ppm)
  Group 3: atoms [2], 1H, alcohol/amine/thiol exchangeable proton, ~2.5 ppm (0.5-5.5 ppm)
```

Inspect proton-bearing atoms:

```bash
python -m nmr_calculator.cli inspect "CCO"
```

Group equivalent proton environments:

```bash
python -m nmr_calculator.cli groups "CCO"
python -m nmr_calculator.cli groups "c1ccccc1"
python -m nmr_calculator.cli groups "CC(=O)C"
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
