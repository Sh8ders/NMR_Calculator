# NMR Calculator

NMR Calculator is a Python project for future 1H NMR, proton NMR, and hydrogen NMR spectrum prediction from molecular structures.

## Phase 5 Status

Phase 5 converts rule-based 1H NMR shift predictions into a clean peak list. The program now parses SMILES with RDKit, groups equivalent proton environments, assigns approximate ppm values and ranges, and emits one sorted peak per predicted proton group.

Multiplicity is currently a placeholder:

- exchangeable proton environments: `br s`
- aromatic proton environments: `m`
- all other environments: `unknown`

True spin-spin splitting will be improved in a later phase. These predictions remain simple rule-based estimates; they are not database-backed, ML-based, or validated against nmrshiftdb2 yet. Spectrum plotting is not implemented in this phase.

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
python -m nmr_calculator.cli peaks "c1ccccc1"
python -m nmr_calculator.cli peaks "CC(=O)C"
```

Example ethanol peak list:

```text
Input SMILES: CCO
Predicted 1H NMR peak list:
  Peak 1: 3.5 ppm, 2H, unknown, heteroatom-adjacent alkyl proton
  Peak 2: 2.5 ppm, 1H, br s, alcohol/amine/thiol exchangeable proton
  Peak 3: 1.2 ppm, 3H, unknown, alkyl CH3
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
