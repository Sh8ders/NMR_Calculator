# NMR Calculator

NMR Calculator is a Python project for future 1H NMR, proton NMR, and hydrogen NMR spectrum prediction from molecular structures.

## Phase 2 Status

Phase 2 adds RDKit-based molecule parsing from SMILES, explicit hydrogen preparation, hydrogen-bearing atom inspection, and basic proton environment labels.

This phase does not predict proton chemical shifts yet.

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

Check Phase 2 prediction readiness:

```bash
python -m nmr_calculator.cli predict "CCO"
```

Expected Phase 2 output:

```text
Input SMILES: CCO
Phase 2 molecule parsing ready.
1H NMR prediction will be implemented in later phases.
```

Inspect proton environments:

```bash
python -m nmr_calculator.cli inspect "CCO"
```

Example output:

```text
Input SMILES: CCO
Hydrogen-bearing atoms:
  Atom 0 (CH3): alkyl CH3
  Atom 1 (CH2): heteroatom-adjacent alkyl proton
  Atom 2 (OH1): alcohol/amine/thiol exchangeable proton
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
