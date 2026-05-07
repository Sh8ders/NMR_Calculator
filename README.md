# NMR Calculator

NMR Calculator is a Python project scaffold for a future 13C NMR chemical shift and spectrum prediction program.

## Phase 1 Status

Phase 1 is complete when the repository has a clean, testable Python package structure, a Typer command-line interface, placeholder modules, and passing tests.

This phase does not implement real chemistry prediction logic yet.

TODO: Add RDKit in Phase 2 for molecule parsing. RDKit is intentionally not required in Phase 1 to keep installation reliable.

## Installation

Use Python 3.11 or newer.

```bash
python -m pip install -e .
```

For development and tests:

```bash
python -m pip install -e ".[dev]"
```

## Run the CLI

```bash
python -m nmr_calculator.cli predict "CCO"
```

Expected Phase 1 output:

```text
Input SMILES: CCO
Phase 1 project setup complete.
13C NMR prediction will be implemented in later phases.
```

## Run Tests

```bash
python -m pytest
```

## Roadmap

- Phase 2: RDKit molecule parsing
- Phase 3: equivalent carbon grouping
- Phase 4: basic rule-based 13C shift prediction
- Phase 5: peak list generation
- Phase 6: spectrum plotting
- Phase 7: molecule image generation
- Phase 8: report generation
- Phase 9: improved rules
- Phase 10: predictor interface
- Phase 11: nmrshiftdb2 database support
- Phase 12: hybrid database plus rule predictor
