# NMR Calculator

NMR Calculator is a Python project for future 1H NMR, proton NMR, and hydrogen NMR spectrum prediction from molecular structures.

## Phase 9 Status

Phase 9 improves the rule-based 1H NMR predictor with better functional group and proton environment detection. The current rules distinguish simple alkyl, heteroatom-adjacent, alpha-to-carbonyl, benzylic, allylic, vinylic, aromatic, aldehyde, terminal alkyne, alcohol, amine, thiol, and carboxylic acid proton environments.

Predictions are still approximate rule-based estimates. Database prediction, ML prediction, nmrshiftdb2 support, improved multiplicity, and true spin-spin splitting are not implemented yet.

The existing workflow still supports prediction tables, peak lists, spectrum PNGs, molecule PNGs, and complete report folders.

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

Try improved rule-based predictions:

```bash
python -m nmr_calculator.cli predict "Cc1ccccc1"
python -m nmr_calculator.cli predict "CC=O"
python -m nmr_calculator.cli predict "CC(=O)O"
```

Generate a complete report folder:

```bash
python -m nmr_calculator.cli report "Cc1ccccc1" --output-dir reports/toluene
```

Other useful commands:

```bash
python -m nmr_calculator.cli peaks "CCO"
python -m nmr_calculator.cli plot "CCO" --output ethanol_1h_nmr.png
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
