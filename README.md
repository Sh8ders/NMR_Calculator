# NMR Calculator

NMR Calculator is a Python project for future 1H NMR, proton NMR, and hydrogen NMR spectrum prediction from molecular structures.

## Phase 12 Status

Phase 12 adds exact-match local CSV database prediction. The `database` prediction method loads a validated local 1H NMR shift database, canonicalizes query and database SMILES with RDKit, finds exact molecule matches, and maps database atom-index assignments onto predicted proton environment groups.

The project can import assigned 1H shifts from nmrshiftdb2 SD/NMReDATA-style files into this local CSV format. Similarity matching, ML prediction, and hybrid fallback are not implemented yet.

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

## Using nmrshiftdb2 data

The full nmrshiftdb2 database file is not committed to this repository. Download it locally, import assigned 1H records into the supported CSV schema, then use exact database prediction with that processed CSV:

```bash
python -m nmr_calculator.cli nmrshiftdb2-download --output data/raw/nmrshiftdb2rawdata.nmredata.sd
python -m nmr_calculator.cli nmrshiftdb2-import data/raw/nmrshiftdb2rawdata.nmredata.sd --output data/processed/nmrshiftdb2_1h.csv
python -m nmr_calculator.cli predict "CCO" --method database --database data/processed/nmrshiftdb2_1h.csv
```

The default download URL is:

```text
https://sourceforge.net/projects/nmrshiftdb2/files/data/nmrshiftdb2rawdata.nmredata.sd/download
```

An alternative snapshot URL is:

```text
https://sourceforge.net/p/nmrshiftdb2/code/HEAD/tree/trunk/snapshots/nmrshiftdb2withsignals.sd?format=raw
```

Exact database prediction only works when the molecule exists in the imported database. Parsing quality depends on whether atom assignments are present in the raw file; records without both an atom index and a 1H shift are skipped.

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
