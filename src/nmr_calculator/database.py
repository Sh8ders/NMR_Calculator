"""Local 1H NMR shift database loading and validation utilities."""

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from nmr_calculator.molecule import canonicalize_smiles, parse_molecule

DATABASE_COLUMNS = [
    "smiles",
    "molecule_name",
    "atom_index",
    "shift_ppm",
    "assignment_label",
    "source",
]

STRING_COLUMNS = [
    "smiles",
    "molecule_name",
    "assignment_label",
    "source",
]


@dataclass(frozen=True)
class NMRDatabaseRecord:
    smiles: str
    molecule_name: str
    atom_index: int
    shift_ppm: float
    assignment_label: str
    source: str


def load_shift_database_csv(path: str) -> pd.DataFrame:
    """Load and validate a local 1H NMR shift database CSV."""
    try:
        df = pd.read_csv(path)
    except Exception as error:
        raise ValueError(f"Could not load shift database CSV: {path}") from error
    return normalize_shift_database(df)


def add_canonical_smiles_column(df: pd.DataFrame) -> pd.DataFrame:
    """Return a copy with a canonical_smiles column added."""
    normalized = normalize_shift_database(df)
    try:
        normalized["canonical_smiles"] = normalized["smiles"].map(canonicalize_smiles)
    except ValueError as error:
        raise ValueError("Could not canonicalize database SMILES values.") from error
    return normalized


def filter_database_by_smiles(df: pd.DataFrame, smiles: str) -> pd.DataFrame:
    """Filter database records to exact canonical SMILES matches."""
    query_canonical_smiles = canonicalize_smiles(smiles)
    if "canonical_smiles" not in df.columns:
        database = add_canonical_smiles_column(df)
    else:
        database = df.copy()

    return database[
        database["canonical_smiles"] == query_canonical_smiles
    ].reset_index(drop=True)


def normalize_shift_database(df: pd.DataFrame) -> pd.DataFrame:
    """Normalize and validate a local 1H NMR shift database table."""
    missing_columns = [column for column in DATABASE_COLUMNS if column not in df.columns]
    if missing_columns:
        raise ValueError(
            "Shift database is missing required columns: "
            f"{', '.join(missing_columns)}."
        )

    normalized = df.dropna(how="all").copy()
    normalized = normalized[DATABASE_COLUMNS]

    for column in STRING_COLUMNS:
        normalized[column] = normalized[column].map(_clean_string_value)
        if normalized[column].eq("").any():
            raise ValueError(f"Shift database column {column!r} contains empty values.")

    normalized["atom_index"] = _convert_atom_index(normalized["atom_index"])
    normalized["shift_ppm"] = _convert_shift_ppm(normalized["shift_ppm"])

    for row_index, smiles in normalized["smiles"].items():
        try:
            mol = parse_molecule(smiles)
        except ValueError as error:
            raise ValueError(
                f"Invalid SMILES in shift database at row {row_index}: {smiles!r}."
            ) from error

        atom_index = int(normalized.at[row_index, "atom_index"])
        if atom_index < 0 or atom_index >= mol.GetNumAtoms():
            raise ValueError(
                "Invalid atom_index in shift database at "
                f"row {row_index}: {atom_index} is outside the molecule atom range."
            )

    return normalized.reset_index(drop=True)


def save_database_cache(df: pd.DataFrame, output_path: str) -> str:
    """Save a normalized 1H NMR shift database CSV cache."""
    normalized = normalize_shift_database(df)
    output_file = Path(output_path)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    normalized.to_csv(output_file, index=False)
    return str(output_file)


def load_database_cache(path: str) -> pd.DataFrame:
    """Load and validate a previously saved database cache CSV."""
    return load_shift_database_csv(path)


def load_nmrshiftdb2_sdf(path: str) -> pd.DataFrame:
    """Placeholder for future nmrshiftdb2 SDF parsing."""
    raise NotImplementedError(
        "SDF parsing for nmrshiftdb2 is not implemented yet. "
        "Use CSV loading for Phase 11."
    )


def _clean_string_value(value) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip()


def _convert_atom_index(values: pd.Series) -> pd.Series:
    try:
        converted = pd.to_numeric(values, errors="raise")
    except Exception as error:
        raise ValueError("Shift database contains an invalid atom_index value.") from error

    if converted.isna().any() or (converted % 1 != 0).any():
        raise ValueError("Shift database atom_index values must be integers.")
    return converted.astype(int)


def _convert_shift_ppm(values: pd.Series) -> pd.Series:
    try:
        converted = pd.to_numeric(values, errors="raise")
    except Exception as error:
        raise ValueError("Shift database contains an invalid shift_ppm value.") from error

    if converted.isna().any():
        raise ValueError("Shift database shift_ppm values must be numeric.")
    return converted.astype(float)
