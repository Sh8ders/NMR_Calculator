"""Local 1H NMR shift database loading and validation utilities."""

from dataclasses import dataclass
from pathlib import Path
import re
from urllib import request
import warnings

import pandas as pd
from rdkit import Chem

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

DEFAULT_NMRSHIFTDB2_URL = (
    "https://sourceforge.net/projects/nmrshiftdb2/files/data/"
    "nmrshiftdb2rawdata.nmredata.sd/download"
)

LIKELY_1H_NMR_PROPERTY_NAMES = {
    "NMREDATA_1D_1H",
    "NMREDATA_ASSIGNMENT",
    "Spectrum 1H",
    "1H NMR",
    "nmrshiftdb2 1H",
}
LIKELY_1H_NMR_PROPERTY_NAMES_UPPER = {
    property_name.upper() for property_name in LIKELY_1H_NMR_PROPERTY_NAMES
}

MOLECULE_NAME_PROPERTIES = [
    "_Name",
    "NAME",
    "Name",
    "Molecule Name",
    "molecule_name",
    "TITLE",
    "Title",
]


@dataclass(frozen=True)
class NMRDatabaseRecord:
    smiles: str
    molecule_name: str
    atom_index: int
    shift_ppm: float
    assignment_label: str
    source: str


@dataclass(frozen=True)
class NMRShiftDB2ImportStats:
    molecules_scanned: int = 0
    molecules_with_1h_like_properties: int = 0
    records_extracted: int = 0
    skipped_records: int = 0


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
    """Load nmrshiftdb2 SD/NMReDATA-style 1H shifts as a normalized DataFrame."""
    database, stats = load_nmrshiftdb2_sdf_with_stats(path)
    database.attrs["nmrshiftdb2_import_stats"] = stats
    return database


def load_nmrshiftdb2_sdf_with_stats(
    path: str,
) -> tuple[pd.DataFrame, NMRShiftDB2ImportStats]:
    """Load nmrshiftdb2 1H shifts and return import diagnostics."""
    supplier = Chem.SDMolSupplier(str(path), removeHs=False)
    records = []
    molecules_scanned = 0
    molecules_with_1h_like_properties = 0
    skipped_records = 0
    for mol_index, mol in enumerate(supplier):
        molecules_scanned += 1
        if mol is None:
            warnings.warn(
                f"Skipping unparseable molecule at SDF record {mol_index}.",
                stacklevel=2,
            )
            continue
        extraction = extract_1h_shift_records_from_mol_with_stats(mol)
        records.extend(extraction["records"])
        skipped_records += extraction["skipped_records"]
        if extraction["has_1h_like_properties"]:
            molecules_with_1h_like_properties += 1

    stats = NMRShiftDB2ImportStats(
        molecules_scanned=molecules_scanned,
        molecules_with_1h_like_properties=molecules_with_1h_like_properties,
        records_extracted=len(records),
        skipped_records=skipped_records,
    )

    if not records:
        raise ValueError(
            f"No assigned 1H shift records could be extracted from {path}. "
            "Inspect the SD properties with: python -m nmr_calculator.cli "
            f"nmrshiftdb2-inspect {path} --max-molecules 5"
        )

    database = normalize_shift_database(pd.DataFrame(records, columns=DATABASE_COLUMNS))
    database.attrs["nmrshiftdb2_import_stats"] = stats
    return database, stats


def inspect_sdf_properties(path: str, max_molecules: int = 5) -> list[dict]:
    """Return property names and short previews for the first molecules in an SD file."""
    supplier = Chem.SDMolSupplier(str(path), removeHs=False)
    inspected = []
    for mol_index, mol in enumerate(supplier):
        if mol_index >= max_molecules:
            break
        if mol is None:
            inspected.append(
                {
                    "molecule_index": mol_index,
                    "molecule_name": "",
                    "canonical_smiles": "",
                    "property_names": [],
                    "property_previews": {},
                }
            )
            continue

        property_names = list(
            mol.GetPropNames(includePrivate=True, includeComputed=False)
        )
        inspected.append(
            {
                "molecule_index": mol_index,
                "molecule_name": _get_molecule_name(mol),
                "canonical_smiles": _mol_to_canonical_smiles(mol),
                "property_names": property_names,
                "property_previews": {
                    property_name: _preview_sdf_property(mol.GetProp(property_name))
                    for property_name in property_names
                },
            }
        )
    return inspected


def extract_1h_shift_records_from_mol(
    mol, source: str = "nmrshiftdb2"
) -> list[dict]:
    """Extract assigned 1H shift records from one RDKit molecule."""
    return extract_1h_shift_records_from_mol_with_stats(mol, source=source)["records"]


def extract_1h_shift_records_from_mol_with_stats(
    mol, source: str = "nmrshiftdb2"
) -> dict:
    """Extract assigned 1H shift records and skip counts from one RDKit molecule."""
    if mol is None:
        return {
            "records": [],
            "skipped_records": 0,
            "has_1h_like_properties": False,
        }

    smiles = _mol_to_database_smiles(mol)
    molecule_name = _get_molecule_name(mol)
    records = []
    seen = set()
    skipped_records = 0
    one_h_signal_labels = _find_nmredata_1h_signal_labels(mol)
    property_names = find_1h_nmr_property_names(mol)

    for property_name in property_names:
        text = mol.GetProp(property_name)
        parsed_records, skipped_count = _parse_property_1h_assignments(
            mol, property_name, text, one_h_signal_labels
        )
        skipped_records += skipped_count
        for atom_index, shift_ppm, assignment_label in parsed_records:
            if atom_index < 0 or atom_index >= mol.GetNumAtoms():
                warnings.warn(
                    "Skipping 1H shift assignment with atom_index outside "
                    f"molecule atom range: {atom_index}.",
                    stacklevel=2,
                )
                skipped_records += 1
                continue
            key = (atom_index, shift_ppm, assignment_label)
            if key in seen:
                continue
            seen.add(key)
            records.append(
                {
                    "smiles": smiles,
                    "molecule_name": molecule_name,
                    "atom_index": atom_index,
                    "shift_ppm": shift_ppm,
                    "assignment_label": assignment_label or f"atom {atom_index}",
                    "source": source,
                }
            )

    return {
        "records": records,
        "skipped_records": skipped_records,
        "has_1h_like_properties": bool(property_names),
    }


def find_1h_nmr_property_names(mol) -> list[str]:
    """Return likely SD property names containing assigned 1H NMR data."""
    property_names = list(mol.GetPropNames(includePrivate=True, includeComputed=False))
    matches = []
    for property_name in property_names:
        upper_name = property_name.upper()
        if upper_name in LIKELY_1H_NMR_PROPERTY_NAMES_UPPER:
            matches.append(property_name)
        elif "1H" in upper_name and "NMR" in upper_name:
            matches.append(property_name)
        elif "NMREDATA" in upper_name and "1H" in upper_name:
            matches.append(property_name)
        elif upper_name == "NMREDATA_ASSIGNMENT":
            matches.append(property_name)
    return matches


def parse_shift_assignment_lines(text: str) -> list[tuple[int, float, str]]:
    """Parse flexible atom/shift/label line formats."""
    records = []
    for raw_line in text.splitlines():
        parsed = _parse_shift_assignment_line(raw_line)
        if parsed is not None:
            records.append(parsed)
    return records


def parse_nmredata_assignment(text: str) -> list[tuple[int, float, str]]:
    """Parse NMReDATA-style assignment text using the flexible line parser."""
    records = []
    for raw_line in text.splitlines():
        for line in raw_line.split("|"):
            parsed = _parse_shift_assignment_line(line)
            if parsed is not None:
                records.append(parsed)
    return records


def parse_nmredata_signal_assignments(
    text: str,
    mol=None,
    one_h_signal_labels: set[str] | None = None,
) -> tuple[list[tuple[int, float, str]], int]:
    """Parse real nmrshiftdb2 NMREDATA_ASSIGNMENT signal-label assignments."""
    records = []
    skipped_records = 0
    for raw_line in text.splitlines():
        clean_line = _clean_nmredata_line(raw_line)
        if not clean_line:
            continue
        match = re.match(
            r"^(?P<label>s\d+)\s*,\s*(?P<shift>-?\d+(?:\.\d+)?)\s*(?:,\s*(?P<atoms>.*))?$",
            clean_line,
            flags=re.IGNORECASE,
        )
        if not match:
            continue

        signal_label = match.group("label").lower()
        if one_h_signal_labels and signal_label not in one_h_signal_labels:
            continue

        atom_tokens = re.findall(r"\d+", match.group("atoms") or "")
        if not atom_tokens:
            skipped_records += 1
            continue

        shift_ppm = float(match.group("shift"))
        atom_indices = []
        for atom_token in atom_tokens:
            atom_index = _resolve_nmredata_atom_index(int(atom_token), mol)
            if atom_index is None:
                skipped_records += 1
                continue
            atom_indices.append(atom_index)

        for atom_index in sorted(set(atom_indices)):
            records.append((atom_index, shift_ppm, signal_label))

    return records, skipped_records


def download_nmrshiftdb2_data(output_path: str, url: str | None = None) -> str:
    """Download the public nmrshiftdb2 SD/NMReDATA file to a local path."""
    source_url = url or DEFAULT_NMRSHIFTDB2_URL
    output_file = Path(output_path)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    try:
        with request.urlopen(source_url) as response:
            output_file.write_bytes(response.read())
    except Exception as error:
        raise RuntimeError(
            f"Could not download nmrshiftdb2 data from {source_url}."
        ) from error

    return str(output_file)


def _parse_property_1h_assignments(
    mol, property_name: str, text: str, one_h_signal_labels: set[str]
) -> tuple[list[tuple[int, float, str]], int]:
    upper_name = property_name.upper()
    if upper_name == "NMREDATA_ASSIGNMENT":
        return parse_nmredata_signal_assignments(
            text, mol=mol, one_h_signal_labels=one_h_signal_labels
        )
    if "NMREDATA" in upper_name:
        return parse_nmredata_assignment(text), 0
    return parse_shift_assignment_lines(text), 0


def _find_nmredata_1h_signal_labels(mol) -> set[str]:
    signal_labels = set()
    for property_name in mol.GetPropNames(includePrivate=True, includeComputed=False):
        upper_name = property_name.upper()
        if not re.match(r"^NMREDATA_1D_1H(?:#\d+)?$", upper_name):
            continue
        for raw_line in mol.GetProp(property_name).splitlines():
            clean_line = _clean_nmredata_line(raw_line)
            match = re.search(r"\bL\s*=\s*(s\d+)\b", clean_line, flags=re.IGNORECASE)
            if match:
                signal_labels.add(match.group(1).lower())
    return signal_labels


def _resolve_nmredata_atom_index(atom_number: int, mol) -> int | None:
    if mol is None:
        return atom_number

    one_based_index = atom_number - 1
    output_index = _resolve_1h_assignment_output_atom_index(mol, one_based_index)
    if output_index is not None:
        return output_index

    output_index = _resolve_1h_assignment_output_atom_index(mol, atom_number)
    if output_index is not None:
        return output_index

    return None


def _resolve_1h_assignment_output_atom_index(mol, atom_index: int) -> int | None:
    if atom_index < 0 or atom_index >= mol.GetNumAtoms():
        return None

    atom = mol.GetAtomWithIdx(atom_index)
    if atom.GetAtomicNum() == 1:
        heavy_neighbors = [
            neighbor.GetIdx()
            for neighbor in atom.GetNeighbors()
            if neighbor.GetAtomicNum() != 1
        ]
        if len(heavy_neighbors) == 1:
            return _heavy_atom_output_index(mol, heavy_neighbors[0])
        return None

    if not _mol_has_explicit_hydrogens(mol) and atom.GetTotalNumHs() > 0:
        return atom_index

    return None


def _mol_has_explicit_hydrogens(mol) -> bool:
    return any(atom.GetAtomicNum() == 1 for atom in mol.GetAtoms())


def _heavy_atom_output_index(mol, atom_index: int) -> int:
    heavy_atom_indices = [
        atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() != 1
    ]
    return heavy_atom_indices.index(atom_index)


def _clean_nmredata_line(line: str) -> str:
    return line.strip().rstrip("\\").strip()


def _mol_to_canonical_smiles(mol) -> str:
    try:
        return Chem.MolToSmiles(mol, canonical=True)
    except Exception:
        return ""


def _mol_to_database_smiles(mol) -> str:
    try:
        return Chem.MolToSmiles(Chem.RemoveHs(mol), canonical=True)
    except Exception:
        return Chem.MolToSmiles(mol, canonical=True)


def _preview_sdf_property(value: str, max_length: int = 240) -> str:
    preview = " | ".join(
        _clean_nmredata_line(line) for line in value.splitlines() if line.strip()
    )
    if len(preview) <= max_length:
        return preview
    return f"{preview[: max_length - 3]}..."


def _get_molecule_name(mol) -> str:
    for property_name in MOLECULE_NAME_PROPERTIES:
        if mol.HasProp(property_name):
            value = mol.GetProp(property_name).strip()
            if value:
                return value
    return Chem.MolToSmiles(mol, canonical=True)


def _parse_shift_assignment_line(line: str) -> tuple[int, float, str] | None:
    clean_line = line.strip()
    if not clean_line or clean_line.startswith(("#", "//")):
        return None

    key_value_record = _parse_key_value_assignment(clean_line)
    if key_value_record is not None:
        return key_value_record

    for pattern in [
        r"^\s*(?:atom\s*[=:]\s*)?(?P<atom>\d+)\s*[,;:]\s*"
        r"(?P<shift>-?\d+(?:\.\d+)?)\s*[,;]?\s*(?P<label>.*)$",
        r"^\s*(?:atom\s*[=:]\s*)?(?P<atom>\d+)\s+"
        r"(?P<shift>-?\d+(?:\.\d+)?)\s*(?P<label>.*)$",
        r"^\s*H(?P<atom>\d+)\s*[,;: ]\s*"
        r"(?P<shift>-?\d+(?:\.\d+)?)\s*[,;]?\s*(?P<label>.*)$",
    ]:
        match = re.search(pattern, clean_line, flags=re.IGNORECASE)
        if match:
            return (
                int(match.group("atom")),
                float(match.group("shift")),
                _clean_assignment_label(match.group("label")),
            )

    return None


def _parse_key_value_assignment(line: str) -> tuple[int, float, str] | None:
    atom_match = re.search(r"\batom(?:_index)?\s*[:=]\s*(\d+)\b", line, re.IGNORECASE)
    if not atom_match:
        return None

    shift_match = re.search(
        r"\b(?:shift|delta|ppm)\s*[:=]\s*(-?\d+(?:\.\d+)?)\b",
        line,
        re.IGNORECASE,
    )
    if shift_match is None:
        shift_match = re.search(r"-?\d+(?:\.\d+)?", line)
    if shift_match is None:
        return None

    label_match = re.search(r"\blabel\s*[:=]\s*([^;,\n]+)", line, re.IGNORECASE)
    label = label_match.group(1) if label_match else ""
    return (
        int(atom_match.group(1)),
        float(shift_match.group(1) if shift_match.lastindex else shift_match.group(0)),
        _clean_assignment_label(label),
    )


def _clean_assignment_label(value: str) -> str:
    label = value.strip(" ,;:")
    return re.sub(r"\s+", " ", label)


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
