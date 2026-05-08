from pathlib import Path

import pandas as pd
import pytest
from rdkit import Chem

from nmr_calculator.database import (
    DATABASE_COLUMNS,
    add_canonical_smiles_column,
    download_nmrshiftdb2_data,
    extract_1h_shift_records_from_mol,
    filter_database_by_smiles,
    find_1h_nmr_property_names,
    inspect_sdf_properties,
    load_database_cache,
    load_nmrshiftdb2_sdf,
    load_shift_database_csv,
    normalize_shift_database,
    parse_nmredata_signal_assignments,
    parse_shift_assignment_lines,
    request,
    save_database_cache,
)

FIXTURE_PATH = Path("tests/fixtures/test_1h_shift_database.csv")
NMRSHIFTDB2_FIXTURE_PATH = Path("tests/fixtures/tiny_nmrshiftdb2_1h.sdf")


def test_load_shift_database_csv_loads_fixture():
    database = load_shift_database_csv(str(FIXTURE_PATH))

    assert list(database.columns) == DATABASE_COLUMNS
    assert len(database) == 6
    assert pd.api.types.is_integer_dtype(database["atom_index"])
    assert pd.api.types.is_float_dtype(database["shift_ppm"])


def test_add_canonical_smiles_column_adds_expected_column():
    database = load_shift_database_csv(str(FIXTURE_PATH))
    database_with_canonical = add_canonical_smiles_column(database)

    assert "canonical_smiles" in database_with_canonical.columns
    assert database_with_canonical["canonical_smiles"].notna().all()


def test_filter_database_by_smiles_returns_exact_canonical_matches():
    database = load_shift_database_csv(str(FIXTURE_PATH))

    ethanol = filter_database_by_smiles(database, "OCC")
    benzene = filter_database_by_smiles(database, "c1ccccc1")
    missing = filter_database_by_smiles(database, "C=C")

    assert len(ethanol) == 3
    assert len(benzene) == 1
    assert missing.empty


def test_missing_required_column_raises_value_error():
    df = pd.DataFrame(
        {
            "smiles": ["CCO"],
            "molecule_name": ["ethanol"],
            "atom_index": [0],
            "shift_ppm": [1.2],
            "assignment_label": ["CH3"],
        }
    )

    with pytest.raises(ValueError, match="missing required columns"):
        normalize_shift_database(df)


def test_invalid_smiles_raises_value_error():
    df = _valid_database_frame()
    df.loc[0, "smiles"] = "not-a-smiles"

    with pytest.raises(ValueError, match="Invalid SMILES"):
        normalize_shift_database(df)


def test_invalid_shift_ppm_raises_value_error():
    df = _valid_database_frame(shift_ppm="not-numeric")

    with pytest.raises(ValueError, match="shift_ppm"):
        normalize_shift_database(df)


def test_invalid_atom_index_raises_value_error():
    df = _valid_database_frame(atom_index="not-an-index")

    with pytest.raises(ValueError, match="atom_index"):
        normalize_shift_database(df)


def test_database_cache_roundtrip(tmp_path):
    database = load_shift_database_csv(str(FIXTURE_PATH))
    cache_path = tmp_path / "test_cache.csv"

    saved_path = save_database_cache(database, str(cache_path))
    loaded_cache = load_database_cache(saved_path)

    assert saved_path == str(cache_path)
    assert cache_path.exists()
    pd.testing.assert_frame_equal(loaded_cache, database)


def test_find_1h_nmr_property_names_finds_fixture_property():
    mol = _first_nmrshiftdb2_fixture_mol()

    property_names = find_1h_nmr_property_names(mol)

    assert "NMREDATA_1D_1H" in property_names


def test_extract_1h_shift_records_from_mol_extracts_fixture_records():
    mol = _first_nmrshiftdb2_fixture_mol()

    records = extract_1h_shift_records_from_mol(mol)

    assert len(records) == 3
    assert records[0] == {
        "smiles": "CCO",
        "molecule_name": "ethanol",
        "atom_index": 0,
        "shift_ppm": 1.23,
        "assignment_label": "CH3",
        "source": "nmrshiftdb2",
    }


def test_parse_shift_assignment_lines_supports_flexible_formats():
    records = parse_shift_assignment_lines(
        "0,1.23,CH3\n"
        "1: 3.65 CH2\n"
        "2 2.50 OH\n"
        "H1, 1.23\n"
        "1.23; atom=0; label=CH3"
    )

    assert records == [
        (0, 1.23, "CH3"),
        (1, 3.65, "CH2"),
        (2, 2.5, "OH"),
        (1, 1.23, ""),
        (0, 1.23, "CH3"),
    ]


def test_parse_nmredata_signal_assignments_supports_real_signal_format():
    mol = Chem.AddHs(Chem.MolFromSmiles("CO"))
    hydrogen_atom_numbers = [
        atom.GetIdx() + 1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1
    ]
    text = (
        f"s0, 3.25, {hydrogen_atom_numbers[0]}, {hydrogen_atom_numbers[1]}\\\n"
        "s1, 52.0, 1\\\n"
        "s2, 1.10\\\n"
    )

    records, skipped_records = parse_nmredata_signal_assignments(
        text, mol=mol, one_h_signal_labels={"s0", "s2"}
    )

    assert records == [(0, 3.25, "s0")]
    assert skipped_records == 1


def test_extract_1h_shift_records_from_mol_supports_real_nmredata_assignments():
    mol = Chem.AddHs(Chem.MolFromSmiles("CO"))
    mol.SetProp("_Name", "methanol")
    hydrogen_atom_numbers = [
        atom.GetIdx() + 1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1
    ]
    mol.SetProp(
        "NMREDATA_1D_1H",
        "Spectrum_Location=molecule/1\\\n"
        "Larmor=400.0\\\n"
        "3.25, L=s0\\\n"
        "1.10, L=s1\\\n",
    )
    mol.SetProp(
        "NMREDATA_ASSIGNMENT",
        f"s0, 3.25, {hydrogen_atom_numbers[0]}, {hydrogen_atom_numbers[1]}\\\n"
        "s9, 52.0, 1\\\n"
        f"s1, 1.10, {hydrogen_atom_numbers[-1]}\\\n",
    )

    records = extract_1h_shift_records_from_mol(mol)

    assert [record["atom_index"] for record in records] == [0, 1]
    assert [record["shift_ppm"] for record in records] == [3.25, 1.1]
    assert [record["assignment_label"] for record in records] == ["s0", "s1"]


def test_inspect_sdf_properties_reports_tiny_fixture_properties():
    inspected = inspect_sdf_properties(str(NMRSHIFTDB2_FIXTURE_PATH), max_molecules=1)

    assert inspected == [
        {
            "molecule_index": 0,
            "molecule_name": "ethanol",
            "canonical_smiles": "CCO",
            "property_names": [
                "_Name",
                "_MolFileInfo",
                "_MolFileComments",
                "_MolFileChiralFlag",
                "NMREDATA_1D_1H",
            ],
            "property_previews": {
                "_Name": "ethanol",
                "_MolFileInfo": "RDKit          2D",
                "_MolFileComments": "",
                "_MolFileChiralFlag": "0",
                "NMREDATA_1D_1H": "0,1.23,CH3 | 1: 3.65 CH2 | 2 2.50 OH",
            },
        }
    ]


def test_load_nmrshiftdb2_sdf_loads_tiny_fixture():
    database = load_nmrshiftdb2_sdf(str(NMRSHIFTDB2_FIXTURE_PATH))

    assert list(database.columns) == DATABASE_COLUMNS
    assert len(database) == 5
    assert set(database["source"]) == {"nmrshiftdb2"}
    pd.testing.assert_frame_equal(normalize_shift_database(database), database)


def test_download_nmrshiftdb2_data_supports_monkeypatched_urlopen(
    monkeypatch, tmp_path
):
    class FakeResponse:
        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc_value, traceback):
            return False

        def read(self):
            return b"tiny sdf"

    def fake_urlopen(url):
        assert url == "https://example.test/nmr.sd"
        return FakeResponse()

    monkeypatch.setattr(request, "urlopen", fake_urlopen)
    output_path = tmp_path / "raw" / "nmr.sd"

    saved_path = download_nmrshiftdb2_data(
        str(output_path), url="https://example.test/nmr.sd"
    )

    assert saved_path == str(output_path)
    assert output_path.read_bytes() == b"tiny sdf"


def _valid_database_frame(atom_index=0, shift_ppm=1.2):
    return pd.DataFrame(
        {
            "smiles": ["CCO"],
            "molecule_name": ["ethanol"],
            "atom_index": [atom_index],
            "shift_ppm": [shift_ppm],
            "assignment_label": ["CH3"],
            "source": ["test"],
        }
    )


def _first_nmrshiftdb2_fixture_mol():
    supplier = Chem.SDMolSupplier(str(NMRSHIFTDB2_FIXTURE_PATH), removeHs=False)
    return supplier[0]
