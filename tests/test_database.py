from pathlib import Path

import pandas as pd
import pytest

from nmr_calculator.database import (
    DATABASE_COLUMNS,
    load_database_cache,
    load_nmrshiftdb2_sdf,
    load_shift_database_csv,
    normalize_shift_database,
    save_database_cache,
)

FIXTURE_PATH = Path("tests/fixtures/test_1h_shift_database.csv")


def test_load_shift_database_csv_loads_fixture():
    database = load_shift_database_csv(str(FIXTURE_PATH))

    assert list(database.columns) == DATABASE_COLUMNS
    assert len(database) == 6
    assert pd.api.types.is_integer_dtype(database["atom_index"])
    assert pd.api.types.is_float_dtype(database["shift_ppm"])


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


def test_load_nmrshiftdb2_sdf_placeholder_raises_not_implemented():
    with pytest.raises(NotImplementedError, match="SDF parsing"):
        load_nmrshiftdb2_sdf("nmrshiftdb2.sdf")


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
