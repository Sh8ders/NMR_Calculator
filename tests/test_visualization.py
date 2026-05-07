import pytest

from nmr_calculator.visualization import (
    get_molecule_image_bytes,
    save_molecule_image,
)


def test_ethanol_molecule_image_creates_png(tmp_path):
    output_path = tmp_path / "ethanol_structure.png"

    saved_path = save_molecule_image("CCO", str(output_path))

    assert saved_path == str(output_path)
    assert output_path.exists()
    assert output_path.stat().st_size > 0


def test_benzene_molecule_image_creates_png(tmp_path):
    output_path = tmp_path / "benzene_structure.png"

    save_molecule_image("c1ccccc1", str(output_path))

    assert output_path.exists()
    assert output_path.stat().st_size > 0


def test_molecule_image_bytes_returns_png_bytes():
    image_bytes = get_molecule_image_bytes("CCO")

    assert isinstance(image_bytes, bytes)
    assert len(image_bytes) > 0
    assert image_bytes.startswith(b"\x89PNG")


def test_molecule_image_without_atom_indices_creates_png(tmp_path):
    output_path = tmp_path / "ethanol_no_indices.png"

    save_molecule_image("CCO", str(output_path), show_atom_indices=False)

    assert output_path.exists()
    assert output_path.stat().st_size > 0


def test_save_molecule_image_raises_value_error_for_invalid_smiles(tmp_path):
    output_path = tmp_path / "invalid.png"

    with pytest.raises(ValueError, match="Invalid SMILES string"):
        save_molecule_image("not-a-smiles", str(output_path))


def test_molecule_image_bytes_raises_value_error_for_invalid_smiles():
    with pytest.raises(ValueError, match="Invalid SMILES string"):
        get_molecule_image_bytes("not-a-smiles")
