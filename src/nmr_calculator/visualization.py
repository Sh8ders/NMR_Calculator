"""Molecule structure image generation utilities."""

from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

from nmr_calculator.molecule import parse_molecule


def save_molecule_image(
    smiles: str, output_path: str, show_atom_indices: bool = True
) -> str:
    """Save a 2D molecule depiction as a PNG file."""
    image_bytes = get_molecule_image_bytes(
        smiles, show_atom_indices=show_atom_indices
    )
    output_file = Path(output_path)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    output_file.write_bytes(image_bytes)
    return str(output_file)


def get_molecule_image_bytes(
    smiles: str, show_atom_indices: bool = True
) -> bytes:
    """Return PNG bytes for a 2D molecule depiction."""
    mol = Chem.Mol(parse_molecule(smiles))
    rdDepictor.Compute2DCoords(mol)

    drawer = rdMolDraw2D.MolDraw2DCairo(500, 350)
    options = drawer.drawOptions()
    options.addAtomIndices = show_atom_indices

    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    return bytes(drawer.GetDrawingText())
