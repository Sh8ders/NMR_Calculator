"""RDKit molecule parsing and proton environment inspection utilities."""

from rdkit import Chem
from rdkit.Chem.rdchem import Atom, BondType, HybridizationType, Mol


def parse_molecule(smiles: str) -> Mol:
    """Parse a SMILES string into an RDKit molecule."""
    normalized_smiles = smiles.strip()
    if not normalized_smiles:
        raise ValueError("Invalid SMILES string: input is empty")

    mol = Chem.MolFromSmiles(normalized_smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: {smiles!r}")
    return mol


def prepare_molecule(smiles: str) -> Mol:
    """Parse a SMILES string and add explicit hydrogens for proton inspection."""
    mol = parse_molecule(smiles)
    return Chem.AddHs(mol)


def get_hydrogen_bearing_atoms(mol: Mol) -> list[dict[str, object]]:
    """Return descriptors for non-hydrogen atoms with at least one attached hydrogen."""
    hydrogen_bearing_atoms = []

    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1:
            continue

        total_hydrogens = _attached_hydrogen_count(atom)
        if total_hydrogens < 1:
            continue

        atom_index = atom.GetIdx()
        heavy_neighbors = [
            neighbor for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() != 1
        ]
        hydrogen_bearing_atoms.append(
            {
                "atom_index": atom_index,
                "atom_symbol": atom.GetSymbol(),
                "total_hydrogens": total_hydrogens,
                "hybridization": str(atom.GetHybridization()),
                "aromatic": atom.GetIsAromatic(),
                "ring_membership": atom.IsInRing(),
                "formal_charge": atom.GetFormalCharge(),
                "neighboring_atom_symbols": [
                    neighbor.GetSymbol() for neighbor in heavy_neighbors
                ],
                "bond_types_to_neighbors": [
                    str(mol.GetBondBetweenAtoms(atom_index, neighbor.GetIdx()).GetBondType())
                    for neighbor in heavy_neighbors
                ],
                "environment_label": label_proton_environment(mol, atom_index),
            }
        )

    return hydrogen_bearing_atoms


def label_proton_environment(mol: Mol, atom_index: int) -> str:
    """Assign a simple descriptive label to hydrogens attached to an atom."""
    atom = mol.GetAtomWithIdx(atom_index)
    if atom.GetAtomicNum() == 1 or _attached_hydrogen_count(atom) < 1:
        return "unknown proton environment"

    symbol = atom.GetSymbol()
    if symbol in {"O", "N", "S"}:
        return _exchangeable_label(symbol)

    if symbol != "C":
        return "unknown proton environment"

    if _is_aldehyde_carbon(atom):
        return "aldehyde proton"

    if atom.GetIsAromatic():
        return "aromatic proton"

    if _has_alkene_bond(atom):
        return "alkene proton"

    if _is_heteroatom_adjacent(atom):
        return "heteroatom-adjacent alkyl proton"

    if atom.GetHybridization() == HybridizationType.SP3:
        hydrogen_count = _attached_hydrogen_count(atom)
        if hydrogen_count == 3:
            return "alkyl CH3"
        if hydrogen_count == 2:
            return "alkyl CH2"
        if hydrogen_count == 1:
            return "alkyl CH"

    return "unknown proton environment"


def _attached_hydrogen_count(atom: Atom) -> int:
    explicit_hydrogens = sum(
        1 for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() == 1
    )
    return explicit_hydrogens + atom.GetNumImplicitHs()


def _exchangeable_label(symbol: str) -> str:
    labels = {
        "O": "alcohol/amine/thiol exchangeable proton",
        "N": "alcohol/amine/thiol exchangeable proton",
        "S": "alcohol/amine/thiol exchangeable proton",
    }
    return labels.get(symbol, "unknown proton environment")


def _is_aldehyde_carbon(atom: Atom) -> bool:
    return any(
        neighbor.GetSymbol() == "O"
        and atom.GetOwningMol().GetBondBetweenAtoms(
            atom.GetIdx(), neighbor.GetIdx()
        ).GetBondType()
        == BondType.DOUBLE
        for neighbor in atom.GetNeighbors()
    )


def _has_alkene_bond(atom: Atom) -> bool:
    return any(
        neighbor.GetSymbol() == "C"
        and atom.GetOwningMol().GetBondBetweenAtoms(
            atom.GetIdx(), neighbor.GetIdx()
        ).GetBondType()
        == BondType.DOUBLE
        for neighbor in atom.GetNeighbors()
    )


def _is_heteroatom_adjacent(atom: Atom) -> bool:
    return any(
        neighbor.GetAtomicNum() not in {1, 6}
        for neighbor in atom.GetNeighbors()
    )
