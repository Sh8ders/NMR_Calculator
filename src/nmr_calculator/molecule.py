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


def get_proton_environment_groups(mol: Mol) -> list[dict[str, object]]:
    """Group hydrogen-bearing atoms into equivalent proton environments."""
    symmetry_classes = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
    grouped_environments: dict[tuple[object, ...], list[dict[str, object]]] = {}

    for environment in get_hydrogen_bearing_atoms(mol):
        atom_index = int(environment["atom_index"])
        key = build_proton_group_key(environment, symmetry_classes[atom_index])
        grouped_environments.setdefault(key, []).append(environment)

    groups = []
    for group_id, environments in enumerate(grouped_environments.values(), start=1):
        atom_indices = sorted(int(environment["atom_index"]) for environment in environments)
        representative_atom_index = atom_indices[0]
        representative_environment = min(
            environments, key=lambda environment: int(environment["atom_index"])
        )
        environment_label = str(representative_environment["environment_label"])

        groups.append(
            {
                "group_id": group_id,
                "atom_indices": atom_indices,
                "proton_count": sum(
                    int(environment["total_hydrogens"]) for environment in environments
                ),
                "representative_atom_index": representative_atom_index,
                "atom_symbols": sorted(
                    {str(environment["atom_symbol"]) for environment in environments}
                ),
                "environment_label": environment_label,
                "symmetry_class": symmetry_classes[representative_atom_index],
                "is_exchangeable": is_exchangeable_proton_environment(environment_label),
                "is_aromatic": all(
                    bool(environment["aromatic"]) for environment in environments
                ),
                "neighbor_summary": summarize_neighbors(representative_environment),
            }
        )

    return groups


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


def is_exchangeable_proton_environment(environment_label: str) -> bool:
    """Return whether a label describes exchangeable OH, NH, or SH protons."""
    return "exchangeable proton" in environment_label


def build_proton_group_key(
    environment: dict[str, object], symmetry_class: int
) -> tuple[object, ...]:
    """Build a stable grouping key for equivalent proton environments."""
    return (
        symmetry_class,
        environment["atom_symbol"],
        environment["environment_label"],
        environment["aromatic"],
        environment["ring_membership"],
        tuple(environment["neighboring_atom_symbols"]),
        tuple(environment["bond_types_to_neighbors"]),
        environment["total_hydrogens"],
        environment["formal_charge"],
    )


def summarize_neighbors(environment: dict[str, object]) -> str:
    """Summarize neighboring heavy atoms and bond types for display."""
    neighboring_symbols = environment["neighboring_atom_symbols"]
    bond_types = environment["bond_types_to_neighbors"]
    if not neighboring_symbols:
        return "no heavy-atom neighbors"

    pairs = [
        f"{bond_type} to {symbol}"
        for symbol, bond_type in zip(neighboring_symbols, bond_types, strict=True)
    ]
    return ", ".join(pairs)


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
