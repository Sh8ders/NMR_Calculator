"""RDKit molecule parsing and proton environment inspection utilities."""

from rdkit import Chem
from rdkit.Chem.rdchem import Atom, BondType, HybridizationType, Mol

HALOGEN_ATOMIC_NUMBERS = {9, 17, 35, 53}


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
        return _exchangeable_label(atom)

    if symbol != "C":
        return "unknown proton environment"

    if is_aldehyde_proton(mol, atom_index):
        return "aldehyde proton"

    if atom.GetIsAromatic():
        return "aromatic proton"

    if is_terminal_alkyne_proton(mol, atom_index):
        return "terminal alkyne proton"

    if is_vinylic(mol, atom_index):
        return "vinylic/alkene proton"

    if is_alpha_to_carbonyl(mol, atom_index):
        return "alkyl alpha to carbonyl"

    attached_heteroatom_type = get_attached_heteroatom_type(mol, atom_index)
    if attached_heteroatom_type:
        return "heteroatom-adjacent alkyl proton"

    if is_benzylic(mol, atom_index):
        return "benzylic proton"

    if is_allylic(mol, atom_index):
        return "allylic proton"

    if atom.GetHybridization() == HybridizationType.SP3:
        hydrogen_count = _attached_hydrogen_count(atom)
        if hydrogen_count == 3:
            return "simple alkyl CH3"
        if hydrogen_count == 2:
            return "simple alkyl CH2"
        if hydrogen_count == 1:
            return "simple alkyl CH"

    return "unknown proton environment"


def is_exchangeable_proton_environment(environment_label: str) -> bool:
    """Return whether a label describes exchangeable OH, NH, or SH protons."""
    return "exchangeable proton" in environment_label


def is_alpha_to_carbonyl(mol: Mol, atom_index: int) -> bool:
    """Return whether an atom is bonded to a carbonyl carbon."""
    atom = mol.GetAtomWithIdx(atom_index)
    if atom.GetAtomicNum() == 1:
        return False

    return any(_is_carbonyl_carbon(neighbor) for neighbor in atom.GetNeighbors())


def is_benzylic(mol: Mol, atom_index: int) -> bool:
    """Return whether an sp3 carbon is directly attached to an aromatic ring."""
    atom = mol.GetAtomWithIdx(atom_index)
    if atom.GetSymbol() != "C" or atom.GetIsAromatic():
        return False
    if atom.GetHybridization() != HybridizationType.SP3:
        return False

    return any(neighbor.GetIsAromatic() for neighbor in atom.GetNeighbors())


def is_allylic(mol: Mol, atom_index: int) -> bool:
    """Return whether an sp3 carbon is directly attached to an alkene carbon."""
    atom = mol.GetAtomWithIdx(atom_index)
    if atom.GetSymbol() != "C" or atom.GetHybridization() != HybridizationType.SP3:
        return False

    return any(_has_alkene_bond(neighbor) for neighbor in atom.GetNeighbors())


def is_vinylic(mol: Mol, atom_index: int) -> bool:
    """Return whether a proton-bearing carbon is part of a C=C bond."""
    atom = mol.GetAtomWithIdx(atom_index)
    if atom.GetSymbol() != "C":
        return False

    return _has_alkene_bond(atom)


def is_terminal_alkyne_proton(mol: Mol, atom_index: int) -> bool:
    """Return whether a proton-bearing carbon is part of a terminal alkyne."""
    atom = mol.GetAtomWithIdx(atom_index)
    if atom.GetSymbol() != "C":
        return False

    return _has_alkyne_bond(atom) and _attached_hydrogen_count(atom) > 0


def is_aldehyde_proton(mol: Mol, atom_index: int) -> bool:
    """Return whether an atom is an aldehyde carbon bearing a proton."""
    atom = mol.GetAtomWithIdx(atom_index)
    return _is_aldehyde_carbon(atom)


def is_carboxylic_acid_proton(mol: Mol, atom_index: int) -> bool:
    """Return whether an oxygen atom is a protonated carboxylic acid oxygen."""
    atom = mol.GetAtomWithIdx(atom_index)
    if atom.GetSymbol() != "O" or _attached_hydrogen_count(atom) < 1:
        return False

    return any(_is_carbonyl_carbon(neighbor) for neighbor in atom.GetNeighbors())


def get_attached_heteroatom_type(mol: Mol, atom_index: int) -> str | None:
    """Return the type of directly attached heteroatom, if present."""
    atom = mol.GetAtomWithIdx(atom_index)
    for neighbor in atom.GetNeighbors():
        atomic_num = neighbor.GetAtomicNum()
        if atomic_num in HALOGEN_ATOMIC_NUMBERS:
            return "halogen"
        if neighbor.GetSymbol() in {"O", "N", "S"}:
            return neighbor.GetSymbol()
    return None


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


def _exchangeable_label(atom: Atom) -> str:
    if atom.GetSymbol() == "O":
        if is_carboxylic_acid_proton(atom.GetOwningMol(), atom.GetIdx()):
            return "carboxylic acid exchangeable proton"
        return "alcohol exchangeable proton"
    if atom.GetSymbol() == "N":
        return "amine exchangeable proton"
    if atom.GetSymbol() == "S":
        return "thiol exchangeable proton"
    return "unknown proton environment"


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


def _has_alkyne_bond(atom: Atom) -> bool:
    return any(
        neighbor.GetSymbol() == "C"
        and atom.GetOwningMol().GetBondBetweenAtoms(
            atom.GetIdx(), neighbor.GetIdx()
        ).GetBondType()
        == BondType.TRIPLE
        for neighbor in atom.GetNeighbors()
    )


def _is_carbonyl_carbon(atom: Atom) -> bool:
    if atom.GetSymbol() != "C":
        return False

    mol = atom.GetOwningMol()
    return any(
        neighbor.GetSymbol() == "O"
        and mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx()).GetBondType()
        == BondType.DOUBLE
        for neighbor in atom.GetNeighbors()
    )
